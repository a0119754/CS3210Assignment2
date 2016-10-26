#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>

/***********************************************************
  Helper functions 
***********************************************************/

//For exiting on error condition
void die(int lineNo);

//For trackinng execution
long long wallClockTime();


/***********************************************************
  Square matrix related functions, used by both world and pattern
***********************************************************/

char** allocateSquareMatrix( int size, char defaultValue );

void freeSquareMatrix( char** );

void printSquareMatrix( char**, int size );


/***********************************************************
   World  related functions
***********************************************************/

#define ALIVE 'X' 
#define DEAD 'O'

char** readWorldFromFile( char* fname, int* size );

int countNeighbours(char** world, int row, int col);

void evolveWorld(char** curWorld, char** nextWorld, int size);


/***********************************************************
   Simple circular linked list for match records
***********************************************************/

typedef struct MSTRUCT {
    int iteration, row, col, rotation;
    struct MSTRUCT *next;
} MATCH;


typedef struct {
    int nItem;
    MATCH* tail;
} MATCHLIST;

MATCHLIST* newList();

void deleteList( MATCHLIST*);

void insertEnd(MATCHLIST*, int, int, int, int);

void printList(MATCHLIST*);


/***********************************************************
   Search related functions
***********************************************************/

//Using the compass direction to indicate the rotation of pattern
#define N 0 //no rotation
#define E 1 //90 degree clockwise
#define S 2 //180 degree clockwise
#define W 3 //90 degree anti-clockwise

char** readPatternFromFile( char* fname, int* size );

void rotate90(char** current, char** rotated, int size);

void searchPatterns(char** world, int wSize, int iteration, 
        char** patterns[4], int pSize, MATCHLIST* list);

void searchPatternInOnlyOneDirection(char** world, int wSize, int iteration, 
        char** patterns[4], int pSize, MATCHLIST* list, int dir);

void searchSinglePattern(char** world, int wSize, int interation,
        char** pattern, int pSize, int rotation, MATCHLIST* list);

/***********************************************************
   Main function
***********************************************************/


int main( int argc, char** argv)
{
    char **curW, **nextW, **temp, dummy[20];
    char **patterns[4];
    int dir, iterations, iter, processes, rank, i, j, noOfResults;
    int size, patternSize;
	int forResults[4];
    long long before, after;
    MATCHLIST*list;
	MATCH* cur;
	MPI_Status mpiStatus;
    
    if (argc < 4 ){
		if (rank == 0)
			fprintf(stderr, "Usage: %s <world file> <Iterations> <pattern file>\n", argv[0]);
        exit(1);
    } 

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	if (processes != 5) {
		fprintf(stderr, "Currently this parallelized version of the code is written for the specific usage of five processes");
		exit(1);
	}
	
	if (rank == 0) {
		// Master will read the world from file
		curW = readWorldFromFile(argv[1], &size);
		nextW = allocateSquareMatrix(size+2, DEAD);
	}

    iterations = atoi(argv[2]);
	
	// Let only master printf to prevent clutter
	if (rank == 0) {
		printf("World Size = %d\n", size);
		printf("Iterations = %d\n", iterations);
	}

	// Master will read the pattern from file and rotate them
	if (rank == 0) {
		patterns[N] = readPatternFromFile(argv[3], &patternSize);
		for (dir = E; dir <= W; dir++){
			patterns[dir] = allocateSquareMatrix(patternSize, DEAD);
			rotate90(patterns[dir-1], patterns[dir], patternSize);
		}
		printf("Pattern size = %d\n", patternSize);
		
#ifdef DEBUG
    printSquareMatrix(patterns[N], patternSize);
    printSquareMatrix(patterns[E], patternSize);
    printSquareMatrix(patterns[S], patternSize);
    printSquareMatrix(patterns[W], patternSize);
#endif
	}
 
    //Start timer
	if (rank == 0)
		before = wallClockTime();
	
	// Master pass slaves size of world
	if (rank == 0) {
		for (i = 1; i < processes; i++) {
			MPI_Send(&size, 1, MPI_INT, i, 1, MPI_COMM_WORLD); // This message has tag of 1
		}
	} else {
		MPI_Recv(&size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &mpiStatus);
		curW = allocateSquareMatrix(size+2, DEAD); // Slaves allocate memory of world
	}
	
	// Master pass slaves size of pattern
	if (rank == 0) {
		for (i = 1; i < processes; i++) {
			MPI_Send(&patternSize, 1, MPI_INT, i, 2, MPI_COMM_WORLD); // This message has tag of 2
		}
	} else {
		MPI_Recv(&patternSize, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &mpiStatus);
		patterns[0] = allocateSquareMatrix(patternSize, DEAD); // Slaves allocate memory for pattern
		// xxx how many patterns should each slave have?
	}
	
	// Distribute patterns to slaves
	if (rank == 0) {
		dir = N;
		for (i = 1; i < processes; i++, dir++) {
			MPI_Send(&(patterns[dir][0][0]), patternSize * patternSize, MPI_CHAR, i, 3, MPI_COMM_WORLD); // This message has tag of 3
		}
	} else {
		MPI_Recv(&(patterns[0][0][0]), patternSize * patternSize, MPI_CHAR, 0, 3, MPI_COMM_WORLD, &mpiStatus);
	}

    //Actual work start
    list = newList();

    for (iter = 0; iter < iterations; iter++){

#ifdef DEBUG
        printf("World Iteration.%d\n", iter);
        printSquareMatrix(curW, size+2);
#endif

		// Distribute world to slaves
		if (rank == 0) {
			for (i = 1; i < processes; i++) {
				MPI_Send(&(curW[0][0]), size * size, MPI_CHAR, i, 4 + iter, MPI_COMM_WORLD); // This message has tag of 4+iter
			}
		} else {
			MPI_Recv(&(curW[0][0]), size * size, MPI_CHAR, 0, 4 + iter, MPI_COMM_WORLD, &mpiStatus);
		}
		
		// Master focus on evolving to next generation
		if (rank == 0) {
			evolveWorld( curW, nextW, size );
		}
		// While slaves search for pattern 
		else {
			searchPatternInOnlyOneDirection(curW, size, iter, patterns, patternSize, list, rank - 1);
		}

		// Only master needs to do this; the slaves only use one world array
		if (rank == 0) {
			// Swap curW and nextW variables once done to move on to the next iteration, and use
			// the other array for world information
			temp = curW;
			curW = nextW;
			nextW = temp;
		}
    }

	// slaves return results to master
	if (rank == 0) {
		for (i = 1; i < processes; i++) {
			MPI_Recv(&noOfResults, 1, MPI_INT, i, 4 + iter, MPI_COMM_WORLD, &mpiStatus);
			
			for (j = 0; j < noOfResults; j++) {
				MPI_Recv(&forResults[0], 1, MPI_INT, i, 5 + iter + j, MPI_COMM_WORLD, &mpiStatus);
				MPI_Recv(&forResults[1], 1, MPI_INT, i, 5 + iter + j, MPI_COMM_WORLD, &mpiStatus);
				MPI_Recv(&forResults[2], 1, MPI_INT, i, 5 + iter + j, MPI_COMM_WORLD, &mpiStatus);
				MPI_Recv(&forResults[3], 1, MPI_INT, i, 5 + iter + j, MPI_COMM_WORLD, &mpiStatus);
				// Again, four consecutive messages with same message tag; hope it is alright.
				
				insertEnd(list, forResults[0], forResults[1], forResults[2], forResults[3]);
			}
		}
		
		// Output results
		printList( list );
		
		//Stop timer
		after = wallClockTime();

		printf("Sequential SETL took %1.2f seconds\n", ((float)(after - before))/1000000000);
	} else {
		noOfResults = list->nItem;
		MPI_Send(&noOfResults, 1, MPI_INT, 0, 4 + iter, MPI_COMM_WORLD);
		cur = list->tail->next; // Personally I have no idea what this is for, but it's in the implementation of printList(), and I don't know exactly how the MATCHLIST structure works exactly and technically, and I'm way too used to object-oriented languages to even wrap my head around the primitive concept of using structs before classes ever existed anyway, so yeah, I am just going to leave it here and pray hard it doesn't do something wrong by skipping the tail
		for (i = 0; i < noOfResults; i++, cur = cur->next) { // Wow didn't know you could do that in a for-loop. Commas, huh? I'll keep that in mind. Guess you DO learn something new everyday
			MPI_Send(&(cur->iteration), 1, MPI_INT, 0, 5 + iter + i, MPI_COMM_WORLD);
			MPI_Send(&(cur->row), 1, MPI_INT, 0, 5 + iter + i, MPI_COMM_WORLD);
			MPI_Send(&(cur->col), 1, MPI_INT, 0, 5 + iter + i, MPI_COMM_WORLD);
			MPI_Send(&(cur->rotation), 1, MPI_INT, 0, 5 + iter + i, MPI_COMM_WORLD);
			// Sending four consecutive messages one after another with the same message tag. Would that be alright?
		}
	}

    //Clean up
    deleteList( list );

    freeSquareMatrix( curW );
	if (rank == 0) // Only master has two world arrays
		freeSquareMatrix( nextW );
	
	// xxx not all slaves are going to create four patterns
    freeSquareMatrix( patterns[0] );
	if (rank == 0) {
		freeSquareMatrix( patterns[1] );
		freeSquareMatrix( patterns[2] );
		freeSquareMatrix( patterns[3] );
	}
	
	MPI_Finalize();

    return 0;
}

/***********************************************************
  Helper functions 
***********************************************************/


void die(int lineNo)
{
    fprintf(stderr, "Error at line %d. Exiting\n", lineNo);
    exit(1);
}

long long wallClockTime( )
{
#ifdef __linux__
    struct timespec tp;
    clock_gettime(CLOCK_REALTIME, &tp);
    return (long long)(tp.tv_nsec + (long long)tp.tv_sec * 1000000000ll);
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (long long)(tv.tv_usec * 1000 + (long long)tv.tv_sec * 1000000000ll);
#endif
}

/***********************************************************
  Square matrix related functions, used by both world and pattern
***********************************************************/

char** allocateSquareMatrix( int size, char defaultValue )
{

    char* contiguous;
    char** matrix;
    int i;

    //Using a least compiler version dependent approach here
    //C99, C11 have a nicer syntax.    
    contiguous = (char*) malloc(sizeof(char) * size * size);
    if (contiguous == NULL) 
        die(__LINE__);


    memset(contiguous, defaultValue, size * size );

    //Point the row array to the right place
    matrix = (char**) malloc(sizeof(char*) * size );
    if (matrix == NULL) 
        die(__LINE__);

    matrix[0] = contiguous;
    for (i = 1; i < size; i++){
        matrix[i] = &contiguous[i*size];
    }

    return matrix;
}

void printSquareMatrix( char** matrix, int size )
{
    int i,j;
    
    for (i = 0; i < size; i++){
        for (j = 0; j < size; j++){
            printf("%c", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void freeSquareMatrix( char** matrix )
{
    if (matrix == NULL) return;

    free( matrix[0] );
}

/***********************************************************
   World  related functions
***********************************************************/

char** readWorldFromFile( char* fname, int* sizePtr )
{
    FILE* inf;
    
    char temp, **world;
    int i, j;
    int size;

    inf = fopen(fname,"r");
    if (inf == NULL)
        die(__LINE__);


    fscanf(inf, "%d", &size);
    fscanf(inf, "%c", &temp);
    
    //Using the "halo" approach
    // allocated additional top + bottom rows
    // and leftmost and rightmost rows to form a boundary
    // to simplify computation of cell along edges
    world = allocateSquareMatrix( size + 2, DEAD );

    for (i = 1; i <= size; i++){
        for (j = 1; j <= size; j++){
            fscanf(inf, "%c", &world[i][j]);
        }
        fscanf(inf, "%c", &temp);
    }

    *sizePtr = size;    //return size
    return world;
    
}

int countNeighbours(char** world, int row, int col)
//Assume 1 <= row, col <= size, no check 
{
    int i, j, count;

    count = 0;
    for(i = row-1; i <= row+1; i++){
        for(j = col-1; j <= col+1; j++){
            count += (world[i][j] == ALIVE );
        }
    }

    //discount the center
    count -= (world[row][col] == ALIVE);

    return count;

}

void evolveWorld(char** curWorld, char** nextWorld, int size)
{
    int i, j, liveNeighbours;

    for (i = 1; i <= size; i++){
        for (j = 1; j <= size; j++){
            liveNeighbours = countNeighbours(curWorld, i, j);
            nextWorld[i][j] = DEAD;

            //Only take care of alive cases
            if (curWorld[i][j] == ALIVE) {

                if (liveNeighbours == 2 || liveNeighbours == 3)
                    nextWorld[i][j] = ALIVE;

            } else if (liveNeighbours == 3)
                    nextWorld[i][j] = ALIVE;
        } 
    }
}

/***********************************************************
   Search related functions
***********************************************************/

char** readPatternFromFile( char* fname, int* sizePtr )
{
    FILE* inf;
    
    char temp, **pattern;
    int i, j;
    int size;

    inf = fopen(fname,"r");
    if (inf == NULL)
        die(__LINE__);


    fscanf(inf, "%d", &size);
    fscanf(inf, "%c", &temp);
    
    pattern = allocateSquareMatrix( size, DEAD );

    for (i = 0; i < size; i++){
        for (j = 0; j < size; j++){
            fscanf(inf, "%c", &pattern[i][j]);
        }
        fscanf(inf, "%c", &temp);
    }
    
    *sizePtr = size;    //return size
    return pattern;
}


void rotate90(char** current, char** rotated, int size)
{
    int i, j;

    for (i = 0; i < size; i++){
        for (j = 0; j < size; j++){
            rotated[j][size-i-1] = current[i][j];
        }
    }
}

void searchPatterns(char** world, int wSize, int iteration, 
        char** patterns[4], int pSize, MATCHLIST* list)
{
    int dir;

    for (dir = N; dir <= W; dir++){
        searchSinglePattern(world, wSize, iteration, 
                patterns[dir], pSize, dir, list);
    }
}

void searchPatternInOnlyOneDirection(char** world, int wSize, int iteration, 
        char** patterns[4], int pSize, MATCHLIST* list, int dir) {
        searchSinglePattern(world, wSize, iteration, 
                patterns[0], pSize, dir, list);
}

void searchSinglePattern(char** world, int wSize, int iteration,
        char** pattern, int pSize, int rotation, MATCHLIST* list)
{
    int wRow, wCol, pRow, pCol, match;


    for (wRow = 1; wRow <= (wSize-pSize+1); wRow++){
        for (wCol = 1; wCol <= (wSize-pSize+1); wCol++){
            match = 1;
#ifdef DEBUGMORE
            printf("S:(%d, %d)\n", wRow-1, wCol-1);
#endif
            for (pRow = 0; match && pRow < pSize; pRow++){
                for (pCol = 0; match && pCol < pSize; pCol++){
                    if(world[wRow+pRow][wCol+pCol] != pattern[pRow][pCol]){
#ifdef DEBUGMORE
                        printf("\tF:(%d, %d) %c != %c\n", pRow, pCol,
                            world[wRow+pRow][wCol+pCol], pattern[pRow][pCol]);
#endif
                        match = 0;    
                    }
                }
            }
            if (match){
                insertEnd(list, iteration, wRow-1, wCol-1, rotation);
#ifdef DEBUGMORE
printf("*** Row = %d, Col = %d\n", wRow-1, wCol-1);
#endif
            }
        }
    }
}

/***********************************************************
   Simple circular linked list for match records
***********************************************************/

MATCHLIST* newList()
{
    MATCHLIST* list;

    list = (MATCHLIST*) malloc(sizeof(MATCHLIST));
    if (list == NULL)
        die(__LINE__);

    list->nItem = 0;
    list->tail = NULL;

    return list;
}

void deleteList( MATCHLIST* list)
{
    MATCH *cur, *next;
    int i;
    //delete items first

    if (list->nItem != 0 ){
        cur = list->tail->next;
        next = cur->next;
        for( i = 0; i < list->nItem; i++, cur = next, next = next->next ) {
            free(cur); 
        }

    }
    free( list );
}

void insertEnd(MATCHLIST* list, 
        int iteration, int row, int col, int rotation)
{
    MATCH* newItem;

    newItem = (MATCH*) malloc(sizeof(MATCH));
    if (newItem == NULL)
        die(__LINE__);

    newItem->iteration = iteration;
    newItem->row = row;
    newItem->col = col;
    newItem->rotation = rotation;

    if (list->nItem == 0){
        newItem->next = newItem;
        list->tail = newItem;
    } else {
        newItem->next = list->tail->next;
        list->tail->next = newItem;
        list->tail = newItem;
    }

    (list->nItem)++;

}

void printList(MATCHLIST* list)
{
    int i;
    MATCH* cur;

    printf("List size = %d\n", list->nItem);    


    if (list->nItem == 0) return;

    cur = list->tail->next;
    for( i = 0; i < list->nItem; i++, cur=cur->next){
        printf("%d:%d:%d:%d\n", 
                cur->iteration, cur->row, cur->col, cur->rotation);
    }
}

