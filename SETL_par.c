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

void printListSorted(MATCHLIST*);


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

void searchSinglePattern(char** world, int wSize, int iteration, char** pattern, int pSize,
		int rotation, MATCHLIST* list, MATCHLIST* listToContinueFinding, int start, int end);
		
void continueSearch(char** world, int wSize, int iteration, char** pattern, int pSize, int rotation,
		MATCHLIST* list, MATCHLIST* listToResumeFrom, MATCHLIST* listToContinueFinding, int start, int end);

/***********************************************************
   Main function
***********************************************************/


int main( int argc, char** argv)
{
	char **curW, **nextW, **temp, dummy[20];
	char **patterns[4];
	int dir, iterations, iter, processes, rank, noOfPatterns, rotation, i, j, k, l, m, n, m2, n2, omg, start, partSize, end, noOfResults;
	int size, patternSize;
	int debug = 1;
	int forResults[4];
	long long before, after;
	MATCHLIST* list;
	MATCHLIST* prevList;
	MATCHLIST* nextList;
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
		if (debug) printf("Passing world size to slaves\n");
		for (i = 1; i < processes; i++) {
			MPI_Send(&size, 1, MPI_INT, i, 1, MPI_COMM_WORLD); // This message has tag of 1
		}
	} else {
		MPI_Recv(&size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &mpiStatus);
		curW = allocateSquareMatrix(size+2, DEAD); // Slaves allocate memory of world
		// xxx if processes > 5, some slaves would not need the entire world
		// However, simply allocating memory enough for entire world for now, since
		// 1. memory usage is rarely our concern, 2. easier to code since their indexes are of the actual ones
	}
	
	// Master pass slaves size of pattern
	if (rank == 0) {
		if (debug) printf("Passing pattern size to slaves\n");
		for (i = 1; i < processes; i++) {
			MPI_Send(&patternSize, 1, MPI_INT, i, 2, MPI_COMM_WORLD); // This message has tag of 2
		}
	} else {
		MPI_Recv(&patternSize, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &mpiStatus);
		
		// --- Slaves allocate memory for pattern ---
		// There are at least four slave processes, and thus we can divide the slaves into
		// groups of four, with each group in charge of a single rotated pattern
		// 2 processes = 4 rotated patterns
		// 3 processes = 2 + 2 rotated patterns
		// 4 processes = 2 + 1 + 1 rotated patterns
		j = (processes >= 5) ? 1 : ((processes == 2) ? 4 : ((rank == 1) ? 2 : ((processes == 3) ? 2 : 1)));
		for (i = 0; i < j; i++) {
			patterns[i] = allocateSquareMatrix(patternSize, DEAD);
		}
	}
	
	// Distribute patterns to slaves
	if (rank == 0) {
		noOfPatterns = 4;
		if (debug) printf("Distributing rotated patterns to slaves\n");
		
		if (processes == 2) {
			// Send all four rotated patterns to the one slave
			for (dir = N; dir <= W; dir++){
				MPI_Send(&(patterns[dir][0][0]), patternSize * patternSize, MPI_CHAR, 1, 3, MPI_COMM_WORLD); // This message has tag of 3
			}
		} else if (processes >= 5) {
			// Send one pattern to processes 1,5,9,..., one to 2,6,10,..., one to 3,7,11,..., one to 4,8,12,...
			for (i = 1, dir = N; i < processes; i++, dir++) {
				MPI_Send(&(patterns[dir % 4][0][0]), patternSize * patternSize, MPI_CHAR, i, 3, MPI_COMM_WORLD); // This message has tag of 3
			}
		} else { // 3 or 4
			// Send two patterns to process 1, and either other two patterns to process 2, or split them between processes 2 and 3
			for (dir = N; dir <= E; dir++){
				MPI_Send(&(patterns[dir][0][0]), patternSize * patternSize, MPI_CHAR, 1, 3, MPI_COMM_WORLD); // This message has tag of 3
			}
			if (processes == 3) {				
				for (dir = S; dir <= W; dir++){
					MPI_Send(&(patterns[dir][0][0]), patternSize * patternSize, MPI_CHAR, 2, 3, MPI_COMM_WORLD); // This message has tag of 3
				}
			} else {
				MPI_Send(&(patterns[S][0][0]), patternSize * patternSize, MPI_CHAR, 2, 3, MPI_COMM_WORLD); // This message has tag of 3
				MPI_Send(&(patterns[W][0][0]), patternSize * patternSize, MPI_CHAR, 3, 3, MPI_COMM_WORLD); // This message has tag of 3
			}
		}
	} else {
		if (processes == 2) {
			// Only slave process; accept all four rotated patterns
			noOfPatterns = 4;
		} else if ((rank == 1) && ((processes == 3) || (processes == 4))) {
			// Process 1 when there are either 3 or 4 processes needs to receive two patterns
			noOfPatterns = 2;
		} else {
			noOfPatterns = 1;
		}
		
		for (i = 0; i < noOfPatterns; i++) {
			MPI_Recv(&(patterns[i][0][0]), patternSize * patternSize, MPI_CHAR, 0, 3, MPI_COMM_WORLD, &mpiStatus);
		}
		
		if (debug) {
			printf("Debug: Printing out patterns after receiving them from master\n");
			printSquareMatrix(patterns[N], patternSize);
			printSquareMatrix(patterns[E], patternSize);
			printSquareMatrix(patterns[S], patternSize);
			printSquareMatrix(patterns[W], patternSize);
			printf("Debug: Printing out patterns after receiving them from master\n");
		}
	}
	
	if ((debug) && (rank == 0)) printf("Executing calculations\n");
	
	j = (processes + 2) / 4; // Maximum number of processes to handle one rotated pattern
	if (debug) printf("Rank %d: j = %d\n", rank, j);
	k = (processes - 1) % 4; // Number of rotated patterns handled by maximum number of processes
	if (debug) printf("Rank %d: k = %d\n", rank, k);
	m = size % j; // Remainder trying to split up the world size
	if (debug) printf("Rank %d: m = %d\n", rank, m);
	n = size / j; // How big should each divided world piece be
	if (debug) printf("Rank %d: n = %d\n", rank, n);
	m2 = (j > 1) ? (size % (j - 1)) : 0; // Remainder trying to split up the world size among (j - 1) processes
	if (debug) printf("Rank %d: m2 = %d\n", rank, m2);
	n2 = (j > 1) ? (size / (j - 1)) : 0; // How much to divide the world up into for (j - 1) processes
	if (debug) printf("Rank %d: n2 = %d\n", rank, n2);
	l = (rank == 0) ? 0 : (rank - 2) / 4;
	if (debug) printf("Rank %d: l = %d\n", rank, l);
	start = (rank == 0) ? 0 : (((k == 0) || (((rank - 1) % 4) <= k)) ? ((l <= m) ? (l * (n + 1)) : (l * n + m)) : ((l <= m2) ? (l * (n2 + 1)) : (l * n2 + m2)));
	if (debug) printf("Rank %d: start = %d\n", rank, start);
	partSize = (rank == 0) ? 0 : (((k == 0) || (((rank - 1) % 4) <= k)) ? (n + ((l < m) ? 1 : 0)) : (n2 + ((l < m2) ? 1 : 0)));
	if (debug) printf("Rank %d: partSize = %d\n", rank, partSize);
	end = (rank == 0) ? 0 : (start + partSize - 1);
	if (debug) printf("Rank %d: end = %d\n", rank, end);
	
	if ((debug) && (rank == 0)) printf("Starting work\n");
	
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
				l = (i - 2) / 4; // Which part of the world this process will receive
				start = (rank == 0) ? 0 : (((k == 0) || (((i - 1) % 4) <= k)) ? ((l <= m) ? (l * (n + 1)) : (l * n + m)) : ((l <= m2) ? (l * (n2 + 1)) : (l * n2 + m2)));
				partSize = (rank == 0) ? 0 : (((k == 0) || (((i - 1) % 4) <= k)) ? (n + ((l < m) ? 1 : 0)) : (n2 + ((l < m2) ? 1 : 0)));
				MPI_Send(&(curW[start][0]), size * partSize, MPI_CHAR, i, 4 + iter, MPI_COMM_WORLD);
			}
		} else {
			MPI_Recv(&(curW[start][0]), size * partSize, MPI_CHAR, 0, 4 + iter, MPI_COMM_WORLD, &mpiStatus);
			
			if (debug) {
				printf("Debug: Printing out world after receiving on iteration %d\n", iter);
				printSquareMatrix(curW, size+2);
				printf("Debug: Printing out world after receiving on iteration %d\n", iter);
			}
		}
		
		// Master focus on evolving to next generation
		if (rank == 0) {
			evolveWorld( curW, nextW, size );
		}
		// While slaves search for pattern 
		else {
			if (processes == 2) {
				// One process search patterns in all four directions
				for (i = 0; i < noOfPatterns; i++) {
					searchSinglePattern(curW, size, iter, patterns[i], patternSize, i, list, 0, 1, size - 1);
				}
			} else if (processes >= 5) {							
				prevList = newList();
				nextList = newList();
				
				rotation = (rank - 2) % 4;
				// One process search for pattern in a specific direction according to its rank
				searchSinglePattern(curW, size, iter, patterns[0], patternSize, rotation, list, nextList, start, end);
				
				if (rank <= 4) {
					// Receive results from previous guy
					MPI_Recv(&omg, 1, MPI_INT, rank - 4, 0, MPI_COMM_WORLD, &mpiStatus);
					
					for (i = 0; i < omg; i++) {
						MPI_Recv(&forResults[0], 1, MPI_INT, rank - 4, 0, MPI_COMM_WORLD, &mpiStatus);
						MPI_Recv(&forResults[1], 1, MPI_INT, rank - 4, 0, MPI_COMM_WORLD, &mpiStatus);
						MPI_Recv(&forResults[2], 1, MPI_INT, rank - 4, 0, MPI_COMM_WORLD, &mpiStatus);
						MPI_Recv(&forResults[3], 1, MPI_INT, rank - 4, 0, MPI_COMM_WORLD, &mpiStatus);
						
						insertEnd(prevList, forResults[0], forResults[1], forResults[2], forResults[3]);
					}
					
					continueSearch(curW, size, iter, patterns[0], patternSize, rotation, list, prevList, nextList, start, end);
				}
				
				if ((rank + 4) < processes) {
					// Send results to next guy
					omg = nextList->nItem;
					MPI_Send(&omg, 1, MPI_INT, rank + 4, 0, MPI_COMM_WORLD);
					cur = list->tail->next;
					for (i = 0; i < omg; i++, cur = cur->next) {
						MPI_Send(&(cur->iteration), 1, MPI_INT, 0, 5 + iter + i, MPI_COMM_WORLD);
						MPI_Send(&(cur->row), 1, MPI_INT, 0, 5 + iter + i, MPI_COMM_WORLD);
						MPI_Send(&(cur->col), 1, MPI_INT, 0, 5 + iter + i, MPI_COMM_WORLD);
						MPI_Send(&(cur->rotation), 1, MPI_INT, 0, 5 + iter + i, MPI_COMM_WORLD);
					}
				}
				
				// Clear prevList and nextList at the end of each iteration
				deleteList( prevList );
				deleteList( nextList );
				
			} else { // 3 or 4 processes
				if (rank == 2) { // Search for patterns in first two directions
					for (i = 0; i < noOfPatterns; i++) {
						searchSinglePattern(curW, size, iter, patterns[i], patternSize, i, list, 0, 1, size - 1);
					}
				} else if (processes == 3) { // Search for patterns in last two directions
					for (i = 0; i < noOfPatterns; i++) {
						searchSinglePattern(curW, size, iter, patterns[i], patternSize, i + 2, list, 0, 1, size - 1);
					}
				} else { // 
					searchSinglePattern(curW, size, iter, patterns[0], patternSize, rank, list, 0, 1, size - 1);
				}
			}
			
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
		printListSorted( list );
		
		//Stop timer
		after = wallClockTime();

		printf("Parallel SETL took %1.2f seconds\n", ((float)(after - before))/1000000000);
	} else {
		if (debug) {
			printf("Debug: Slave process %d is printing list\n", rank);
			printList(list);
			printf("Debug: Slave process %d has finished printing list\n", rank);
		}
		
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
	
	// Free memory used for patterns
	for (i = 0; i < noOfPatterns; i++) {
		freeSquareMatrix( patterns[i] );
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
                patterns[dir], pSize, dir, list, 0, 1, wSize - 1);
    }
}

void searchSinglePattern(char** world, int wSize, int iteration, char** pattern, int pSize,
		int rotation, MATCHLIST* list, MATCHLIST* listToContinueFinding, int start, int end)
{
	int wRow, wCol, pRow, pCol, match;
	int cTerminate = wSize - pSize + 1;
	int wTerminate = (end < cTerminate) ? end : cTerminate;

	for (wRow = start; wRow <= wTerminate; wRow++){
		for (wCol = 1; wCol <= cTerminate; wCol++){
			match = 1;
//#ifdef DEBUGMORE
			printf("S:(%d, %d)\n", wRow-1, wCol-1);
//#endif
			for (pRow = 0; match && pRow < pSize; pRow++){
				for (pCol = 0; match && pCol < pSize; pCol++){
					if(world[wRow+pRow][wCol+pCol] != pattern[pRow][pCol]){
//#ifdef DEBUGMORE
						printf("\tF:(%d, %d) %c != %c\n", pRow, pCol,
						world[wRow+pRow][wCol+pCol], pattern[pRow][pCol]);
//#endif
						match = 0;    
					}
					if ((wRow + pRow) > (end + 1)) {
						printf("%d + %d is more than end %d\n", wRow, pRow, end);
						pRow = pSize;
						pCol = pSize;
						match = 2;
					}
				}
			}
			
			printf("(%d, %d) gets match = %d\n", wRow - 1, wCol - 1, match);
			
			if (match == 1){
				insertEnd(list, iteration, wRow-1, wCol-1, rotation);
//#ifdef DEBUGMORE
printf("*** Row = %d, Col = %d\n", wRow-1, wCol-1);
//#endif
			} else if (match == 2) {
				if (listToContinueFinding) {
					insertEnd(listToContinueFinding, iteration, wRow - 1, wCol - 1, rotation);
				}
			}
		}
	}
}

void continueSearch(char** world, int wSize, int iteration, char** pattern, int pSize, int rotation,
		MATCHLIST* list, MATCHLIST* listToResumeFrom, MATCHLIST* listToContinueFinding, int start, int end) {
    int i, match, pRow, pCol;
    MATCH* cur;
	
    if (list->nItem == 0) return;

    cur = list->tail->next;
    for( i = 0; i < list->nItem; i++, cur=cur->next) {
		match = 1;
		
		for (pRow = start - cur->row + 1; match && pRow < pSize; pRow++){
			for (pCol = 0; match && pCol < pSize; pCol++){
				if(world[cur->row + pRow + 1][cur->col + pCol + 1] != pattern[pRow][pCol]){
					match = 0;    
				}
				if ((cur->row + pRow + 1) > end) {
					pRow = pSize;
					pCol = pSize;
					match = 2;
				}
			}
		}
		if (match == 1){
			insertEnd(list, iteration, cur->row, cur->col, rotation);
		} else if (match == 2) {
			if (listToContinueFinding) {
				insertEnd(listToContinueFinding, iteration, cur->row, cur->col, rotation);
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

int isSmaller(int iter1, int rot1, int r1, int c1, int iter2, int rot2, int r2, int c2) {
	if (iter1 < iter2) return 1;
	if (iter1 > iter2) return 0;
	if (rot1 < rot2) return 1;
	if (rot1 > rot2) return 0;
	if (r1 < r2) return 1;
	if (r1 > r2) return 0;
	if (c1 < c2) return 1;
	return 0;
}

void printListSorted(MATCHLIST* list) {
    int i, x, iter, rot, r, c, takeCurrent;
    MATCH* cur;
    MATCH* designated;

    printf("List size = %d\n", list->nItem);    

    if (list->nItem == 0) return;

	for (x = 0; x < list->nItem; x++) {
		cur = list->tail->next;
		iter = -1;
		takeCurrent = 0;
		for( i = 0; i < list->nItem; i++, cur=cur->next){
			if (cur->iteration != -1) {
				if ((iter == -1) || (isSmaller(cur->iteration, cur->rotation, cur->row, cur->col, iter, rot, r, c))) {
					iter = cur->iteration;
					rot = cur->rotation;
					r = cur->row;
					c = cur->col;
					designated = cur;
				}
			}
			
			designated->iteration = -1;
			printf("%d:%d:%d:%d\n", iter, r, c, rot);
		}
	}
}

