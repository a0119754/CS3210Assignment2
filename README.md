# CS3210Assignment2

Personal repository to facilitate transfer of data between multiple machines for the second assignment of CS3210 Parallel Computing.

+ To compile SETL.c: `gcc -o SETL SETL.c`
+ To execute SETL: `./SETL [filename of world file] [number of iterations] [filename of pattern file]`
+ To compile SETL_par.c: `mpicc SETL_par.c -o SETL_par`
+ To execute SETL_par: `mpirun -machinefile machinefile.lab -rankfile rankfile.lab -np [number_of_processes] ./SETL_par [filename of world file] [number of iterations] [filename of pattern file]`
+ To send a file through SCP: `scp [filename] [account_name]@[host_name]:~`

## Rough idea on what to do:

1. Distribute the four rotated positions of the single pattern to all processes
  * A naive approach would be to have master do the rotating, then send all 4 to each process
  * May want to spend some time on this stage experimenting and figuring out which way would be faster to rotate while sending (maybe send to 1, who will send to 2 who send to 3 who send to 4, etc., meanwhile 1 rotates and send to 2 to 3 to 4 etc., then repeat for two more rotations?)
  * Master does not need to know how the rotated patterns look like (or does it? Do we want to involve the master in searching? Or nah, master should evolve while the rest search)
2. For each iteration:
  1. Split up the world into n portions for n processes ~~(I recommend splitting them horizontally, so process 1 takes top x%, process 2 takes next x%, etc., and process n takes bottom x%, where x=100/n)~~ (How about first checking the number of rows and columns, then act accordingly? ie. if number of columns is more than number of rows then split vertically, else split horizontally?)
  2. Give each process a portion of the world
  3. Master evolves current iteration into the next, while the slaves search for aliens (if have time, sort results obtained in point 2.7?)
  4. Slave starts searching starting each cell which most top-left belongs to them (or is this efficient? Bottom would have lots of waiting, but it's meaningless for two processors to start searching for the same halves when one has already invalidated it). Incompleted results is sent to the next slave
  5. Slave only returns results to master, as well as send results of incomplete searches to the next slave, after:
    1. Searching his own map (for aliens which most top-left cell lies in their own map)
    2. Receiving results for incomplete searches from previous slave and continuing on the search for possible aliens
  6. Master receives n results from n processes
  7. As master has finished evolving the world, we proceed to the next iteration (results need to be sorted as well, perhaps sort it in point 2.3?)
3. Print results, print time taken to complete, job done \o/

P.S.: May want to print out to inform the completion of a job, along with timestamp, to have a rough idea of which process finished its work and is idling, while the others are still struggling to complete theirs. Gaining insight is important.

- To compile: `mpicc [filename].c -o [filename]`
- To run: `mpirun -np [number_of_processes] ./[filename]`
- To run on on remote node: `mpirun -H soctf-pdc-[number] -np [number_of_processes] ./[filename]` (Note: You need the compiled binary in that remote node)

The machine configuraiton file (machinefile.lab) should specify which nodes we want the processes to execute in. mpirun cycles through the listed nodes in round-robin fashion by default. For example, to specify running on nodes 001 and 002:
- soctf-pdc-001
- soctf-pdc-002
(Execute using: mpirun -machinefile machinefile.lab -np [number_of_processes] ./[filename])

The following code allows you to determine which node which process is running on:
```
char hostname[256];
int sc_status;
memset(hostname, 0, sizeof(hostname));
sc_status = gethostname(hostname, sizeof(hostname)-1);
if (sc_status)
{
	perror("gethostname");
	return sc_status;
}
printf("Process %d is on hostname %s\n", rank, hostname);
```

Rank file (combined with a machine file) allows for more precise control of the mapping of processes to nodes:
```
rank [number]=[hostname] slot=[socket_range]:[core_range]
rank 0=soctf-pdc-001 slot=0:0
rank 1=soctf-pdc-002 slot=0:0
rank 2=soctf-pdc-001 slot=0:0
rank 3=soctf-pdc-002 slot=0:0-2
rank 4=soctf-pdc-002 slot=0:0-2
rank 5=soctf-pdc-002 slot=0:0-2
```
(To execute: mpirun -machinefile machinefile.lab -rankfile rankfile.lab -np [number_of_processes] ./[filename])

- MPI_Send() and MPI_Recv() are blocking (caller is stopped until message is successfully sent/received)
- MPI_Isend() and MPI_Irecv() are non-blocking (caller continues executing instructions) (MPI_Test() checks if communication has finished, and MPI_Wait() blocks until it has finished)

## Lab system information

* Jetson - 4 CPUs
* i-7 - 4 cores x 2 threads each = 8 CPUs
* i-5 - 4 CPUs

## To do:

1. When processes >= 5, results passing on from previous guy and to next guy seems erratic
2. When using more than one cluster, results seem erratic as well
