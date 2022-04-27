
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define TAG      1
#define SEED     921
#define NUM_ITER 1000000000

int main(int argc, char* argv[])
{
    int local_count = 0;
    double x, y, z, pi;

    int rank, size, provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // unique seed for each rank
    srand(rank * SEED);
    
    double start_time, stop_time, elapsed_time;
	start_time = MPI_Wtime();

    // Calculate PI following a Monte Carlo method
    for (int iter = 0; iter < NUM_ITER / size; iter++)
    {
        // Generate random (X,Y) points
        x = (double)random() / (double)RAND_MAX;
        y = (double)random() / (double)RAND_MAX;
        z = sqrt((x*x) + (y*y));
        
        // Check if point is in unit circle
        if (z <= 1.0) local_count++;
    }    

    int d = 0; // depth in the tree
    while(size >>= 1) {
        // receiver are the ranks multiple of 2^(d + 1)
        // then rank of sender fo a receiver `rank` is `rank + 2^d`
        int r = rank % (2 << d);

        if(r)
            MPI_Send(&local_count, 1, MPI_INT, rank - r, TAG, MPI_COMM_WORLD);
        else {
            int val, sender_rank = rank + (1 << d);
            MPI_Recv(&val, 1, MPI_INT, sender_rank, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            local_count += val;
        }
        d++;
    }

    stop_time = MPI_Wtime();
    
    if(rank == 0) {
        pi = ((double)local_count / (double)NUM_ITER) * 4.0;

        printf("pi=%f, t=%f\n", pi, stop_time - start_time);
    }

    MPI_Finalize();
    return 0;
}