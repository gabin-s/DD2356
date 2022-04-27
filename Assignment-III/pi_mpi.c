
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

    if(rank == 0) {
        int val;
        long total_count = local_count;

        for(int i = 1; i < size; i++) {
            MPI_Recv(&val, 1, MPI_INT, /* src */ i, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_count += val;
        }

        // Estimate Pi and display the result
        pi = ((double)total_count / (double)NUM_ITER) * 4.0;

    }
    else {
        MPI_Send(&local_count, 1, MPI_INT, /* dest */ 0, TAG, MPI_COMM_WORLD);
    }

    stop_time = MPI_Wtime();
    
    if(rank == 0)
        printf("pi=%f, t=%f\n", pi, stop_time - start_time);        

    MPI_Finalize();
    return 0;
}