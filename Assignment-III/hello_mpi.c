#include <mpi.h>
#include <stdio.h>

#define TAG_A 0

int main(int argc, char* argv[]) {
    
    int rank, size, i, provided;
    float A[10];

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0) {
        for(i = 0; i < 10; i++) A[i] = i;

        for(i=1; i < size; i++)
            MPI_Send(&A, 10, MPI_FLOAT, i, TAG_A, MPI_COMM_WORLD);
    }
    else {
        MPI_Recv(&A, 10, MPI_FLOAT, 0, TAG_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    for(i = 0; i < 10; i++)
        A[i] *= rank;

    printf("My rank is %d\n", rank);
    printf("My values for A: ");
    for(i = 0; i < 10; i++) printf("%.2f, ", A[i]);
    printf("\n");

    MPI_Finalize();
    return 0;
}