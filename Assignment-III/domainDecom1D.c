
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TAG 123

int main(int argc, char *argv[]){

    int rank, size, i, provided;
    
    // number of cells (global)
    int nxc = 128; // make sure nxc is divisible by size
    double L = 2*3.1415; // Length of the domain
    
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // number of nodes (local to the process): 0 and nxn_loc-1 are ghost cells 
    int nxn_loc = nxc/size + 3; // number of nodes is number cells + 1; we add also 2 ghost cells
    double L_loc = L/((double) size);
    double dx = L / ((double) nxc);
    
    // define out function
    double *f    = malloc(nxn_loc * sizeof(double));
    double *dfdx = malloc(nxn_loc * sizeof(double));

    for (i=1; i<(nxn_loc-1); i++)
        f[i] = sin(L_loc*rank + (i-1) * dx);
    
    // communicate and fill ghost cells f[0] and f[nxn_loc-1]
    int next = (rank == size - 1) ? 0 : (rank + 1);
    int prev  = (rank == 0) ? (size - 1):(rank - 1);

    // send f[2] to f[nxn_loc - 1] of the previous rank (send left, receive right)
    MPI_Sendrecv(&f[2], 1, MPI_DOUBLE, prev, TAG, 
                 &f[nxn_loc-1], 1, MPI_DOUBLE, next, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // send f[nxn_loc-3] to f[0] of the next rank (send right, receive left)
    MPI_Sendrecv(&f[nxn_loc-3], 1, MPI_DOUBLE, next, TAG, 
                 &f[0], 1, MPI_DOUBLE, prev, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // here we finish the calculations

    // calculate first order derivative using central difference
    // here we need to correct value of the ghost cells!

    for (i=1; i<(nxn_loc-1); i++)
        dfdx[i] = (f[i+1] - f[i-1])/(2*dx);

    if (rank==0) // print only rank 0 for convenience
        fprintf(stderr, "My rank %d of %d\n", rank, size );
    
    // print x and y, stop at n-3 to avoid duplicate points
    for (i=1; i<nxn_loc-2; i++)
	    printf("%f %f\n", L_loc*rank + (i-1) * dx, dfdx[i]);
    
    free(f);
    free(dfdx);

    MPI_Finalize();
}






