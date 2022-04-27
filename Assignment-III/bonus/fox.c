#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define TAG_SENDRECV 123
#define TAG_SYNC 456

#define MPI_NUMBER_T MPI_DOUBLE
typedef double number_t;

void fail(int code) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == 0)
        fprintf(stderr, "FATAL ERROR: %d\n", code);

    MPI_Finalize();
    exit(code);
}

void generate_id(int m, number_t *A) {
    for(int i = 0; i < m*m; i++)
        A[i] = 0;

    // fill diagonal with ones
    for(int i = 0; i < m*m; i += m+1)
        A[i] = 1;
}

void generate_rand(int m, number_t *A) {
    for(int i = 0; i < m*m; i++)
        A[i] = (int) (10.0 * rand() / RAND_MAX);
}

void matmul(int m, number_t *A, number_t *B, number_t *dst) {
    for(int mi = 0; mi < m*m; mi += m) {
        for(int j = 0; j < m; j++) {
            dst[mi + j] = 0;
            for(int k = 0; k < m; k++)
                dst[mi + j] += A[mi + k] * B[m*k + j];
        }
    }
}

void matadd(int m, number_t *A, number_t *B, number_t *dst) {
    for(int mi = 0; mi < m*m; mi += m) {
        for(int j = 0; j < m; j++) {
            dst[mi + j] = A[mi + j] + B[mi + j];
        }
    }
}


// compute AB + C and stores it in C
void matmuladd(int m, number_t *A, number_t *B, number_t *C) {
    for(int mi = 0; mi < m*m; mi += m) {
        for(int j = 0; j < m; j++) {
            for(int k = 0; k < m; k++)
                C[mi + j] += A[mi + k] * B[m*k + j];
        }
    }
}

void print_matrix(int m, number_t *A) {
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < m; j++)
            printf("%.2f ", A[m*i + j]);
        printf("\n");
    }
}

/**
 * @brief Copy the tile `rank` of size m x m from `src` to `dst`.
 * 
 * @param m the size of the tile
 * @param rank the index of the tile in `src` counted row-wise
 * @param dst the m x m destination
 * @param src the m*m x m*m source
 */
void tile_copy(int m, int s_N, int rank, number_t *dst, number_t *src) {
    div_t qr = div(rank, s_N);
    int offset = m*(qr.rem + s_N*m*qr.quot);

    for(int i = 0; i < m; i++)
        memcpy(dst + m*i, src + offset + m*s_N*i, m*sizeof(number_t));
}

// replace the tile for rank `rank` in `dst`
void tile_replace(int m, int s_N, int rank, number_t *dst, number_t *src) {
    div_t qr = div(rank, s_N);
    int offset = m*(qr.rem + s_N*m*qr.quot);

    for(int i = 0; i < m; i++)
        memcpy(dst + offset + m*s_N*i, src + m*i, m*sizeof(number_t));
}

void mpi_matmul(int M, int s_N, int rank, number_t *A, number_t *B) {
    const int m = M / s_N;
    const int N = s_N * s_N;


    // allocate buffers
    number_t* A_tile = (number_t*)malloc(m * m * sizeof(number_t));
    number_t* B_tile = (number_t*)malloc(m * m * sizeof(number_t));
    number_t* C_tile = (number_t*)calloc(m * m, sizeof(number_t));
    number_t* tmp    = (number_t*)malloc(m * m * sizeof(number_t));

    if(A_tile == NULL || B_tile == NULL || C_tile == NULL || tmp == NULL) 
        fail(1);

    // get my tiles from A and B
    tile_copy(m, s_N, rank, A_tile, A);
    tile_copy(m, s_N, rank, B_tile, B);

    printf("r%d, m=%d, s_N=%d, A=%f, B=%f, C=%f\n", rank, m, s_N, *A_tile, *B_tile, *C_tile);
    print_matrix(m, A_tile);

    // create a new communicator for the row
    div_t colorkey = div(rank, s_N);
    int color = colorkey.quot;
    int key   = colorkey.rem;
    
    MPI_Comm commRow;
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &commRow);

    // rank for above and below
    int above = (rank - s_N) % N; if(above < 0) above += N;
    int below = (rank + s_N) % N;

    // for each diagonal
    for(int i=0; i < s_N; i++) {        
        // 1. broadcast A row-wise for the i-th diagonal
        int root = (i + color) % s_N;
        if(key == root)
            memcpy(tmp, A_tile, m*m*sizeof(number_t));

        // the received A is saved in tmp
        MPI_Bcast(tmp, m*m, MPI_NUMBER_T, root, commRow);

        // multiply received A_tile (in tmp) with my own B_tile
        // and add the result to the calculated C_tile
        matmuladd(m, tmp, B_tile, C_tile);

        // 2. send B to above and receive B from below
        MPI_Sendrecv_replace(
            B_tile, m*m, MPI_NUMBER_T, above, TAG_SENDRECV, below, TAG_SENDRECV, 
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    free(A_tile);
    free(B_tile);
    free(tmp);

    MPI_Send(C_tile, m*m, MPI_NUMBER_T, /* dst */ 0, TAG_SYNC, MPI_COMM_WORLD);
    MPI_Comm_free(&commRow);

    free(C_tile);
}

int main(int argc, char* argv[]) {
    number_t *A, *B; // input matrices
    int N;   // number of processes
    int M;   // number of rows/columns in the input matrices
    int s_N; // number of rows/columns in the process grid
    int m; // number of rows/columns in each tile (== M / s_n)

    int rank, provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
    MPI_Comm_size(MPI_COMM_WORLD, &N);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(argc < 2) {
        if(rank == 0)
            fprintf(stderr, "usage: %s M\n", argv[0]);
        fail(2);
    }

    /* ensure the number of processes is a square */
    float sqrt_size = sqrt(N);
    s_N = (int) sqrt_size;
    
    if(sqrt_size != s_N) {
        if(rank == 0)
            fprintf(stderr, "You must run a number of processes which is a square, %d is not.\n", N);
        
        MPI_Finalize();
        exit(0);
    }
    

    // -- generate random input matrices --
    M = atoi(argv[1]);
    A = (number_t*) malloc(M * M * sizeof(number_t));
    B = (number_t*) malloc(M * M * sizeof(number_t));
    if(A == NULL || B == NULL) fail(3);

    generate_rand(M, A);
    generate_rand(M, B);
    // -- -- 

    if(M % s_N || M < s_N) {
        if(rank == 0)
            fprintf(stderr, "invalid matrix size: M=%d is not "
                "compatible with a processor grid width of %d.\n", M, s_N);
        
        MPI_Finalize();
        exit(0);
    }
    m = M / s_N;

    double start_time, stop_time, elapsed_time;

    number_t *C_tiles[N];
    MPI_Request *requests;

    if(rank == 0) {
        requests = (MPI_Request*) malloc(N*sizeof(MPI_Request));
        if(requests == NULL) fail(8);

        for(int src = 0; src < N; src++) {
            C_tiles[src] = (number_t*) malloc(m*m*sizeof(number_t));
            if(C_tiles[src] == NULL) fail(7);

            MPI_Irecv(C_tiles[src], m*m, MPI_NUMBER_T, src, TAG_SYNC, MPI_COMM_WORLD, &requests[src]);
        }
    }

    start_time = MPI_Wtime();
    mpi_matmul(M, s_N, rank, A, B);

    if(rank == 0) {
        printf("--- A\n");
        print_matrix(M, A);

        // output matrix
        number_t *C = (number_t*) malloc(M*M*sizeof(number_t));
        if(C == NULL) fail(5);

        MPI_Waitall(N, requests, MPI_STATUSES_IGNORE);

        for(int src = 0; src < N; src++) {
            // replace tile in C
            tile_replace(m, s_N, src, C, C_tiles[src]);
            free(C_tiles[src]);
        }

        stop_time = MPI_Wtime();

        elapsed_time = stop_time - start_time;
        printf("elapsed: %fs\n", elapsed_time);

        number_t *C_brutal = (number_t*) malloc(M*M*sizeof(number_t));
        // brutal multiplication
        start_time = MPI_Wtime();
        matmul(M, A, B, C_brutal); // reuse C_tiles
        stop_time = MPI_Wtime();

        elapsed_time = stop_time - start_time;
        printf("elaspsed (brutal): %fs\n", elapsed_time);

        int cmp = memcmp(C, C_brutal, M*M*sizeof(number_t));
        if(cmp)
            fprintf(stderr, "bad computation\n");
       
        free(C_brutal);
        free(C);
    }

    free(A);
    free(B);

    MPI_Finalize();
    return 0;
}