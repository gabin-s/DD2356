#include <stdio.h>
#include <string.h> /* For memset */ 
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>

#define NO_OUT
#define N_RUNS 100

#define DIM 2
#define X 0
#define Y 1
#define G 6.67430e-11
#define delta_t 0.05

typedef double vect_t[DIM];
double* masses;

vect_t *pos, *forces, *vel;

void compute_force(int q) {
    double x_diff, y_diff, dist, dist_cubed;
    vect_t force_qk;

    for(int k = 0; k < q; k++) {
        x_diff = pos[q][X] - pos[k][X]; 
        y_diff = pos[q][Y] - pos[k][Y]; 
        dist = sqrt(x_diff*x_diff + y_diff*y_diff); 
        dist_cubed = dist*dist*dist; 
        force_qk[X] = G*masses[q]*masses[k]/dist_cubed * x_diff; 
        force_qk[Y] = G*masses[q]*masses[k]/dist_cubed * y_diff;
        forces[q][X] += force_qk[X]; 
        forces[q][Y] += force_qk[Y]; 
        forces[k][X] -= force_qk[X]; 
        forces[k][Y] -= force_qk[Y]; 
    }

    pos[q][X] += delta_t*vel[q][X]; 
    pos[q][Y] += delta_t*vel[q][Y]; 
    vel[q][X] += delta_t/masses[q]*forces[q][X]; 
    vel[q][Y] += delta_t/masses[q]*forces[q][Y];
}

void init_bodies(int n) {
    pos      = malloc(n*sizeof(vect_t));
    vel      = malloc(n*sizeof(vect_t));
    forces   = malloc(n*sizeof(vect_t));
    masses   = malloc(n*sizeof(double));

    // initialize bodies
    for(int q = 0; q < n; q++) {
        pos[q][X] = (rand() / (double)(RAND_MAX)) * 2 - 1;
        pos[q][Y] = (rand() / (double)(RAND_MAX)) * 2 - 1;

        vel[q][X] = (rand() / (double)(RAND_MAX)) * 2 - 1;
        vel[q][Y] = (rand() / (double)(RAND_MAX)) * 2 - 1;

        masses[q] = fabs((rand() / (double)(RAND_MAX)) * 2 - 1);
    }
}

void dump_pos_vel(int n) {
    #ifndef NO_OUT
    for (int q = 0; q < n; q++)
        printf("%f,%f,%f,%f,", pos[q][X], pos[q][Y], vel[q][X], vel[q][Y]);
               
    printf("\n");
    #endif
}

void nsolve(int n, int n_steps) {
    // initialization
    init_bodies(n);

    // main loop
    for (int step = 1; step <= n_steps; step++) { 
        forces = memset(forces, 0, n*sizeof(vect_t));

        for (int q = 0; q < n-1; q++) 
            compute_force(q);

        dump_pos_vel(n);
    }
}

// function with timer                                                             
double mysecond(){
  struct timeval tp;
  struct timezone tzp;
  int i;

  i = gettimeofday(&tp,&tzp);
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

int main(int argc, char* argv[]) {
    int n = 10;
    int n_steps = 100;

    if(argc > 1) n = atoi(argv[1]);
    if(argc > 2) n_steps = atoi(argv[2]);
    if(argc > 3) {
        printf("usage: %s [n [n_steps]]\n", argv[0]);
        return -1;
    }

    // cold start
    nsolve(/* n */ n, /* n_steps */ n_steps);

    double t0 = mysecond();
    double times[N_RUNS];
    for(int i = 0; i < N_RUNS; i++) {
        nsolve(/* n */ n, /* n_steps */ n_steps);
        times[i] = mysecond();
    }
    
    double std = 0, mean = 0;
    double min = 1000, max = 0;

    double t = t0;
    for(int i = 0; i < N_RUNS; i++) {
        double dt = times[i] - t;
        t = times[i];
        printf("%f\n", dt);

        if(dt < min) min = dt;
        if(dt > max) max = dt;

        mean += dt;
        std  += dt*dt;
    }

    std  = sqrt((std - mean) / N_RUNS); 
    mean = mean / N_RUNS;

    printf("avg=%f, std=%f, min=%f, max=%f\n", mean, std, min, max);

    return 0;
}