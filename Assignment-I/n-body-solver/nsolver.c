#include <stdio.h>
#include <string.h> /* For memset */ 
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

#define DIM 2
#define X 0
#define Y 1
#define G 6.673e-11

typedef double vect_t[DIM];

// function with timer                                                             
double mysecond(){
  struct timeval tp;
  struct timezone tzp;
  int i;

  i = gettimeofday(&tp,&tzp);
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

/**
 * @brief Compute and stores in `forces` the force vectors for all particles.
 * 
 * @param N number of particles
 * @param forces array where calculated force vectors will be stored
 * @param pos position vector of all particles
 * @param masses masses of all particles
 */
void compute_forces(int N, vect_t* forces, vect_t *pos, double* masses) {
    double x_diff, y_diff, dist, dist_cubed;

    // reset forces
    memset(forces, 0, sizeof(double) * N);

    vect_t force_qk;

    for(int q = 0; q < N; q++) {
#ifdef REDUCED
        for(int k = q + 1; k < N; k++) {
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
#else
        for(int k = 0; k < N; k++) {
            if(k == q) continue;

            x_diff = pos[q][X] - pos[k][X];
            y_diff = pos[q][Y] - pos[k][Y];
            dist = sqrt(x_diff*x_diff + y_diff*y_diff);
            dist_cubed = dist*dist*dist;
            forces[q][X] -= G*masses[q]*masses[k]/dist_cubed * x_diff;
            forces[q][Y] -= G*masses[q]*masses[k]/dist_cubed * y_diff;
        }
#endif
    }
}

/**
 * @brief Move all particles
 * 
 * @param N total number of particles
 * @param delta_t the time delta used in the simulation
 * @param vel velocity vectors for all particles
 * @param forces force vectors for all particles
 * @param masses masses of all particles
 */
void move_particles(int N, double dt, vect_t* pos,vect_t* vel, vect_t* forces, double* masses) {
    for(int q = 0; q < N; q++) {
        pos[q][X] += dt*vel[q][X]; 
        pos[q][Y] += dt*vel[q][Y]; 
        vel[q][X] += dt/masses[q]*forces[q][X]; 
        vel[q][Y] += dt/masses[q]*forces[q][Y];
    }
}

void init_random(int N, vect_t *pos, vect_t *vel, double *masses) {
    for(int q = 0; q < N; q++) {
        pos[q][X] = (rand() / (double)(RAND_MAX)) * 2 - 1;
        pos[q][Y] = (rand() / (double)(RAND_MAX)) * 2 - 1;

        vel[q][X] = (rand() / (double)(RAND_MAX)) * 2 - 1;
        vel[q][Y] = (rand() / (double)(RAND_MAX)) * 2 - 1;

        masses[q] = fabs((rand() / (double)(RAND_MAX)) * 2 - 1);
    }
}

void print_usage(char* prgm_name) {
    fprintf(stderr, "usage: %s <N> <dt> <n_steps>\n", prgm_name);
}

int main(int argc, char* argv[]) {
    int N, n_steps;
    double dt; // time step
    double t1, t2; // used for runtime measurement

    // parse arguments
    if(argc != 4) {
        print_usage(argv[0]);
        return 1;
    }
    N  = atoi(argv[1]);
    dt = atof(argv[2]);
    n_steps = atoi(argv[3]);

    // allocate memory
    vect_t* forces = malloc(N * sizeof(vect_t));
    vect_t* vel    = malloc(N * sizeof(vect_t));
    vect_t* pos    = malloc(N * sizeof(vect_t));
    double* masses = malloc(N * sizeof(double));
    
    if(forces == NULL || vel == NULL || pos == NULL || masses == NULL) {
        fprintf(stderr, "unable to allocate memory\n");
        return 2;
    }

    t1 = mysecond();

    // random initialization
    init_random(N, pos, vel, masses);

    // actual simulation
    for(int step = 0; step <= n_steps; step++) {
        compute_forces(N, forces, pos, masses);
        move_particles(N, dt, pos, vel, forces, masses);
    }

    t2 = mysecond();
    printf("Runtime: %f\n", t2 - t1);

    // free memory
    free(forces);
    free(vel);
    free(pos);
    free(masses);

    return 0;
}
