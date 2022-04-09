#include <stdlib.h>
#include <omp.h>
#include <stdio.h>
#include <sys/time.h>
#include <math.h>

#define NRUNS 100

// function with timer                                                             
double mysecond(){
  struct timeval tp;
  struct timezone tzp;
  int i;

  i = gettimeofday(&tp,&tzp);
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

void generate_random(double *input, size_t size)
{
  for (size_t i = 0; i < size; i++) {
    input[i] = rand() / (double)(RAND_MAX);
  }
}

double serial_sum(double *x, size_t size)
{
    double sum_val = 0.0;

    for (size_t i = 0; i < size; i++) {
        sum_val += x[i];
    }

    return sum_val;
}


double omp_sum(double *x, size_t size) {
    double res;

    #pragma omp parallel for
    for(size_t i = 0; i < size; i++) {
        res += x[i];
    }

    return res;
}

double omp_critical_sum(double *x, size_t size) {
    double res;

    #pragma omp parallel for
    for(size_t i = 0; i < size; i++) {
        #pragma omp critical
        res += x[i];
    }

    return res;
}


double omp_local_sum_best(double *x, size_t size) {
    const int max_threads = omp_get_max_threads();

    /* Compute partial sums (max_threads partial sums) */
    double res;

    // create an union to pad local sum in memory, to avoid false sharing
    union {
        double val;
        char pad[128];
    } local_sum[max_threads];

    #pragma omp parallel shared(local_sum)
    {
        int id = omp_get_thread_num();
        local_sum[id].val = 0;
        
        #pragma omp for
        for(size_t i = 0; i < size; i++) {
            local_sum[id].val += x[i];
        }
    }

    /* Compute the full sum from partial results */
    res = 0;
    for(size_t i = 0; i < max_threads; i++) {
        res += local_sum[i].val;
    }

    return res;
}

double omp_local_sum(double *x, size_t size) {
    const int max_threads = omp_get_max_threads();

    /* Compute partial sums (max_threads partial sums) */
    double local_sum[max_threads], res;

    #pragma omp parallel shared(local_sum)
    {
        int id = omp_get_thread_num();
        local_sum[id] = 0;
        
        #pragma omp for
        for(size_t i = 0; i < size; i++) {
            local_sum[id] += x[i];
        }
    }

    /* Compute the full sum from partial results */
    res = 0;
    for(size_t i = 0; i < max_threads; i++) {
        res += local_sum[i];
    }

    return res;
}

void evaluate_fct(double (*f)(double*, size_t), double* x, size_t n) {
    /* Measure performance of serial code */
    f(x, n); // avoid cold run

    double ts[NRUNS], res;
    double t0 = mysecond();
    for(int i = 0; i < NRUNS; i++) {
        res = f(x, n);
        ts[i] = mysecond();
    }
    
    double s, s2, t, dt;
    s = s2 = 0;
    t = t0;
    for(int i = 0; i < NRUNS; i++) {
        dt = ts[i] - t;
        t = ts[i];
        
        s  += dt;
        s2 += dt*dt;
    }

    double mean = s / NRUNS;
    double std  = sqrt(s2 / NRUNS - mean*mean);

    printf("res=%f, mean=%f, std=%f\n", res, mean, std);
}

int main(int argc, char* argv[]) {
    size_t n = 1024;
    if(argc > 1) n = atoi(argv[1]);
    if(argc > 2) {
        printf("usage: %s [size]");
        return 1;
    }

    printf("size: %d, max threads: %d\n", n, omp_get_max_threads());
    
    double *x, res;

    x = malloc(n * sizeof(double));
    generate_random(x, n);

    printf("serial_sum:    ");
    evaluate_fct(serial_sum, x, n);

    printf("omp_sum:       ");
    evaluate_fct(omp_sum, x, n);

    printf("omp_critical_sum: ");
    evaluate_fct(omp_critical_sum, x, n);

    printf("omp_local_sum: ");
    evaluate_fct(omp_local_sum, x, n);

    printf("omp_local_sum_best: ");
    evaluate_fct(omp_local_sum_best, x, n);

    free(x);
    return 0;
}