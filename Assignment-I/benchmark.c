#define N 5000
#define NRUNS 1000
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h> 

double mysecond();

int main(){
  int i, j;
  double t1, t2; // timers                                                         
  double a[N], b[N], c[N]; // arrays  
  
  // init arrays                                                                   
  for (i = 0; i < N; i++){
    a[i] = 47.0;
    b[i] = 3.1415;
  }

  // prevent cold start
  for(i = 0; i < N; i++)
    c[i] = a[i]*b[i];

  // measure performance
  t1 = mysecond();
  for(j = 0; j < NRUNS; j++) {
    for(i = 0; i < N; i++) {
      c[i] = a[i]*b[i];
    }
  }
  t2 = mysecond();

  // trick smart compiler to avoid marking code above as dead
  for(i = 0; i < N; i++)
    if(c[i] * c[i] == 0) printf("coucou");

  printf("Execution time: %11.8f s\n", (t2 - t1) / NRUNS);
  return 0;
}

// function with timer                                                             
double mysecond(){
  struct timeval tp;
  struct timezone tzp;
  int i;

  i = gettimeofday(&tp,&tzp);
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}
