#include <omp.h>
#include <stdio.h>

int main(int argc, char* argv[]) {
    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        printf("Hello world from thread %d!\n", id);
    }
    return 0;
}