#!/bin/sh

for nthreads in 1 2 4 8 12 16 20 24 28 32; do
    for i in $(seq 10); do
        echo -n "$nthreads "
        OMP_NUM_THREADS=$nthreads ./dftw | grep seconds | awk '{print $4}'
    done
done