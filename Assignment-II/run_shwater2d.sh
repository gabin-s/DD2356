#!/bin/sh

for nthreads in 1 2 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64; do
    echo -n "$nthreads "
    OMP_NUM_THREADS=$nthreads ./shwater2d 1000 | grep 'Solver took' | awk '{ print $3 }'
done