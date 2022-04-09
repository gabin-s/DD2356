#!/bin/sh

for num_threads in 1 2 4 8 16 24 28 32; do
    echo $num_threads
    OMP_NUM_THREADS=$num_threads ./sum 10000000
done