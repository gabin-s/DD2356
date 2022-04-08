#!/bin/sh

for num_threads in 1 2 4 8 12 16 20 24 28 32; do
    for i in 1 2 3 4 5; do
        echo -n "$num_threads "
        OMP_NUM_THREADS=$num_threads ./stream | grep 'Copy:' | awk '{ print $2 }'
    done
done