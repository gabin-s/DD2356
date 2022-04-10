#!/bin/sh

for i in $(seq 100); do
    OMP_NUM_THREADS=32 ./dftw | grep seconds | awk '{print $4}'
done