#!/bin/sh

N_RUNS=50

# constant parameters for ./nsolver
dt=0.05
n_steps=100

for npart in 500 1000 2000; do
    echo -n "$npart"
    for i in $(seq $N_RUNS); do
        out=$(./nsolver $npart $dt $n_steps)
        echo -n " ${out#Runtime: }"
    done
    echo
done