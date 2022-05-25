#!/bin/sh

N_RUNS=10

# constant parameters for ./nsolver
dt=0.05
n_steps=100

for npart in 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000; do
    echo -n "$npart"
    for i in $(seq $N_RUNS); do
        out=$(./nsolver $npart $dt $n_steps)
        echo -n " ${out#Runtime: }"
    done
    echo
done
