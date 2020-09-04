#!/bin/bash

for n_proc in {1..15}
do
    for ((i=0; i<10; i++))
    do
        OMP_NUM_THREADS=$n_proc ./omp-convex-hull < circ100k.in >> omp-strong.txt
    done
done
