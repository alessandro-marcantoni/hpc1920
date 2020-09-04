#!/bin/bash

for n_proc in {1..15}
do
    for ((i=0; i<10; i++))
    do
        mpirun -n $n_proc mpi-convex-hull < circ100k.in >> mpi-strong.txt
    done
    echo "Done with $n_proc threads"
done
