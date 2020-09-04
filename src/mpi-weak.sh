#!/bin/bash

for n_proc in {1..15}
do
    n_points=$(($n_proc * 10000))
    ./make-circle-input $n_points > input.in
    for ((i=0; i<10; i++))
    do
        mpirun -n $n_proc mpi-convex-hull < input.in >> mpi-weak.txt
    done
    cat mpi-weak.txt
    echo "Done with $n_proc threads and $n_points points"
done
