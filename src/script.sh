#!/bin/bash

for input in inputs/*.in
do
    mpirun -n 1 mpi-convex-hull < $input > /dev/null
    mpirun -n 2 mpi-convex-hull < $input > /dev/null
done
