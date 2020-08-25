#!/bin/bash

EXE="riemann.out"

if [[ -f "$EXE" ]]
then
    mpirun -s knl -np 1 riemann.out
    #srun -p knl -n 1 --constraint=cache,a2a ./${EXE} > _result.txt
    #cat _result.txt
    #cat _result.txt >> result.txt
    #rm -f _result.txt
fi
