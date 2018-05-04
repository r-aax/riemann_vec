#!/bin/bash

EXE="riemann.out"

if [[ -f "$EXE" ]]
then
    srun -p knl -n 1 ./${EXE} > _result.txt
    cat _result.txt
    cat _result.txt >> result.txt
    rm -f _result.txt
fi
