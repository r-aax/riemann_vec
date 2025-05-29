#!/bin/bash

FLAGS="-DINTEL -O3 -xmic-avx512"
INFO_FLAGS="-qopt-report=5"
EXE="riemann.out"

rm -f $EXE

mpiicc \
    -I src \
    test/*.cpp src/*.cpp \
    $FLAGS \
    -DOPENMP_CHUNKS \
    -DTEST_MODE=1 \
    -DREPEATS_ORIG=1 \
    -DREPEATS_OPT=10 \
    -DINNER_REPEATS=100 \
    -o ${EXE} \
    -lm -fopenmp
