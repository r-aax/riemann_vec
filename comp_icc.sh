#!/bin/bash

FLAGS="-DINTEL -O3 -xmic-avx512"
INFO_FLAGS="-qopt-report=5"
EXE="riemann.out"

rm -f $EXE

#mpiicc \
#    -I src \
#    test/*.cpp src/*.cpp \
#    -DTEST_MODE=1 \
#    -DREPEATS=3 \
#    -DINNER_REPEATS=10 \
#    $FLAGS $INFO_FLAGS \
#    -lm -fopenmp \
#    -S

mpiicc \
    -I src \
    test/*.cpp src/*.cpp \
    $FLAGS \
    -DOPENMP_INTERLEAVE \
    -DTEST_MODE=1 \
    -DREPEATS=3 \
    -DINNER_REPEATS=10 \
    -o ${EXE} \
    -lm -fopenmp
