#!/bin/bash

FLAGS="-DINTEL -O3 -xmic-avx512 -qopt-prefetch -fno-inline"
INFO_FLAGS="-qopt-report=5"
EXE="riemann.out"

rm -f $EXE

mpiicc \
    -DDEBUG \
    -I src \
    test/*.cpp src/*.cpp \
    $FLAGS $INFO_FLAGS \
    -lm -fopenmp \
    -S

mpiicc \
    -DDEBUG \
    -I src \
    test/*.cpp src/*.cpp \
    $FLAGS \
    -o $EXE \
    -lm -fopenmp
