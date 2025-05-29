#!/bin/sh

set -x

rm -f riemann.out

g++ \
    ./src/*.cpp ./test/test.cpp \
    -DOPENMP_CHUNKS \
    -DTEST_MODE=0 \
    -DREPEATS_ORIG=3 \
    -DREPEATS_OPT=3 \
    -DINNER_REPEATS=10 \
    -I./src \
    -O3 \
    -lm -fopenmp \
    -o riemann.out

./riemann.out
