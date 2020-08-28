#!/bin/sh

set -x

rm -f riemann.out

g++ \
    ./src/*.cpp ./test/test.cpp \
    -DOPENMP_RACE \
    -DTEST_MODE=0 \
    -DREPEATS=3 \
    -DINNER_REPEATS=10 \
    -I./src \
    -O3 \
    -lm -fopenmp \
    -o riemann.out

./riemann.out
