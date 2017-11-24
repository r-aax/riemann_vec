#!/bin/sh

set -x

rm -f riemann.out

g++ \
    ./src/*.cpp ./test/test.cpp \
    -I./src \
    -lm -fopenmp \
    -o riemann.out

./riemann.out