#!/bin/sh

set -x

rm -f riemann.out
rm -f gmon.out
rm -f *.gcda
rm -f *.gcno
rm -f *.gcov
rm -f src/*.gcda
rm -f src/*.gcno
rm -f src/*.gcov

g++ \
    ./src/*.cpp ./test/test.cpp \
    -I./src \
    -g -pg \
    -fprofile-arcs -ftest-coverage \
    -lm -fopenmp \
    -o riemann.out

./riemann.out

gprof riemann.out gmon.out -p > prof_p.txt
gprof riemann.out gmon.out -q > prof_q.txt
gprof riemann.out gmon.out -A > prof_A.txt

mv riemann.gcda src
mv riemann.gcno src

gcov ./src/riemann.c