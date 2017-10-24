#!/bin/sh

rm -f riemann.out
g++ ./src/*.cpp -I./src -o riemann.out
./riemann.out