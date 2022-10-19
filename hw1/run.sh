#!/bin/bash
mpicc $1.c -o $1.out && mpiexec -f ../hosts -n $2 ./$1.out ${@:3}
