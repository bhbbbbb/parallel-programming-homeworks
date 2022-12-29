#/src/bin/bash
gcc -g -Wall -fopenmp $1 -o ./bin/${1%%.*}.out && ./bin/${1%%.*}.out ${@:2}
