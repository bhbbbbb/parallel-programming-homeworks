#/src/bin/bash
g++ -g -Wall -fopenmp $1 -o ./bin/${1%%.*}.out && ./bin/${1%%.*}.out ${@:2}
