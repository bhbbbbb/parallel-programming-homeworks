#/src/bin/bash
# g++ -g -Wall -fopenmp $1 -o ./bin/${1%%.*}.out && ./bin/${1%%.*}.out ${@:2}
mpic++ -g -Wall -fopenmp $1 -o ./bin/${1%%.*}.out && mpiexec -f ../hosts -n $2 ./bin/${1%%.*}.out ${@:3}
