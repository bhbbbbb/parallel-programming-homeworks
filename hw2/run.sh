#/src/bin/bash
mpic++ $1.cpp -o ./bin/$1 && mpiexec -f ../hosts -n $2 ./bin/$1 ${@:3}
