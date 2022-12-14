#/src/bin/bash
g++ $1.cpp -o ./bin/$1.out -lpthread && ./bin/$1.out ${@:2}
