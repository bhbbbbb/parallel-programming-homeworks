#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

void hello();

int main(int argc, char** argv) {

    int num_threads = argc < 2 ? 1 : atoi(argv[1]);

    #pragma omp parallel num_threads(num_threads)
        hello();
    hello();

    return 0;
}

void hello() {
    int my_rank = omp_get_thread_num();
    int thread_count = omp_get_num_threads();

    printf("Hello from thread %d / %d\n", my_rank, thread_count);

    return;
}
