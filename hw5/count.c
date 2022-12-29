#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <assert.h>

#define NUM_TO_SORT 100000
#define MAX_RAND_INT 10

void count_sort(const int* a, int n, int* out_arr, int num_threads) {

    int count;

    #pragma omp parallel for default(none) private(count) shared(a, out_arr, n) num_threads(num_threads)
    for (int i = 0; i < n; i++) {

        count = 0;

        for (int j = 0; j < n; j++)
            count += ( (a[j] < a[i]) || (a[j] == a[i] && j < i) );
        
        out_arr[count] = a[i];
    }

    return;
}

int cmp(const void* a, const void* b) {
    return *((int*)a) - *((int*)b);
}

int main(int argc, char** argv) {

    srand(0xAAAA);

    int num_threads = argc < 2 ? omp_get_max_threads() : atoi(argv[1]);

    int arr[NUM_TO_SORT], out[NUM_TO_SORT];

    for (int i = 0; i < NUM_TO_SORT; i++) arr[i] = rand() % MAX_RAND_INT;

    double count_time, qsort_time; 

    count_time = omp_get_wtime();
    count_sort(arr, sizeof(out) / sizeof(int), out, num_threads);
    count_time = omp_get_wtime() - count_time;

    qsort_time = omp_get_wtime();
    qsort(arr, sizeof(arr) / sizeof(int), sizeof(int), cmp);
    qsort_time = omp_get_wtime() - qsort_time;

    printf("num_threads = %d, n = %d\n", num_threads, NUM_TO_SORT);
    printf("count_sort: %lf s.\nqsort: %lf s.\n", count_time, qsort_time);

    for (int i = 0; i < NUM_TO_SORT; i++) {
        assert(arr[i] == out[i]);
    }

    printf("\n");

    return 0;
}
