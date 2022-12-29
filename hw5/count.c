#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <assert.h>

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
    int num_to_sort = argc < 3 ? 100 : atoi(argv[2]);
    int max_num = argc < 4 ? 10 : atoi(argv[3]);

    // int arr[num_to_sort], out[num_to_sort];
    int* arr = (int*)malloc(sizeof(int) * num_to_sort);
    int* out = (int*)malloc(sizeof(int) * num_to_sort);

    for (int i = 0; i < num_to_sort; i++) arr[i] = rand() % max_num;

    double count_time, qsort_time; 

    count_time = omp_get_wtime();
    count_sort(arr, num_to_sort, out, num_threads);
    count_time = omp_get_wtime() - count_time;

    qsort_time = omp_get_wtime();
    qsort(arr, num_to_sort, sizeof(int), cmp);
    qsort_time = omp_get_wtime() - qsort_time;

    printf("num_threads = %d, n = %d\n", num_threads, num_to_sort);
    printf("count_sort: %lf s.\nqsort: %lf s.\n", count_time, qsort_time);

    for (int i = 0; i < num_to_sort; i++) {
        assert(arr[i] == out[i]);
        printf("%d, ", arr[i]);
    }

    printf("\n");
    free(arr);
    free(out);

    return 0;
}
