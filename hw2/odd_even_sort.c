// #define DEBUG
// #define BLOCKING
#define OUTPUT
#define MAXINT RAND_MAX

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

#ifndef BLOCKING
#define BOUND(chunk_size) (chunk_size - 2)
#define EVEN_START 2
#else
#define BOUND(chunk_size) (chunk_size - 1)
#define EVEN_START 0
#endif

/**
 * @param next_id be negative if there is no next
 * @param prev_id be negative if there is no previous
*/
int odd_even_sort(int* chunk_arr, int chunk_size, int next_id, int prev_id, int tag);

static inline int randint() {
    return rand() % MAXINT;
}

static inline void swap(int *arr, int x, int y) {
    int tmp = arr[x];
    arr[x] = arr[y];
    arr[y] = tmp;
}

#ifdef DEBUG
void print_chunk_arr(int id, int* chunk_arr, int chunk_size);
#endif

int main(int argc, char *argv[]) {


    int id, comm_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    // init random seeds
    srand(0XAAAA + id);

    int n;
    if (id == 0) {
        printf("Please enter an integer N to sort: ");
        fflush(stdout);
        scanf("%d", &n);
    }

    double start_time = 0.0;
    if (id == 0) start_time = MPI_Wtime();

    MPI_Bcast(&n, 1, MPI_INT32_T, 0, MPI_COMM_WORLD);

    int *chunk_sizes = (int*)malloc(comm_size * sizeof(int));
    int *displacements = (int*)malloc(comm_size * sizeof(int));
    int *chunk_arr = (int*)malloc(n * sizeof(int));

    int d = 0;
    for (int i = 0; i < comm_size; i++) {
        chunk_sizes[i] = (n / comm_size) + (int)(i < n % comm_size);
        displacements[i] = d;
        d += chunk_sizes[i];
    }
    
    for (int i = 0; i < chunk_sizes[id]; i++) chunk_arr[i] = randint();

    #ifdef DEBUG
    for (int i = 0; i < comm_size; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (id == i)
            print_chunk_arr(id, chunk_arr, chunk_sizes[id]);
    }
    #endif

    int prev_id = (id - 1 >= 0 && chunk_sizes[id - 1] > 0) ? id - 1 : -1;
    int next_id = (id + 1 < comm_size && chunk_sizes[id + 1] > 0) ? id + 1 : -1;
    
    int sorted = 0;
    int sorted_recv_buf;
    int t = 0;
    while (!sorted) {

        sorted = 1;
        if (chunk_sizes[id] >= 1)
            sorted = odd_even_sort(chunk_arr, chunk_sizes[id], prev_id, next_id, t++);
        
        MPI_Allreduce(&sorted, &sorted_recv_buf, 1, MPI_INT32_T, MPI_BAND, MPI_COMM_WORLD);
        sorted = sorted_recv_buf;

    }

    #ifdef DEBUG
    for (int i = 0; i < comm_size; i++) {
        MPI_Barrier(MPI_COMM_WORLD);

        if (id == i) print_chunk_arr(id, chunk_arr, chunk_sizes[id]);
    }
    #endif

    int *sorted_arr = (int*)malloc(n * sizeof(int));
    MPI_Gatherv(
        chunk_arr, chunk_sizes[id], MPI_INT32_T,
        sorted_arr, chunk_sizes, displacements, MPI_INT32_T, 0, MPI_COMM_WORLD
    );
    MPI_Barrier(MPI_COMM_WORLD);
    if (id == 0) {
        fprintf(stdout, "\n--------- results(in %f sec.) ---------------\n", MPI_Wtime() - start_time);
        fflush(stdout);
#ifdef OUTPUT
        for (int i = 0; i < n; i++) printf("%d, ", sorted_arr[i]);
#endif
        printf("\n");
    }

    free(sorted_arr);
    free(chunk_sizes);
    free(chunk_arr);
    MPI_Finalize();
    return 0;
}

// #endif

/**
 * @param next_id be negative if there is no next
 * @param prev_id be negative if there is no previous
*/
int odd_even_sort(int* chunk_arr, int chunk_size, int prev_id, int next_id, int tag) {

    int sorted = 1;
    int tem;
    int next_val;
    int prev_val;

    // check boundaries
    #ifndef BLOCKING

    MPI_Request req_buff[4];
    if (prev_id >= 0) {
        MPI_Isend(&chunk_arr[0], 1, MPI_INT32_T, prev_id, tag, MPI_COMM_WORLD, &req_buff[0]);
        MPI_Irecv(&prev_val, 1, MPI_INT32_T, prev_id, tag, MPI_COMM_WORLD, &req_buff[1]);
    }

    if (next_id >= 0) {
        MPI_Isend(&chunk_arr[chunk_size - 1], 1, MPI_INT32_T, next_id, tag, MPI_COMM_WORLD, &req_buff[2]);
        MPI_Irecv(&next_val, 1, MPI_INT32_T, next_id, tag, MPI_COMM_WORLD, &req_buff[3]);
    }

    #else

    if (prev_id >= 0) MPI_Send(&chunk_arr[0], 1, MPI_INT32_T, prev_id, tag, MPI_COMM_WORLD);
    if (next_id >= 0) MPI_Send(&chunk_arr[chunk_size - 1], 1, MPI_INT32_T, next_id, tag, MPI_COMM_WORLD);

    if (prev_id >= 0) {
        MPI_Recv(&prev_val, 1, MPI_INT32_T, prev_id, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (prev_val > chunk_arr[0]) chunk_arr[0] = prev_val;
    }
    if (next_id >= 0) {
        MPI_Recv(&next_val, 1, MPI_INT32_T, next_id, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (chunk_arr[chunk_size - 1] > next_val) chunk_arr[chunk_size - 1] = next_val;
    }

    #endif
    
    // odd step
    for (int i = 1; i < BOUND(chunk_size); i += 2) {

        if (chunk_arr[i] > chunk_arr[i + 1]) {

            swap(chunk_arr, i, i + 1);
            sorted = 0;
        }
    }

    // even step
    for (int i = EVEN_START; i < BOUND(chunk_size); i += 2) {

        if (chunk_arr[i] > chunk_arr[i + 1]) {

            swap(chunk_arr, i, i + 1);
            sorted = 0;
        }
    }

    #ifndef BLOCKING
    if (prev_id >= 0) {

        for (int i = 0; i < 2; i++) MPI_Wait(&req_buff[i], MPI_STATUS_IGNORE);

        if (chunk_arr[0] < prev_val) {
            chunk_arr[0] = prev_val;
            sorted = 0;
        }

    }
    if (chunk_arr[0] > chunk_arr[1]) {

        swap(chunk_arr, 0, 1);
        sorted = 0;
    }

    if (next_id >= 0) {
        for (int i = 2; i < 4; i++) MPI_Wait(&req_buff[i], MPI_STATUS_IGNORE);

        if (chunk_arr[chunk_size - 1] > next_val) {
            chunk_arr[chunk_size - 1] = next_val;
            sorted = 0;
        }
    }
    if (chunk_arr[chunk_size - 2] > chunk_arr[chunk_size - 1]) {

        swap(chunk_arr, chunk_size - 2, chunk_size - 1);
        sorted = 0;
    }
    #endif

    return sorted;
}

#ifdef DEBUG
void print_chunk_arr(int id, int* chunk_arr, int chunk_size) {
    char *list = (char*)malloc(20 * chunk_size * sizeof(char));

    int j = 0;
    for (int i = 0; i < chunk_size; i++) {
        j += sprintf(&list[j], "%d, ", chunk_arr[i]);
    }

    fprintf(stdout, "%d) size=%d | %s\n", id, chunk_size, chunk_size > 0 ? list : "");
    fflush(stdout);

    free(list);
}
#endif
