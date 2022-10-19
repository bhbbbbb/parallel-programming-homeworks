#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <mpi.h>

// #define max (USHRT_max << 14)

#define UNIFORM() ((double)(rand() - 1) / RAND_MAX)

int main(int argc, char *argv[]) {

    typedef unsigned long long ull;

    int comm_size, id;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    srand(0xAAAAA + id);

    int n = 1;
    if (id == 0) {
        printf("Please enter an integer N to make each process toss in (0xFFFF << N) times: ");
        fflush(stdout);
        scanf("%d", &n);
    }
    double a_time = MPI_Wtime();
    for (int i = 1; i < comm_size; i <<= 1) {

        if (id < i && id + i < comm_size)
            MPI_Send(&n, 1, MPI_INT32_T, (id + i), 0, MPI_COMM_WORLD);
        
        else if (id < (i << 1)) {
            MPI_Recv(&n, 1, MPI_INT32_T, (id - i), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            fprintf(stdout, "%d) recieved n = %d in %f\n", id, n, (MPI_Wtime() - a_time));
            fflush(stdout);
        }
    }

    ull max = USHRT_MAX << n;


    ull num_hit = 0;
    double x, y;

    double start_time = MPI_Wtime();
    for (ull i = 0; i < max; i++) {
        x = UNIFORM();
        y = UNIFORM();
        num_hit += (int)((x * x + y * y) <= 1.0);
    }
    double spent_time = MPI_Wtime() - start_time;


    ull num_hit_recv;
    int reduction;
    for (int i = comm_size; i > 1 && id < i; ) {
        reduction = (i / 2) + (i & 1);
        if (id >= reduction)
            MPI_Send(&num_hit, 1, MPI_UINT64_T, id - reduction, 0, MPI_COMM_WORLD);
        else if (id + (i & 1) != reduction) {
            MPI_Recv(&num_hit_recv, 1, MPI_UINT64_T, id + reduction, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            num_hit += num_hit_recv;
            fflush(stdout);
        }
        i = reduction;
    }

    fprintf(stdout, "Process %d spent %f us in calculation.\n", id, spent_time * 1e6);
    fflush(stdout);

    if (id == 0) 
    {
        double pi = ((double)(4ULL * num_hit)) / ((ull)max * (ull)comm_size);
        printf("num_hit = %lld, pi = %.12f.\n", num_hit, pi);
    }
    MPI_Finalize();
    return 0;
}
