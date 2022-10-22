#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <mpi.h>

#define MAX (USHRT_MAX << 14)

#define UNIFORM() ((double)(rand() - 1) / RAND_MAX)

int main(int argc, char *argv[]) {

    typedef unsigned long long ull;

    int comm_size, id;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    srand(0xAAAAA + id);

    ull num_hit = 0;
    double x, y;

    double start_time = MPI_Wtime();
    for (ull i = 0; i < MAX; i++) {
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
            MPI_Send(&num_hit, 1, MPI_INT32_T, id - reduction, 0, MPI_COMM_WORLD);
        else if (id + (i & 1) != reduction) {
            MPI_Recv(&num_hit_recv, 1, MPI_INT32_T, id + reduction, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            num_hit += num_hit_recv;
            fflush(stdout);
        }
        i = reduction;
    }

    fprintf(stdout, "Process %d spent %f us in checking circuit.\n", id, spent_time * 1e6);
    fflush(stdout);

    if (id == 0) 
    {
        double pi = ((double)(4ULL * num_hit)) / ((ull)MAX * (ull)comm_size);
        printf("num_hit = %lld, pi = %.12f.\n", num_hit, pi);
        printf("Process collected data in %f us.\n",(MPI_Wtime() - start_time - spent_time) * 1e6);
    }
    MPI_Finalize();
    return 0;
}
