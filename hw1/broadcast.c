#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
    int comm_size, id;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    fprintf(stdout, "\nid = %d\n", id);
    fflush(stdout);

    int n = 1;
    if (id == 0) {
        // printf("Please enter an integer N to make each process toss in (0xFFFF << N) times: ");
        // fflush(stdout);
        // scanf("%d", &n);
        n = 2;
    }
    double start_time = MPI_Wtime();
    if (id == 0)
        for (int i = 1; i < comm_size; i++)
            MPI_Send(&n, 1, MPI_INT32_T, i, 0, MPI_COMM_WORLD);

    else {
        MPI_Recv(&n, 1, MPI_INT32_T, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        fprintf(stdout, "%d) recieved n = %d in %f\n", id, n, (MPI_Wtime() - start_time));
        fflush(stdout);
    }
}