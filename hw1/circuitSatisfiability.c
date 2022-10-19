/* circuitSatifiability.c solves the Circuit Satisfiability
 *  Problem using a brute-force sequential solution.
 *
 *   The particular circuit being tested is "wired" into the
 *   logic of function 'check_circuit'. All combinations of
 *   inputs that satisfy the circuit are printed.
 *
 *   16-bit version by Michael J. Quinn, Sept 2002.
 *   Extended to 32 bits by Joel C. Adams, Sept 2013.
 */

#include <stdio.h>     // printf()
#include <limits.h>    // UINT_MAX
#include <mpi.h>

int check_circuit(int, int);

int main(int argc, char *argv[]) {
    int count = 0;        /* number of solutions */
    int comm_size, id;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);



    fprintf(stdout, "Process %d/%d started.\n", id, comm_size);
    fflush(stdout);

    double start_time = MPI_Wtime();
    for (int i = id; i <= USHRT_MAX; i += comm_size) {
        count += check_circuit(id, i);
    }
    double spent_time = MPI_Wtime() - start_time;


    int count_recv;
    int reduction;
    for (int i = comm_size; i > 1 && id < i; ) {
        reduction = (i / 2) + (i & 1);
        if (id >= reduction)
            MPI_Send(&count, 1, MPI_INT32_T, id - reduction, 0, MPI_COMM_WORLD);
        else if (id + (i & 1) != reduction) {
            MPI_Recv(&count_recv, 1, MPI_INT32_T, id + reduction, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            count += count_recv;
            fflush(stdout);
        }
        i = reduction;
    }

    fprintf(stdout, "Process %d spent %f us in checking circuit.\n", id, spent_time * 1e6);
    fflush(stdout);

    if (id == 0) 
    {
        printf("A total of %d solutions were found.\n", count);
        printf("Process collected data in %f us.\n",(MPI_Wtime() - start_time - spent_time) * 1e6);
    }
    MPI_Finalize();
    return 0;
}

/* EXTRACT_BIT is a macro that extracts the ith bit of number n.
 *
 * parameters: n, a number;
 *             i, the position of the bit we want to know.
 *
 * return: 1 if 'i'th bit of 'n' is 1; 0 otherwise 
 */

#define EXTRACT_BIT(n,i) ( (n & (1<<i) ) ? 1 : 0)


/* check_circuit() checks the circuit for a given input.
 * parameters: id, the id of the process checking;
 *             bits, the (long) rep. of the input being checked.
 *
 * output: the binary rep. of bits if the circuit outputs 1
 * return: 1 if the circuit outputs 1; 0 otherwise.
 */

#define SIZE 16

int check_circuit (int id, int bits) {
    int v[SIZE];        /* Each element is a bit of bits */
    int i;

    for (i = 0; i < SIZE; i++)
        v[i] = EXTRACT_BIT(bits,i);

    if (  (v[0] || v[1]) && (!v[1] || !v[3]) && (v[2] || v[3])
        && (!v[3] || !v[4]) && (v[4] || !v[5])
        && (v[5] || !v[6]) && (v[5] || v[6])
        && (v[6] || !v[15]) && (v[7] || !v[8])
        && (!v[7] || !v[13]) && (v[8] || v[9])
        && (v[8] || !v[9]) && (!v[9] || !v[10])
        && (v[9] || v[11]) && (v[10] || v[11])
        && (v[12] || v[13]) && (v[13] || !v[14])
        && (v[14] || v[15])  )
    {
        printf ("%d) %d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d \n", id,
            v[15],v[14],v[13],v[12],
            v[11],v[10],v[9],v[8],v[7],v[6],v[5],v[4],v[3],v[2],v[1],v[0]);
        fflush (stdout);
        return 1;
    } else {
        return 0;
    }
}

