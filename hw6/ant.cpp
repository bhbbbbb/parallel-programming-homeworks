#include <vector>
#include <memory>
#include <omp.h>
#include <climits>
#include <mpi.h>

#include "ant.h"
#include "config.h"

int Pheromone::lambda = LAMBDA;
int Pheromone::mu = MU;
WeightType Pheromone::q = Q;
int Ant::alpha = ALPHA;
int Ant::beta = BETA;

void exange_local_best(std::unique_ptr<Route>& local_best_route, int p_id, int num_cities);
int parse_args(int p_id, int argc, char** argv);

/**
 * @param argv <___, num_processes, num_threads, alpha, beta, lambda, mu, q>
*/
int main(int argc, char** argv) {

    int comm_size, p_id;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &p_id);

    int num_threads = parse_args(p_id, argc, argv);

    Matrix dis_table("cities/gr17_d.txt");
    // Matrix dis_table("cities/fri26_d.txt");
    // Matrix dis_table("cities/dantzig42_d.txt");
    // Matrix dis_table("cities/att48_d.txt");

    int num_cities = dis_table.length();
    int num_ants = NUM_ANTS;
    int num_colonies = NUM_COLONIES;

    // this is the local of processes, global for their threads
    std::unique_ptr<Route> local_best_route;

    std::vector<Ant> ants(num_ants, Ant(num_cities));
    Pheromone pheromone(num_cities, 1);

    double start_time = omp_get_wtime();
    #pragma omp parallel default(none)\
        shared(local_best_route, num_colonies, dis_table, p_id, num_cities)\
        firstprivate(ants, pheromone)
    {
        int t_id = omp_get_thread_num();
        int current_local_best = INT32_MAX;

        for (int i = 0; i < num_colonies; i++) {
            int found_best = thread_task(ants, pheromone, dis_table, local_best_route);

            // if (found_best < current_local_best) current_local_best = found_best;
            // if (i % 128 == 0 && current_local_best > local_best_route->get_length()) pheromone.dropout();

            #pragma omp master
            if (i % GLOBAL_EXANGE_INTERVAL == 0 || i == num_colonies - 1) {
                exange_local_best(local_best_route, p_id, num_cities);
                MPI_Barrier(MPI_COMM_WORLD);
            }
        }

        #pragma omp critical (print)
        {
            // std::printf("pheronome at thread %d\n", omp_get_thread_num());
            // pheromone.print();
            // std::printf("---------------------------------------------------\n");
        }

        #pragma omp barrier
    }
    double spent_time = omp_get_wtime() - start_time;

    for (int i = 0; i < comm_size; i++) {

        if (p_id == i) {
            std::fprintf(stdout, "P%d) num_threads = %d\n", p_id, num_threads);
            std::fflush(stdout);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    // while (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {

    //     std::printf("not success in p%d", p_id);

    // }

    if (p_id == 0) {
        // std::printf("P%d) num_threads = %d\n", p_id, num_threads);
        std::printf("num_ants = %d\n", num_ants);
        std::printf("num_colonies_per_thread = %d\n", num_colonies);
        std::printf("alpha = %d, beta = %d\n", ALPHA, BETA);
        std::printf("lambda = %d, mu = %d, ", LAMBDA, MU);
        std::printf("Q = %ld\n", Q);
        std::printf("spent_time = %lf\n", spent_time);
        std::printf("best route length: \n");
        local_best_route->print();
        std::printf("\n");
    }

    MPI_Finalize();

    return 0;
}

void exange_local_best(std::unique_ptr<Route>& local_best_route, int p_id, int num_cities) {

    int local_send_buff[2] = {local_best_route->get_length(), p_id};
    int recv_buff[2];

    MPI_Allreduce(local_send_buff, recv_buff, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);

    MPI_Bcast(local_best_route->data(), num_cities, MPI_INT32_T, recv_buff[1], MPI_COMM_WORLD);

    if (p_id == recv_buff[1] && recv_buff[0] < local_best_route->get_length()) {

        local_best_route->update_length(recv_buff[0]);

        std::printf("P%d) new best route found: ", p_id);
        local_best_route->print();
        std::printf("\n");
    }

    return;
}

/**
 * @param argv <___, num_threads_for_p0, num_threads_for_p1, ...>
 * 
 * set -1 to automatically use max number of threads
 * unset to use DEFAULT_NUM_THREADS
 * 
*/
int parse_args(int p_id, int argc, char** argv) {

    if (argc >= p_id + 2) {

        int num_thread = std::atoi(argv[p_id + 1]);
        
        if (num_thread >= 0) {
            omp_set_num_threads(num_thread);
            
            return num_thread;
        }

        return omp_get_max_threads();

    }

    omp_set_num_threads(DEFAULT_NUM_THREADS);

    return DEFAULT_NUM_THREADS;
}
