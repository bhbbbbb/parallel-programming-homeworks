#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cassert>
#include <memory>
#include <omp.h>
#include <climits>
// #include <mpi.h>

using WeightType = unsigned long;

constexpr int NUM_THREADS = 4;
constexpr int ALPHA = 1;
constexpr int BETA = 1;
constexpr int LAMBDA = 1;
constexpr int NUM_ANTS = 1 << 7;
constexpr int MU = 4;
constexpr int NUM_COLONIES = 1 << 18;
constexpr int QPOW = (sizeof(int) * 8 + 2);
constexpr WeightType Q = (1UL << QPOW);
constexpr WeightType QQ = Q >> 8;


class Matrix {

public:

    Matrix(const char* file_name) {

        std::ifstream fin(file_name);
        std::string fristline;

        std::getline(fin, fristline);
        std::stringstream ss(fristline);

        int tem;
        ss >> tem; // first int must be 0, so can be ignore

        while (ss >> tem) data.push_back(tem);

        n = data.size() + 1;

        data.resize(n * (n - 1) / 2); // avoid further copy
        probes.resize(n);
        probes[0] = 0;

        for (int i = 1; i < n - 1; i++) {

            probes[i] = probe(i, n);

            for (int j = 0; j < n; j++) {

                if (i >= j) fin >> tem; // ignore left-half

                else fin >> at(i, j);
            }

        }
        probes[n - 1] = probe(n - 1, n);

        fin.close();

        return;
    }

    Matrix(int _n, WeightType init_value): n{_n}, data(_n * (_n - 1) / 2, init_value), probes(_n) {

        for (int i = 0; i < n ; i++) probes[i] = probe(i, n);

    }

    inline WeightType& at(int i, int j) {
        if (i > j) std::swap(i, j);
        return data[probes[i] + j - i - 1];
    }

    inline const WeightType& at(int i, int j) const {
        if (i > j) std::swap(i, j);
        return data[probes[i] + j - i - 1];
    }

    inline auto begin() {
        return data.begin();
    }
    inline auto end() {
        return data.end();
    }
    inline auto cbegin() const {
        return data.cbegin();
    }
    inline auto cend() const {
        return data.cend();
    }

    inline int length() const {
        return n;
    }

    void print() const {

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) std::printf("%3ld ", i == j ? 0 : at(i, j));

            std::printf("\n");
        }
    }

private:

    static inline int probe(int i, int n) {
        return (2 * n - i - 1) * i / 2;
    }

    int n;
    std::vector<WeightType> data;
    std::vector<int> probes;
};

class Route {

public:

    Route(int num_cities): route(num_cities), length{0} {}

    inline const int& operator[](int i) const {
        return route[i];
    }
    inline int& operator[](int i) {
        return route[i];
    }

    inline bool better(const Route& rhs) const {
        return length < rhs.length;
    }

    void reset() {
        route.clear();
        route.resize(route.capacity());
        length = 0;
    }

    int update_length(const Matrix& dis_table) {
        
        length = dis_table.at(route[0], route[route.size() - 1]);

        for (int i = 0; i < static_cast<int>(route.size()) - 1; i++)
            length += dis_table.at(route[i], route[i + 1]);
        
        return length;
    }

    int get_length() const {
        return length;
    }

    void print() const {

        std::printf("L = %d\n", length);
        for (int c : route) std::printf("%2d, ", c);
        std::printf("\n");

    }

private:

    std::vector<int> route;
    int length;
};


class Pheromone : public Matrix {

public:

    Pheromone(int _n, WeightType init_value) : Matrix(_n, init_value) {}

    void update(const Route& route, const Matrix& dis_table) {

        at(route[0], route[length() - 1]) += q / dis_table.at(route[0], route[length() - 1]);

        for (int i = 0; i < length() - 1; i++)
            at(route[i], route[i + 1]) += q / dis_table.at(route[i], route[i + 1]);

    }

    void evaporate() {
        for (auto& t : *this) t = t * lambda / mu;
    }

    void dropout() {

        for (auto& t : *this) {
            if (t == 0) {
                if ((std::rand() & 0X3) == 0) t = QQ;
            }
            else if (std::rand() & 2 == 0) t << 3;
        }

    }
    
    inline static void set_params(int lambda, int mu, int q) {
        Pheromone::lambda = lambda;
        Pheromone::mu = mu;
        Pheromone::q = q;
    }

    void print() const {

        for (int i = 0; i < length(); i++) {
            
                
            for (int j = 0; j < length(); j++) {
                int count = 0;
                if (i != j) {
                    auto v = at(i, j);

                    while (v > 0) {
                        v >>= 1;
                        count++;
                    }
                }
                std::printf("%2d ", count);
            }

            std::printf("\n");
        }
    }

private:
    
    /**
     * (1 - rho) = lambda / mu, in which lambda <= mu;
     * tau <- (lambda / mu) * tau;
    */
    static int lambda, mu;
    static WeightType q;
};


class Ant {

public:

    Ant(int num_cities): route_size{0}, route(num_cities), visited(num_cities, false) {}

    const Route& traverse(const Pheromone& pheromone, const Matrix& dis_table) {

        // reset route
        if (route_size > 0) {
            route.reset();
            route_size = 0;
            visited.clear();
            visited.resize(visited.capacity());
        }


        // init initial city
        route[0] = std::rand() % dis_table.length();
        visited[route[route_size++]] = true;

        for (int i = 1; i < dis_table.length(); i++) {

            int next = random_select_next(pheromone, dis_table);

            visited[next] = true;

            route[route_size++] = next;
        }

        route.update_length(dis_table);

        return route;
    }

    inline const Route& get_route() const {
        return route;
    }

    inline static void set_params(int alpha, int beta) {
        Ant::alpha = alpha;
        Ant::beta = beta;
    }

private:

    int random_select_next(const Pheromone& pheromone, const Matrix& dis_table) const {

        auto cumulative_weights = get_residual_cities(pheromone, dis_table);

        WeightType max_value = 0;
        for (
            auto b = cumulative_weights.rbegin(), e = cumulative_weights.rend();
            b != e;
            b++
        )
            if (*b) {
                max_value = *b;
                break;
            }
        
        assert(max_value != 0);

        WeightType rand_val = std::rand() % max_value;


        int i;
        for (i = 0; i < dis_table.length(); i++) {

            if (rand_val < cumulative_weights[i]) break;
        }

        return i;
    }

    inline std::vector<WeightType> get_residual_cities(
        const Pheromone& pheromone, const Matrix& dis_table) const {

        assert(route_size > 0);

        std::vector<WeightType> cumulative_weights(dis_table.length(), 0);
        WeightType last = 0;

        for (int i = 0; i < dis_table.length(); i++) {

            if (visited[i]) continue;

            WeightType w = std::pow(pheromone.at(route[route_size - 1], i), alpha);

            w /= std::pow(dis_table.at(route[route_size - 1], i), beta);

            if (w <= 0) w = 1;

            last += w;
            cumulative_weights[i] = last;
        }

        return cumulative_weights;
    }


    int route_size;
    Route route;
    std::vector<bool> visited;
    static int alpha, beta;
};

int thread_task(
    std::vector<Ant>& ants,
    Pheromone& pheromone,
    const Matrix& dis_table,
    std::unique_ptr<Route>& best_route
) {


    const Route* local_best_route = nullptr;
    int found_best = 0;

    for (auto& ant : ants) {

        const auto& route = ant.traverse(pheromone, dis_table);

        if (!local_best_route) {
            local_best_route = &route;
            continue;
        }

        else if (route.better(*local_best_route)) {

            local_best_route = &route;
            found_best = route.get_length();
        }

    }

    if (best_route == nullptr || local_best_route->better(*best_route)) {
        #pragma omp critical (exchange)
        {
            if (best_route == nullptr || local_best_route->better(*best_route)) {
                best_route.reset(new Route(*local_best_route));
                best_route->print();
            }
        }
    }

    pheromone.evaporate();

    for (const auto& ant : ants) {
        
        pheromone.update(ant.get_route(), dis_table);

    }
    return found_best;
}

int Pheromone::lambda = LAMBDA;
int Pheromone::mu = MU;
WeightType Pheromone::q = Q;
int Ant::alpha = ALPHA;
int Ant::beta = BETA;

void parse_args(int argc, char** argv) {

    int num_threads = (argc >= 2) ? std::atoi(argv[1]) : NUM_THREADS;
    
    if (argc >= 4) Ant::set_params(std::atoi(argv[2]), std::atoi(argv[3]));

    else if (argc >= 7) Pheromone::set_params(
        std::atoi(argv[4]),
        std::atoi(argv[5]),
        std::atoi(argv[6])
    );

    omp_set_num_threads(num_threads);
}

/**
 * @param argv <___, num_processes, num_threads, alpha, beta, lambda, mu, q>
*/
int main(int argc, char** argv) {

    // int comm_size, p_id;
    // MPI_Init(&argc, &argv);
    // MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    // MPI_Comm_rank(MPI_COMM_WORLD, &p_id);

    parse_args(argc, argv);

    // Matrix dis_table("cities/gr17_d.txt");
    // Matrix dis_table("cities/fri26_d.txt");
    // Matrix dis_table("cities/dantzig42_d.txt");
    Matrix dis_table("cities/att48_d.txt");

    int num_ants = NUM_ANTS;
    int num_colonies = NUM_COLONIES;

    // this is the local of processes, global for their threads
    std::unique_ptr<Route> local_best_route;

    std::vector<Ant> ants(num_ants, Ant(dis_table.length()));
    Pheromone pheromone(dis_table.length(), 1);

    double start_time = omp_get_wtime();
    #pragma omp parallel default(none)\
        shared(local_best_route, num_colonies, dis_table) firstprivate(ants, pheromone)
    {
        int t_id = omp_get_thread_num();
        int current_local_best = INT32_MAX;
        // #pragma omp for
        for (int i = 0; i < num_colonies / omp_get_num_threads(); i++) {
            int found_best = thread_task(ants, pheromone, dis_table, local_best_route);
            if (found_best < current_local_best) current_local_best = found_best;
            if (i % 128 == 0 && current_local_best > local_best_route->get_length()) pheromone.dropout();
        }

        #pragma omp critical (print)
        {
            std::printf("pheronome at thread %d\n", omp_get_thread_num());
            pheromone.print();
            std::printf("---------------------------------------------------\n");
        }

        #pragma omp barrier
    }
    double spent_time = omp_get_wtime() - start_time;

    std::printf("num_threads = %d\n", NUM_THREADS);
    std::printf("num_ants = %d\n", num_ants);
    std::printf("num_colonies = %d\n", num_colonies);
    std::printf("spent_time = %lf\n", spent_time);
    std::printf("best route length: %d\n", local_best_route->get_length());
    local_best_route->print();
    local_best_route.release();
    std::printf("\n");

    // MPI_Finalize();

    return 0;
}
