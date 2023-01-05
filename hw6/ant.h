#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cassert>
#include <memory>
#include <omp.h>


#include "config.h"

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

    inline int update_length(int new_length) {
        return length = new_length;
    }

    auto data() {
        return route.data();
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
