#include <iostream>
#include <cstdlib>
#include <fstream>
#include <queue>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <omp.h>

#include <dirent.h>
#include <unistd.h>

constexpr int LINE_BUFF_MAX = 1024;
constexpr bool PRINT_RESULT = false;

class Queue {

public:
    Queue(): queue() {
        omp_init_lock(&lock);
    }

    void push(std::string& str) {
        omp_set_lock(&lock);
        queue.push(str);
        omp_unset_lock(&lock);
    }

    inline bool emtpy() const {
        return queue.empty();
    }

    inline std::size_t size() const {
        return queue.size();
    }

    std::string pop() {
        std::string str;
        #pragma omp critical (try_pop)
        {
            omp_set_lock(&lock);
            str = queue.front();
            queue.pop();
            omp_unset_lock(&lock);
        }
        return str;
    }


    ~Queue() {
        omp_destroy_lock(&lock);
    }

private:
    std::queue<std::string> queue;
    omp_lock_t lock;
};


std::vector<std::string> get_files_in_dir(const char* dirname) {
    constexpr auto EXT_TXT = ".txt";
    auto dir = opendir(dirname);
    std::vector<std::string> file_list;
    char file_path_buff[256];
    if (dir == NULL) perror("");

    for (auto ent = readdir(dir); ent != NULL; ent = readdir(dir)) {

        auto pch = strstr(ent->d_name, EXT_TXT);
        if (pch != NULL && strlen(pch) == 4) {
            strcpy(file_path_buff, dirname);
            strcat(file_path_buff, "/");
            strcat(file_path_buff, ent->d_name);
            file_list.emplace_back(file_path_buff);
        }
        // strncpy(ext, )
    }
    closedir(dir);
    return file_list;
}

void read_input_file(const std::string& file_name, Queue& queue) {

    std::ifstream fin(file_name);

    char buff[LINE_BUFF_MAX];
    while (fin.getline(buff, LINE_BUFF_MAX)) {
        
        auto str = std::string(buff);
        queue.push(str);
    }

    return;
}


void tokenize(Queue& queue, std::unordered_map<std::string, int>& buckets);

int main(int argc, char** argv) {

    const char* dir_name = argc < 3 ? "input" : argv[2];
    auto file_list = get_files_in_dir(dir_name);
    std::printf("found following files...\n");
    for (auto& str : file_list) std::cout << str << std::endl;

    int num_threads = argc < 2 ? omp_get_max_threads() : atoi(argv[1]);
    int num_files = file_list.size();
    int num_consumers = std::min(num_files, num_threads);

    std::printf("--------------------\n");
    std::printf("num_threads = %d, num_consumers = %d\n", num_threads, num_consumers);

    Queue queue;
    std::unordered_map<std::string, int> buckets;

    
    auto start_time = omp_get_wtime();
    double producer_time;
    #pragma omp parallel num_threads(num_threads)
    {

        // producer
        #pragma omp for
        for (int i = 0; i < num_consumers; i++) {
            read_input_file(file_list[i], queue);
        }


        #pragma omp single
        {
            producer_time = omp_get_wtime() - start_time;
            // std::printf("q_size = %ld\n", queue.size());
        }

        // consumer
        #pragma omp task
        tokenize(queue, buckets);
        
    }

    auto spent_time = omp_get_wtime() - start_time;

    std::printf("spent_time = %lf s, producer_time = %lf s\n", spent_time, producer_time);

    if (PRINT_RESULT)
        for (const auto& bucket : buckets)
            std::printf("(%s, %d)\n", bucket.first.c_str(), bucket.second);


    return 0;
}

void tokenize(Queue& queue, std::unordered_map<std::string, int>& buckets) {

    const std::string BLANK("++ ++");
    std::string line, str;
    std::stringstream ss;
    std::unordered_map<std::string, int> local_buckets;

    while (true) {
        
        #pragma omp critical (size)
        {
            if (!queue.emtpy()) {
                line = queue.pop();
            }
            else
                line = BLANK;
        }
        if (line == BLANK) break;

        ss.str(line);
        ss.clear();
        while (ss >> str) {
            local_buckets[str]++;
        }
    }

    #pragma omp critical (merge_buckets)
    for (const auto& bucket : local_buckets) buckets[bucket.first] += bucket.second;
}
