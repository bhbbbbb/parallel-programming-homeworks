#include <iostream>
#include <string>
#include <vector>
#include <pthread.h>

int thread_count;

void* hello(void* rank);

int main(int argc, char** argv) {

    // pthread_t* thread_handles;

    thread_count = argc < 2 ? 1 : std::atoi(argv[1]);

    std::cout << "thread_count = " << thread_count << std::endl;

    std::vector<pthread_t> thread_handles(thread_count);
    // thread_handles = new pthread_t[thread_count];

    for (long t = 0; t < thread_count; t++)
        pthread_create(&thread_handles[t], NULL, hello, (void*)t);
    
    std::cout << "Hello from the main thread\n";

    for (long t = 0; t < thread_count; t++)
        pthread_join(thread_handles[t], NULL);
    
    // delete[] thread_handles;
}

void* hello(void* rank) {

    auto my_rank = (long)rank;

    std::printf("hello from %ld of %d\n", my_rank, thread_count);

    return NULL;
}