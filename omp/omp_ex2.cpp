#include <iostream>
#include <omp.h>

// This is an example of a ~false sharing~ program

#define NUM_THREADS 4
static long num_steps = 100000;

int main() {

    int    i, nthreads;
    double step = 1.0 / (double)num_steps;
    double pi;
    double sum[NUM_THREADS];

    omp_set_num_threads(NUM_THREADS);

#pragma omp parallel
    {
        int    id, nthrds;
        double x;

        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();

        if (id == 0) // arbitrary choice of thread 0 to store the value of nthrds for when  it leaves the parallel region
            nthreads = nthrds;

        sum[id] = 0;
        for (int i = id; i < num_steps; i += nthrds) {
            x = (i + 0.5) * step;
            sum[id] += 4.0 / (1.0 + x * x);
        }
    }

    for (int i = 0; i < nthreads; i++) {
        pi += sum[i] * step;
    }

    std::cout << pi << std::endl;

    return 0;
}