#include <iostream>
#include <omp.h>

// In this version we solve the ~false sharing~ issue by using
// syncronization methods "critical"

#define NUM_THREADS 1
static long num_steps = 100000;

int main() {

    int    i, nthreads;
    double step = 1.0 / (double)num_steps;
    double pi;

    omp_set_num_threads(NUM_THREADS);

#pragma omp parallel
    {
        int    i, id, nthrds;
        double x, sum;

        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();

        if (id == 0) // arbitrary choice of thread 0 to store the value of nthrds for when  it leaves the parallel region
            nthreads = nthrds;

        for (i = id, sum = 0; i < num_steps; i += nthrds) {
            x = (i + 0.5) * step;
            sum += 4.0 / (1.0 + x * x);
        }
#pragma omp critical
        pi += sum * step;
    }

    std::cout << pi << std::endl;

    return 0;
}