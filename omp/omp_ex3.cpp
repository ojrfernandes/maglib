#include <iostream>
#include <omp.h>

// In this version we solve the ~false sharing~ issue by adding an
// extra dimention to the 'sum' array. By padding the array we enssure
// each piece of information is stored in a different cache line

#define NUM_THREADS 8
#define PAD 8 // assume 64 byte L1 cache line size
static long num_steps = 100000;

int main() {

    int    i, nthreads;
    double step = 1.0 / (double)num_steps;
    double pi;
    double sum[NUM_THREADS][PAD];

    omp_set_num_threads(NUM_THREADS);

#pragma omp parallel
    {
        int    id, nthrds;
        double x;

        id = omp_get_thread_num();
        nthrds = omp_get_num_threads();

        if (id == 0) // arbitrary choice of thread 0 to store the value of nthrds for when  it leaves the parallel region
            nthreads = nthrds;

        sum[id][0] = 0;
        for (int i = id; i < num_steps; i += nthrds) {
            x = (i + 0.5) * step;
            sum[id][0] += 4.0 / (1.0 + x * x);
        }
    }

    for (int i = 0; i < nthreads; i++) {
        pi += sum[i][0] * step;
    }

    std::cout << pi << std::endl;

    return 0;
}