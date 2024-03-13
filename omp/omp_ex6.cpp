#include <iostream>
#include <omp.h>

static long num_steps = 100000;
double      step;

int main() {
    int    i;
    double x, pi, sum = 0.0;

    step = 1.0 / (double)num_steps;

#pragma omp parallel
    {
        double x;                  // It is fundamental that each thread has a copy of x
#pragma omp for reduction(+ : sum) // reduction clause used to parallelize incrementation
        for (i = 0; i < num_steps; i++) {
            x = (i + 0.5) * step;
            sum = sum + 4.0 / (1.0 + x * x);
        }

        pi = step * sum;
    }
    std::cout << pi << std::endl;

    return 0;
}