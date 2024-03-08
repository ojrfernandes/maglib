#include <omp.h>

void subdomain(float* x, int istart, int ipoints){
    int i;

    for (i = 0; i < ipoints; i++)
        x[istart+i] = 123.456;
}

void sub(float* x, int npoints){
    int iam, nt, ipoints, istart;

#pragma omp parallel default(shared) private(iam, nt, ipoints, istart)
    {
        iam = omp_get_thread_num();
        nt = omp_get_num_threads();
        ipoints = npoints / nt;
    }
}