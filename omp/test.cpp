#include "../sode/sode.h"
#include <fusion_io.h>
#include <iostream>
#include <omp.h>

std::vector<fio_source *> src;
std::vector<fio_field *> mag_field;
fio_option_list opt;
sode solver(SODE_RK56_CK, 2);
std::vector<fio_hint> hint;

int main() {
    const char *source_path = "/home/jfernandes/Software/fusion-io/examples/data/m3dc1/C1.h5";
    const int source_type = FIO_M3DC1_SOURCE;
    const int timeslice = 1;
    const int nthreads = 10;

    // reserve nthreads spaces for vectors and initialize them
    src.reserve(nthreads);
    mag_field.reserve(nthreads);
    hint.reserve(nthreads);

    for (int i = 0; i < nthreads; i++) {
        // load source from file
        fio_open_source(&src[i], source_type, source_path);
        // set options for fields obtained from this source
        src[i]->get_field_options(&opt);
        opt.set_option(FIO_TIMESLICE, timeslice);
        opt.set_option(FIO_PART, FIO_TOTAL);
        // get magnetic field from source
        src[i]->get_field(FIO_MAGNETIC_FIELD, &mag_field[i], &opt);
    }

    double x[3] = {0.5, -0.2, 0.0};
    double B[10][3];

    omp_set_num_threads(nthreads);

#pragma omp parallel for
    for (int i = 0; i < 10; i++) {
        int j = omp_get_thread_num();
        src[j]->allocate_search_hint(&hint[j]);
        mag_field[j]->eval(x, B[i], hint[j]);
        src[j]->deallocate_search_hint(&hint[j]);
    }

    for (int i = 0; i < 10; i++) {
        std::cout << "(" << B[i][0] << ", " << B[i][1] << ", " << B[i][2] << ")" << std::endl;
    }

    return 0;
}
