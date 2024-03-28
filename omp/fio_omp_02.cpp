#include <fusion_io.h>
#include <iostream>
#include <omp.h>
#include <sstream>

int main() {

    // Path to the source file
    const char *source_path = "/home/jfernandes/Software/maglib_local/maglit/our_m3dc1_data/i_coils/n03/C1.h5";

    double B[10][3];
    double x[10][3];

#pragma omp parallel
    {
        fio_source *src;
        int         result = fio_open_source(&src, FIO_M3DC1_SOURCE, source_path);
        if (result != FIO_SUCCESS) {
            std::cout << "Error opening file " << source_path << std::endl;
            delete src;
            exit(-1);
        }

        fio_option_list opt;
        src->get_field_options(&opt);
        opt.set_option(FIO_TIMESLICE, 0);
        opt.set_option(FIO_PART, FIO_TOTAL);

        fio_field *mag_field;
        result = src->get_field(FIO_MAGNETIC_FIELD, &mag_field, &opt);
        if (result != FIO_SUCCESS) {
            std::cout << "Error obtaining magnetic field from " << source_path << std::endl;
            delete src;
            exit(-1);
        }

#pragma omp for
        for (int i = 0; i < 10; i++) {
            x[i][0] = 0.6;
            x[i][1] = 0.3;
            x[i][2] = 0.1;
            result = (*mag_field).eval(x[i], B[i]);

            std::stringstream ss;
            ss << "BR, BPhi, BZ = " << B[i][0] << ", " << B[i][1] << ", " << B[i][2] << " from thread " << omp_get_thread_num() << std::endl;
            std::cout << ss.str();
        }

        delete src;
    }

    return 0;
}