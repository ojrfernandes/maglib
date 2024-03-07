#include <fusion_io.h>
// #include <omp.h>
#include <iostream>

int main() {
    fio_source *src;
    const char *source_path = "/home/jose/Software/maglib_backup/maglit/our_m3dc1_data/cp_coils/symmetric/n03/C1.h5";

    int result = fio_open_source(&src, FIO_M3DC1_SOURCE, source_path);
    if (result != FIO_SUCCESS) {
        std::cout << "Error opening file " << source_path << std::endl;
        delete (src);
        return -1;
    }

    fio_option_list opt;
    (*src).get_field_options(&opt);
    opt.set_option(FIO_TIMESLICE, 0);
    opt.set_option(FIO_PART, FIO_TOTAL);

    fio_field *mag_field;
    result = (*src).get_field(FIO_MAGNETIC_FIELD, &mag_field, &opt);
    if (result != FIO_SUCCESS) {
        std::cout << "Error opening file " << source_path << std::endl;
        mag_field = 0;
        return -1;
    }

    double B[3];
    double x[3];
    x[0] = 0.6; // R
    x[1] = 0.0; // Phi
    x[2] = 0.1; // Z
    result = (*mag_field).eval(x, B);
    std::cout << "BR, BPhi, BZ = " << B[0] << ", " << B[1] << ", " << B[2] << std::endl;
    // std::cout << "Printed by thread: " << omp_get_thread_num() << std::endl;

    return 0;
}