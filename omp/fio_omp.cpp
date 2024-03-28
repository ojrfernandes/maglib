#include <fusion_io.h>
#include <iostream>
#include <omp.h>

int main() {

    fio_source *src;
    const char *source_path = "/home/jfernandes/Software/maglib_local/maglit/our_m3dc1_data/i_coils/n03/C1.h5";

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

    double B[10][3];
    double x[10][3];
    int    h_array[10];

    for (int i = 0; i < 10; i++) {
        x[i][0] = 0.6;
        x[i][1] = 0.4;
        x[i][2] = 0.1;
    }

    int i;

    size_t hint_size = (*src).sizeof_search_hint();
    std::cout << "Size of search hint: " << hint_size << std::endl;
    void *hint = calloc(10, hint_size);

    // #pragma omp parallel for
    for (i = 0; i < 10; i++) {
        void *h = ((char *)hint);
        result = (*mag_field).eval(x[i], B[i], h);
        h_array[i] = *((int *)h);
    }

    for (int i = 0; i < 10; i++) {
        std::cout << "\n"
                  << h_array[i] << std::endl;
        std::cout << "R, Phi, Z = " << x[i][0] << ", " << x[i][1] << ", " << x[i][2] << std::endl;
        std::cout << "BR, BPhi, BZ = " << B[i][0] << ", " << B[i][1] << ", " << B[i][2] << std::endl;
    }

    return 0;
}