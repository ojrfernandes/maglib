#include <fusion_io.h>
#include <m3dc1_source.h>
#include <iostream>
#include <omp.h>

int main() {

    fio_source *src;
    const char *source_path = "/home/jfernandes/Software/fusion-io/examples/data/m3dc1/C1.h5";

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

    double B[100][3];
    double x[100][3];
    int    h_array[100];

    for (int i = 0; i < 100; i++) {
        x[i][0] = 0.6;
        x[i][1] = 0.4;
        x[i][2] = 0.1;
    }   

    #pragma omp parallel for reduction(+:h_array[:100])
    for (int i = 0; i < 100; i++) {
        void* h;
        result = (*src).allocate_search_hint(&h);
        if (h == nullptr) {
            std::cerr << "Failed to allocate memory for hint." << std::endl;
        }
        else {
            result = (*mag_field).eval(x[i], B[i]);
            h_array[i] = *((int *)h);
            (*src).deallocate_search_hint(&h);
        }
    }
    
    for (int i = 0; i < 100; i++) {
        std::cout << "\n" << h_array[i] << std::endl;
        std::cout << "R, Phi, Z = " << x[i][0] << ", " << x[i][1] << ", " << x[i][2] << std::endl;
        std::cout << "BR, BPhi, BZ = " << B[i][0] << ", " << B[i][1] << ", " << B[i][2] << std::endl;
        }

    return 0;

}