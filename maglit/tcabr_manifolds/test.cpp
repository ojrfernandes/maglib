#include "../maglit.h"
#include "findx.h"

int main() {
    char source_path[] = "/home/jfernandes/Software/maglib_local/maglit/our_m3dc1_data/i_coils/n03/C1.h5";

    maglit tracer(source_path, FIO_M3DC1_SOURCE);
    double dphi_init = 0.01;
    double dphi_min = 1e-5;
    double dphi_max = 0.2;
    tracer.configure(dphi_init, dphi_min, dphi_max);

    double Phi0 = 0;
    double R_xpoint = 0.46;
    double Z_xpoint = -0.21;

    bool found = x_point(tracer, R_xpoint, Z_xpoint, Phi0, 1e-12, 100);
    if (found) {
        std::cout << "Found X-point at R = " << R_xpoint << ", Z = " << Z_xpoint << ", phi = " << Phi0 << std::endl;
    } else {
        std::cerr << "No X-point found" << std::endl;
        return 1;
    }

    int    crossings = 10;
    double phi_max = 2 * M_PI;

    FILE  *f1 = fopen("test_xpoint_2.dat", "a");
    int    status = SODE_CONTINUE_GOOD_STEP;
    double Phi = 0;
    for (int i = 0; i < crossings; i++) {
        tracer.reset();
        tracer.alloc_hint();
        Phi = Phi0;
        do {
            status = tracer.step(R_xpoint, Z_xpoint, Phi, phi_max, 0);
        } while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP);

        if (status == SODE_SUCCESS_TIME) {
            fprintf(f1, "%.16f  %.16f  %.16f\n", R_xpoint, Z_xpoint, Phi);
        } else {
            printf("SODE_ERROR\n");
        }
        Phi = Phi0;
        tracer.clear_hint();
    }

    fclose(f1);

    return 0;
}