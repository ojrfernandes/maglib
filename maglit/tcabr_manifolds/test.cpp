#include "../maglit.h"
#include "findx.h"

int main() {
    char source_path[] = "/home/jfernandes/Software/maglib_local/maglit/our_m3dc1_data/i_coils/n03/C1.h5";

    maglit tracer(source_path, FIO_M3DC1_SOURCE);
    double dphi_init = 0.01;
    double dphi_min = 1e-5;
    double dphi_max = 0.2;
    tracer.configure(dphi_init, dphi_min, dphi_max);

    double Phi0 = 0.2 * M_PI;
    double R_xpoint = 0.48;
    double Z_xpoint = -0.2;

    int    crossings = 100;
    double R0 = 0.4698660634246061;
    double Z0 = -0.2134422706264537;
    double phi_max = 2 * M_PI;

    FILE  *f1 = fopen("crossings.dat", "a");
    int    status = SODE_CONTINUE_GOOD_STEP;
    double Phi = 0;
    for (int i = 0; i < crossings; i++) {
        tracer.reset();
        tracer.alloc_hint();
        Phi = Phi0;
        do {
            status = tracer.step(R0, Z0, Phi, phi_max, 0);
            // printf("status: %d\n", status);
        } while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP);
        Phi = Phi0;
        if (status == SODE_SUCCESS_TIME) {
            fprintf(f1, "%f  %f  %f\n", R0, Z0, Phi0);
            printf("SODE_SUCCESS_TIME");
        } else {
            printf("SODE_ERROR");
        }

        tracer.clear_hint();
    }
    fclose(f1);

    bool found = newton_raphson(tracer, R_xpoint, Z_xpoint, Phi0);
    if (found) {
        std::cout << "Found X-point at R = " << R_xpoint << ", Z = " << Z_xpoint << ", phi = " << Phi0 << std::endl;
    } else {
        std::cout << "No X-point found" << std::endl;
    }

    return 0;
}