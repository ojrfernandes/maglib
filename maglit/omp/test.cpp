#include "../maglit.h"
#include "../tcabr_footprints/tcabr_collider.h"
#include <omp.h>

// computes the connection_length
double dist(double R0, double Z0, double phi0, double R1, double Z1, double phi1) {
    double dx = R1 * cos(phi1) - R0 * cos(phi0);
    double dy = R1 * sin(phi1) - R0 * sin(phi0);
    double dz = Z1 - Z0;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

omp_lock_t lock;

int main() {
    int num_threads = 5;
    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);
    char source_path[] = "/home/jfernandes/Software/maglib_local/maglit/our_m3dc1_data/i_coils/n03/C1.h5";
    char shape_path[] = "/home/jfernandes/Software/maglib/maglit/tcabr_footprints/tcabr_first_wall_m3dc1.txt";
    int timeslice = 1;

    // load file containing vessel shape
    tcabr_shape *shape = new tcabr_shape;
    if (!load_shape(shape_path, shape)) {
        std::cerr << "Error on loading vessel shape.\n"
                  << std::endl;
        return 1;
    }

    omp_init_lock(&lock);
    //  Parallelize the loop using OpenMP
#pragma omp parallel
    {
        //  each thread has its own tracer
        omp_set_lock(&lock);
        maglit tracer(source_path, FIO_M3DC1_SOURCE, timeslice);
        omp_unset_lock(&lock);

        // configure tracer parameters
        tracer.set_inside(shape, tcabr_inside);
        double dphi_init = 0.01;
        double dphi_min = 1e-5;
        double dphi_max = 0.2;
        tracer.configure(dphi_init, dphi_min, dphi_max);
#pragma omp barrier

// Evaluate the magnetic field in parallel
#pragma omp for
        for (int i = 0; i < 10; i++) {
            double R = 0.5;
            double Z = 0.1;
            double Phi = 0;
            double phi_max = 2 * M_PI;
            double arc = 0.0;
            int R0, Z0, phi0;
            int status = SODE_CONTINUE_GOOD_STEP;
            tracer.reset();
            tracer.alloc_hint(); // Hint is thread-specific now
            do {
                R0 = R;
                Z0 = Z;
                phi0 = Phi;
                status = tracer.step(R, Z, Phi, phi_max, 0);
                if (status == SODE_CONTINUE_GOOD_STEP) {
                    arc += dist(R0, Z0, phi0, R, Z, Phi);
                }
            } while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP);
            printf("Thread %d: arc length = %f\n", omp_get_thread_num(), arc);
            tracer.clear_hint();
        }
    }
    omp_destroy_lock(&lock);

    return 0;
}
