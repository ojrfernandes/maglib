
#include "test.h"

int main() {
    char   source_path[] = "/home/jose/Software/maglib/maglit/our_m3dc1_data/i_coils/n03/C1.h5";
    char   shape_path[] = "tcabr_first_wall_m3dc1";
    char   output_path[] = "i_coils_n3.dat";
    double Rmin = 0.471;
    double Rmax = 0.478;
    int    nR = 10;
    int    nPhi = 20;

    // Check if the output file name is unique
    if (access(output_path, F_OK) != -1) {
        printf("There is a file with the same name saved to the chosen directory. Please change the output file name to avoid overwriting your data. \n");
        return 1;
    }

    tcabr_shape *shape = new tcabr_shape;
    if (!load_shape(shape_path, shape))
        return -1;

    maglit    tracer(source_path, FIO_M3DC1_SOURCE);
    auxfields aux_field(source_path, FIO_M3DC1_SOURCE);
    tracer.set_inside(shape, tcabr_inside);
    double dphi_init = 0.01;
    double dphi_min = 1e-5;
    double dphi_max = 0.2;
    tracer.configure(dphi_init, dphi_min, dphi_max);
    // tracer.set_verb(); // activate messages

    double Zfloor = -0.24;

    map_scalars scalars;
    FILE       *f0 = fopen(output_path, "w");
    printf("Running wall grid\n");
    fprintf(f0, "#R0 Z0 phi0 R1 Z1 phi1 deltaPhi length psiMin\n");
    for (int i = 0; i < nPhi / 3; i++) {
        double Z0 = Zfloor;
        double Z1;
        double phi0 = 2 * M_PI * i / nPhi;
        double phi1;
        for (int j = 0; j < nR; j++) {
            double R0 = Rmin + (Rmax - Rmin) * j / nR;
            double R1;
            map_wall(tracer, aux_field, R0, Z0, phi0, R1, Z1, phi1, scalars);
            fprintf(f0, "%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n", R0, Z0, phi0, R1, Z1, phi1, scalars.deltaPhi, scalars.length, scalars.psimin);
            printf("%f%%\n", 100.0 * (i * nR + j) / (nR * nPhi));
        }
    }

    fclose(f0);
    free_shape(shape);
    return 0;
}

void map_wall(maglit &tracer, auxfields &aux_field, double R0, double Z0, double phi0, double &R1, double &Z1, double &phi1, map_scalars &scalars) {
    int    status;
    double phi_max = 100 * 2 * M_PI;
    R1 = R0;
    Z1 = Z0;
    phi1 = phi0;
    tracer.reset();
    double arc = 0;
    double psin[1] = {5};
    scalars.psimin = *psin;
    do {
        R0 = R1;
        Z0 = Z1;
        phi0 = phi1;
        status = tracer.step(R1, Z1, phi1, phi_max, -1);
        if (status == SODE_CONTINUE_GOOD_STEP) {
            arc += dist(R0, Z0, phi0, R1, Z1, phi1);
            aux_field.psin_eval(R1, phi1, Z1, psin);
        }
        if (*psin < scalars.psimin)
            scalars.psimin = *psin;
    } while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP);
    scalars.deltaPhi = phi1 - phi0;
    scalars.length = arc;
    printf("\n%f\n", scalars.psimin);
    status_printer(status);
}

void status_printer(int status) {
    switch (status) {
    case SODE_OK:
        printf("SODE_OK\n");
        break;
    case SODE_SUCCESS_TIME:
        printf("SODE_SUCCESS_TIME\n");
        break;
    case SODE_SUCCESS_MONITOR:
        printf("SODE_SUCCESS_MONITOR\n");
        break;
    case SODE_CONTINUE_GOOD_STEP:
        printf("SODE_CONTINUE_GOOD_STEP\n");
        break;
    case SODE_CONTINUE_BAD_STEP:
        printf("SODE_CONTINUE_BAD_STEP\n");
        break;
    case SODE_FAILED:
        printf("SODE_FAILED\n");
        break;
    case SODE_BAD_FUNC:
        printf("SODE_BAD_FUNC\n");
        break;
    default:
        printf("Unknown Status\n");
        break;
    }
}

double
dist(double R0, double Z0, double phi0, double R1, double Z1, double phi1) {
    double dx = R1 * cos(phi1) - R0 * cos(phi0);
    double dy = R1 * sin(phi1) - R0 * sin(phi0);
    double dz = Z1 - Z0;
    return sqrt(dx * dx + dy * dy + dz * dz);
}