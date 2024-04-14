#include "test.h"

int main() {
    // read paths input file
    std::string pathsFile = "paths.txt";
    char       *source_path;
    char       *shape_path;
    char       *output_path;
    readPaths(pathsFile, source_path, shape_path, output_path);

    // define grid parameters
    double Rmin = 0.471;
    double Rmax = 0.478;
    int    nR = 10;
    int    nPhi = 20;

    // check if the output file name is unique
    if (access(output_path, F_OK) != -1) {
        std::cerr << "There is a file with the same name saved to the chosen directory. Please change the output file name to avoid overwriting your data.\n"
                  << std::endl;
        return -1;
    }

    // load file containing vessel shape
    tcabr_shape *shape = new tcabr_shape;
    if (!load_shape(shape_path, shape)) {
        std::cerr << "Error on loading vessel shape.\n"
                  << std::endl;
        return -1;
    }

    // define tracer maglit and auxfields objects
    maglit tracer(source_path, FIO_M3DC1_SOURCE);

    // configure tracer parameters
    tracer.set_inside(shape, tcabr_inside);
    double dphi_init = 0.01;
    double dphi_min = 1e-5;
    double dphi_max = 0.2;
    tracer.configure(dphi_init, dphi_min, dphi_max);
    // tracer.set_verb(); // activate messages

    // vessel floor Z cordinate
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
            map_wall(tracer, R0, Z0, phi0, R1, Z1, phi1, scalars);
            fprintf(f0, "%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n", R0, Z0, phi0, R1, Z1, phi1, scalars.deltaPhi, scalars.length, scalars.psimin);
            printf("%f%%\n", 100.0 * (i * nR + j) / (nR * nPhi));
        }
    }

    fclose(f0);
    free_shape(shape);

    delete[] source_path;
    delete[] shape_path;
    delete[] output_path;

    return 0;
}

void map_wall(maglit &tracer, double R0, double Z0, double phi0, double &R1, double &Z1, double &phi1, map_scalars &scalars) {
    int    status;
    double phi_max = 100 * 2 * M_PI;
    R1 = R0;
    Z1 = Z0;
    phi1 = phi0;
    tracer.reset();
    double  arc = 0;
    double  psin_value = 5.0;
    double *psin = &psin_value;
    scalars.psimin = *psin;
    tracer.alloc_hint();
    do {
        R0 = R1;
        Z0 = Z1;
        phi0 = phi1;
        status = tracer.step(R1, Z1, phi1, phi_max, -1);
        if (status == SODE_CONTINUE_GOOD_STEP) {
            arc += dist(R0, Z0, phi0, R1, Z1, phi1);
            tracer.psin_eval(R1, phi1, Z1, psin);
        }
        if (*psin < scalars.psimin)
            scalars.psimin = *psin;
    } while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP);
    tracer.clear_hint();
    scalars.deltaPhi = phi1 - phi0;
    scalars.length = arc;
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

double dist(double R0, double Z0, double phi0, double R1, double Z1, double phi1) {
    double dx = R1 * cos(phi1) - R0 * cos(phi0);
    double dy = R1 * sin(phi1) - R0 * sin(phi0);
    double dz = Z1 - Z0;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

void readPaths(const std::string &readingPath, char *&source_path, char *&shape_path, char *&output_path) {

    std::ifstream pathFile(readingPath);

    if (pathFile.is_open()) {
        std::string source_path_str, shape_path_str, output_path_str;
        int         line_index = 0;
        std::string line;

        while (std::getline(pathFile, line) && line_index < 3) {
            if (!line.empty() && line[0] != '#') {
                if (line_index == 0) {
                    source_path_str = line;
                } else if (line_index == 1) {
                    shape_path_str = line;
                } else if (line_index == 2) {
                    output_path_str = line;
                }
                line_index++;
            }
        }
        pathFile.close();

        source_path = new char[source_path_str.length() + 1];
        std::strcpy(source_path, source_path_str.c_str());

        shape_path = new char[shape_path_str.length() + 1];
        std::strcpy(shape_path, shape_path_str.c_str());

        output_path = new char[output_path_str.length() + 1];
        std::strcpy(output_path, output_path_str.c_str());

    } else {
        std::cerr << "Unable to open file: " << readingPath << std::endl;
    }
}