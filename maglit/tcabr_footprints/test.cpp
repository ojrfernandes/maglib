#include "test.h"
#include <omp.h>

int main() {
    // read params from input file
    std::string pathsFile = "params.txt";
    char *source_path;
    char *shape_path;
    char *output_path;
    int plate;
    int timeslice;
    double gridMin;
    double gridMax;
    int nGrid;
    int nPhi;
    readParams(pathsFile, source_path, shape_path, output_path, plate, timeslice, gridMin, gridMax, nGrid, nPhi);

    // check if the output file name is unique
    if (access(output_path, F_OK) != -1) {
        std::cerr << "There is a file with the same name saved to the chosen directory. Please change the output file name to avoid overwriting your data.\n"
                  << std::endl;
        return 1;
    }

    // load file containing vessel shape
    tcabr_shape *shape = new tcabr_shape;
    if (!load_shape(shape_path, shape)) {
        std::cerr << "Error on loading vessel shape.\n"
                  << std::endl;
        return 1;
    }

    // declare scalars type
    map_scalars scalars;

    // allocate memory to store results
    std::vector<std::string> dataWrite;
    dataWrite.reserve(nGrid * nPhi);

    int num_threads = 8;
    omp_set_dynamic(0);
    omp_set_num_threads(num_threads);
    omp_lock_t lock;
    omp_init_lock(&lock);

#pragma omp parallel
    {
        // define tracer maglit object
        omp_set_lock(&lock);
        maglit tracer(source_path, FIO_M3DC1_SOURCE, timeslice);
        tracer.psin_init();
        omp_unset_lock(&lock);

        // configure tracer parameters
        tracer.set_inside(shape, tcabr_inside);
        double dphi_init = 0.01;
        double dphi_min = 1e-5;
        double dphi_max = 0.2;
        tracer.configure(dphi_init, dphi_min, dphi_max);

        std::vector<std::string> localDataWrite;
        dataWrite.reserve(nGrid);

        // run the grid
#pragma omp barrier

        if (plate == 0) {
            // print the parameters
            // std::cout << "Plate: Floor\n";
            // std::cout << "Rmin: " << gridMin << std::endl;
            // std::cout << "Rmax: " << gridMax << std::endl;
            // std::cout << "nR: " << nGrid << std::endl;
            // std::cout << "nPhi: " << nPhi << std::endl;

            // vessel floor Z cordinate
            double Zfloor = -0.24;

            floor_grid(nPhi, nGrid, gridMin, gridMax, Zfloor, tracer, scalars, dataWrite);
        } else {
            // std::cout << "Plate: Wall\n";
            // std::cout << "Zmin: " << gridMin << std::endl;
            // std::cout << "Zmax: " << gridMax << std::endl;
            // std::cout << "nZ: " << nGrid << std::endl;
            // std::cout << "nPhi: " << nPhi << std::endl;

            // vessel wall R cordinate
            double Rfloor = 0.435;

            // invert map
            tracer.inverse_map(true);

            wall_grid(nPhi, nGrid, gridMin, gridMax, Rfloor, tracer, scalars, localDataWrite);
        }
#pragma omp critical
        {
            dataWrite.insert(dataWrite.end(), localDataWrite.begin(), localDataWrite.end());
        }
    }
    omp_destroy_lock(&lock);

    // create output file
    std::ofstream f0(output_path);
    if (!f0.is_open()) {
        std::cerr << "Failed to open file at" << output_path << std::endl;
        return 1;
    }

    // write output file header
    f0 << "#R0 Z0 phi0 R1 Z1 phi1 deltaPhi length psiMin\n";
    // write output file data
    for (const auto &line : dataWrite) {
        f0 << line;
    }
    dataWrite.shrink_to_fit();

    // close output file
    f0.close();

    // free allocated memory
    free_shape(shape);
    delete[] source_path;
    delete[] shape_path;
    delete[] output_path;

    return 0;
}

void floor_grid(int nPhi, int nR, double Rmin, double Rmax, double Zfloor, maglit &tracer, map_scalars &scalars, std::vector<std::string> &dataWrite) {
    // loop over the grid
#pragma omp for
    for (int i = 0; i < nPhi; i++) {
        double Z0 = Zfloor;
        double Z1;
        double phi0 = 2 * M_PI * i / nPhi;
        double phi1;
        for (int j = 0; j < nR; j++) {
            double R0 = Rmin + (Rmax - Rmin) * j / nR;
            double R1;
            evolve_lines(tracer, R0, Z0, phi0, R1, Z1, phi1, scalars);
            std::ostringstream stringStream;
            stringStream << std::fixed << std::setprecision(3) << R0 << " " << Z0 << " " << phi0 << " " << R1 << " " << Z1 << " " << phi1 << " " << scalars.deltaPhi << " " << scalars.length << " " << scalars.psimin << "\n";
            dataWrite.push_back(stringStream.str());
            printf("Thread %d: progress = %f\n", omp_get_thread_num(), 100.0 * (i * nR + j) / (nR * nPhi));
        }
    }
}

void wall_grid(int nPhi, int nZ, double Zmin, double Zmax, double Rfloor, maglit &tracer, map_scalars &scalars, std::vector<std::string> &dataWrite) {
// loop over the grid
#pragma omp for
    for (int i = 0; i < nPhi; i++) {
        double R0 = Rfloor;
        double R1;
        double phi0 = 2 * M_PI * i / nPhi;
        double phi1;
        for (int j = 0; j < nZ; j++) {
            double Z0 = Zmin + (Zmax - Zmin) * j / nZ;
            double Z1;
            evolve_lines(tracer, R0, Z0, phi0, R1, Z1, phi1, scalars);
            std::ostringstream stringStream;
            stringStream << std::fixed << std::setprecision(3) << R0 << " " << Z0 << " " << phi0 << " " << R1 << " " << Z1 << " " << phi1 << " " << scalars.deltaPhi << " " << scalars.length << " " << scalars.psimin << "\n";
            dataWrite.push_back(stringStream.str());
            printf("Thread %d: progress = %f\n", omp_get_thread_num(), 100.0 * (i * nZ + j) / (nZ * nPhi));
        }
    }
}

// follow field lines and return scalar values when crossing the walls
void evolve_lines(maglit &tracer, double R0, double Z0, double phi0, double &R1, double &Z1, double &phi1, map_scalars &scalars) {
    R1 = R0;
    Z1 = Z0;
    phi1 = phi0;
    int status;
    double phi_max = 100 * 2 * M_PI;
    double arc = 0;
    double psin_value = 5.0;
    double *psin = &psin_value;
    scalars.psimin = *psin;
    tracer.reset();
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
    // status_printer(status);
}

// print the status of the integrator
void status_printer(int status) {
    switch (status) {
    case SODE_OK:
        std::cout << "SODE_OK\n";
        break;
    case SODE_SUCCESS_TIME:
        std::cout << "SODE_SUCCESS_TIME\n";
        break;
    case SODE_SUCCESS_MONITOR:
        std::cout << "SODE_SUCCESS_MONITOR\n";
        break;
    case SODE_CONTINUE_GOOD_STEP:
        std::cout << "SODE_CONTINUE_GOOD_STEP\n";
        break;
    case SODE_CONTINUE_BAD_STEP:
        std::cout << "SODE_CONTINUE_BAD_STEP\n";
        break;
    case SODE_FAILED:
        std::cout << "SODE_FAILED\n";
        break;
    case SODE_BAD_FUNC:
        std::cout << "SODE_BAD_FUNC\n";
        break;
    default:
        std::cout << "Unknown Status\n";
        break;
    }
}

// computes the connection_length
double dist(double R0, double Z0, double phi0, double R1, double Z1, double phi1) {
    double dx = R1 * cos(phi1) - R0 * cos(phi0);
    double dy = R1 * sin(phi1) - R0 * sin(phi0);
    double dz = Z1 - Z0;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

// read paths to the source, shape and output from a text file along with the initial variables grid parameters
void readParams(const std::string &readingPath, char *&source_path, char *&shape_path, char *&output_path, int &plate, int &timeslice, double &gridMin, double &gridMax, int &nGrid, int &nPhi) {

    std::ifstream pathFile(readingPath);

    if (pathFile.is_open()) {
        std::string source_path_str, shape_path_str, output_path_str, plate_str, timeslice_str, gridMin_str, gridMax_str, nGrid_str, nPhi_str;
        int line_index = 0;
        std::string line;

        while (std::getline(pathFile, line) && line_index < 9) {
            if (!line.empty() && line[0] != '#') {
                if (line_index == 0) {
                    source_path_str = line;
                } else if (line_index == 1) {
                    shape_path_str = line;
                } else if (line_index == 2) {
                    output_path_str = line;
                } else if (line_index == 3) {
                    plate_str = line;
                } else if (line_index == 4) {
                    timeslice_str = line;
                } else if (line_index == 5) {
                    gridMin_str = line;
                } else if (line_index == 6) {
                    gridMax_str = line;
                } else if (line_index == 7) {
                    nGrid_str = line;
                } else if (line_index == 8) {
                    nPhi_str = line;
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

        plate = std::stoi(plate_str);
        timeslice = std::stoi(timeslice_str);
        gridMin = std::stod(gridMin_str);
        gridMax = std::stod(gridMax_str);
        nGrid = std::stoi(nGrid_str);
        nPhi = std::stoi(nPhi_str);

    } else {
        std::cerr << "Unable to open file: " << readingPath << std::endl;
    }
}