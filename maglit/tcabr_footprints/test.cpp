#include "test.h"

class input_values {
  public:
    // Constructor to initialize reading_path
    input_values(const std::string &readingPath);
    // Destructor to clean up dynamically allocated memory
    ~input_values();
    // read paths to the source, shape and output from a text file along with the initial grid parameters
    bool readInputFile();

    // Variables to store the paths and parameters
    char *source_path;
    char *shape_path;
    char *output_path;
    int plate;
    int timeslice;
    double gridMin;
    double gridMax;
    int nGrid;
    int nPhi;

  private:
    std::string reading_path; // File path
};

// Constructor
input_values::input_values(const std::string &readingPath) : source_path(nullptr), shape_path(nullptr), output_path(nullptr), reading_path(readingPath) {}

// Destructor to free allocated memory
input_values::~input_values() {
    delete[] source_path;
    delete[] shape_path;
    delete[] output_path;
}

// read paths to the source, shape and output from a text file along with the initial grid parameters
bool input_values::readInputFile() {
    std::ifstream pathFile(this->reading_path);

    if (!pathFile.is_open()) {
        std::cerr << "Unable to open file: " << this->reading_path << std::endl;
        return false;
    }

    std::string line;
    int line_index = 0;

    while (std::getline(pathFile, line)) {
        if (!line.empty() && line[0] != '#') { // Ignore empty lines and comments starting with #
            switch (line_index) {
            case 0:
                // Allocate memory for source_path and copy the string
                source_path = new char[line.length() + 1];
                std::strcpy(source_path, line.c_str());
                break;
            case 1:
                // Allocate memory for shape_path and copy the string
                shape_path = new char[line.length() + 1];
                std::strcpy(shape_path, line.c_str());
                break;
            case 2:
                // Allocate memory for output_path and copy the string
                output_path = new char[line.length() + 1];
                std::strcpy(output_path, line.c_str());
                // Check if the output file already exists
                if (access(output_path, F_OK) != -1) {
                    std::cerr << "There is a file with the same name saved to the chosen directory. "
                              << "Please change the output file name to avoid overwriting your data."
                              << std::endl;
                    return false; // Exit the function to avoid overwriting
                }
                break;
            case 3:
                plate = std::stoi(line);
                break;
            case 4:
                timeslice = std::stoi(line);
                break;
            case 5:
                gridMin = std::stod(line);
                break;
            case 6:
                gridMax = std::stod(line);
                break;
            case 7:
                nGrid = std::stoi(line);
                break;
            case 8:
                nPhi = std::stoi(line);
                break;
            default:
                std::cerr << "Unexpected line index: " << line_index << std::endl;
                return false;
            }
            line_index++;
        }
    }

    if (line_index < 9) {
        std::cerr << "File does not contain enough data. Expected 9 values, found " << line_index << "." << std::endl;
        return false;
    }

    pathFile.close();
    return true;
}

int main() {
    // read params from input file
    std::string pathsFile = "params.txt";

    input_values input(pathsFile);
    bool rStatus = input.readInputFile();
    if (!rStatus) {
        std::cerr << "Error reading input file." << std::endl;
        return 1;
    }

    // load file containing vessel shape
    tcabr_shape *shape = new tcabr_shape;
    if (!load_shape(input.shape_path, shape)) {
        std::cerr << "Error on loading vessel shape.\n"
                  << std::endl;
        return 1;
    }

    // allocate memory to store results
    std::vector<std::string> dataWrite(input.nGrid * input.nPhi);

    int num_threads = omp_get_max_threads();
    omp_set_dynamic(1);
    omp_set_num_threads(num_threads);
    omp_lock_t lock;
    omp_init_lock(&lock);

#pragma omp parallel
    {
        // declare scalars type
        map_scalars scalars;
        // define tracer maglit object
        omp_set_lock(&lock);
        maglit tracer(input.source_path, FIO_M3DC1_SOURCE, input.timeslice);
        omp_unset_lock(&lock);

        // configure tracer parameters
        tracer.set_inside(shape, tcabr_inside);
        double dphi_init = 0.01;
        double dphi_min = 1e-5;
        double dphi_max = 0.2;
        tracer.configure(dphi_init, dphi_min, dphi_max);

        // run the grid
#pragma omp barrier

        if (input.plate == 0) {
            // vessel floor Z cordinate
            double Zfloor = -0.24;

            floor_grid(input.nPhi, input.nGrid, input.gridMin, input.gridMax, Zfloor, tracer, scalars, dataWrite);
        } else {
            // vessel wall R cordinate
            double Rfloor = 0.435;

            // invert map
            tracer.inverse_map(true);

            wall_grid(input.nPhi, input.nGrid, input.gridMin, input.gridMax, Rfloor, tracer, scalars, dataWrite);
        }
    }
    omp_destroy_lock(&lock);

    // create output file
    std::ofstream f0(input.output_path);
    if (!f0.is_open()) {
        std::cerr << "Failed to open file at" << input.output_path << std::endl;
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
    input.~input_values();

    return 0;
}

void floor_grid(int nPhi, int nR, double Rmin, double Rmax, double Zfloor, maglit &tracer, map_scalars &scalars, std::vector<std::string> &dataWrite) {
    // loop over the grid
#pragma omp for schedule(dynamic)
    for (int i = 0; i < nPhi; i++) {
        for (int j = 0; j < nR; j++) {
            double R0 = Rmin + (Rmax - Rmin) * j / nR;
            double phi0 = 2 * M_PI * i / nPhi;
            double Z0 = Zfloor;
            double R1, phi1, Z1;

            // evolve lines
            tracer.alloc_hint();
            evolve_lines(tracer, R0, Z0, phi0, R1, Z1, phi1, scalars);
            tracer.clear_hint();

            // Calculate the index for dataWrite
            int index = i * nR + j;

            // Create the output string
            std::ostringstream stringStream;
            stringStream << std::fixed << std::setprecision(3)
                         << R0 << " " << Z0 << " " << phi0 << " "
                         << R1 << " " << Z1 << " " << phi1 << " "
                         << scalars.deltaPhi << " " << scalars.length << " " << scalars.psimin << "\n";

            // Assign directly to the preallocated vector
            dataWrite[index] = stringStream.str();

            // Optional: Progress indication
            printf("Thread %d: progress = %f\n", omp_get_thread_num(), 100.0 * (i * nR + j) / (nR * nPhi));
        }
    }
}

void wall_grid(int nPhi, int nZ, double Zmin, double Zmax, double Rfloor, maglit &tracer, map_scalars &scalars, std::vector<std::string> &dataWrite) {
    // loop over the grid
#pragma omp for schedule(dynamic)
    for (int i = 0; i < nPhi; i++) {
        for (int j = 0; j < nZ; j++) {
            double R0 = Rfloor;
            double phi0 = 2 * M_PI * i / nPhi;
            double Z0 = Zmin + (Zmax - Zmin) * j / nZ;
            double R1, phi1, Z1;

            // evolve lines
            tracer.alloc_hint();
            evolve_lines(tracer, R0, Z0, phi0, R1, Z1, phi1, scalars);
            tracer.clear_hint();

            // Calculate the index for dataWrite
            int index = i * nZ + j;

            // Create the output string
            std::ostringstream stringStream;
            stringStream << std::fixed << std::setprecision(3)
                         << R0 << " " << Z0 << " " << phi0 << " "
                         << R1 << " " << Z1 << " " << phi1 << " "
                         << scalars.deltaPhi << " " << scalars.length << " " << scalars.psimin << "\n";

            // Assign directly to the preallocated vector
            dataWrite[index] = stringStream.str();

            // Optional: Progress indication
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
    double psin0 = 5.0;
    double *psin1 = &psin0;
    scalars.psimin = *psin1;
    tracer.reset();
    do {
        R0 = R1;
        Z0 = Z1;
        phi0 = phi1;
        status = tracer.step(R1, Z1, phi1, phi_max, -1);
        if (status == SODE_CONTINUE_GOOD_STEP) {
            arc += dist(R0, Z0, phi0, R1, Z1, phi1);
            tracer.psin_eval(R1, phi1, Z1, psin1);
            if (*psin1 < scalars.psimin) {
                scalars.psimin = *psin1;
            }
        }
    } while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP);
    scalars.deltaPhi = phi1 - phi0;
    scalars.length = arc;
}

// computes the connection_length
double dist(double R0, double Z0, double phi0, double R1, double Z1, double phi1) {
    double dx = R1 * cos(phi1) - R0 * cos(phi0);
    double dy = R1 * sin(phi1) - R0 * sin(phi0);
    double dz = Z1 - Z0;
    return sqrt(dx * dx + dy * dy + dz * dz);
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