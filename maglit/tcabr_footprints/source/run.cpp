#include "../../maglit.h"
#include "footprint.h"
#include "input_values.h"
#include "tcabr_collider.h"

int main() {
    // read params from input file
    std::string pathsFile = "params.txt";

    input_values input(pathsFile);
    bool readStatus = input.readInputFile();
    if (!readStatus) {
        std::cerr << "Error reading input file." << std::endl;
        return 1;
    }

    // load file containing vessel shape
    tcabr_shape *shape = new tcabr_shape;
    bool shapeStatus = load_shape(input.shape_path, shape);
    if (!shapeStatus) {
        std::cerr << "Error on loading vessel shape.\n"
                  << std::endl;
        return 1;
    }

    // create footprint object
    footprint footprint(input.plate, input.gridMin, input.gridMax, input.nGrid, input.nPhi);

    // set omp parameters
    omp_lock_t lock;
    omp_set_dynamic(1);
    omp_set_num_threads(input.num_theads);
    omp_init_lock(&lock);

#pragma omp parallel
    {
        // define tracer maglit object
        omp_set_lock(&lock); // lock to avoid concurrent hdf5 read
        maglit tracer(input.source_path, FIO_M3DC1_SOURCE, input.timeslice);
        omp_unset_lock(&lock); // unlock

        // configure tracer parameters
        tracer.set_inside(shape, tcabr_inside);
        tracer.configure(0.01, 1e-5, 0.2);

#pragma omp barrier

        footprint.runGrid(tracer);
    }
    omp_destroy_lock(&lock);

    // create output file
    std::cout << "\n Saving output file at " << input.output_path << std::endl;
    std::ofstream f0(input.output_path);
    if (!f0.is_open()) {
        std::cerr << "Failed to open file at " << input.output_path << std::endl;
        return 1;
    }

    // write output file header
    f0 << "#R0" << std::string(17, ' ') << "Z0" << std::string(17, ' ')
       << "phi0" << std::string(15, ' ') << "length" << std::string(13, ' ')
       << "psiMin\n";

    // Write output file data
    for (const auto &row : footprint.outputData) {
        f0 << std::fixed << std::setprecision(16)
           << row[0] << " " << row[1] << " " << row[2]
           << " " << row[3] << " " << row[4] << "\n";
    }

    footprint.outputData.shrink_to_fit();

    // close output file
    f0.close();

    // free allocated memory
    free_shape(shape);

    return 0;
}