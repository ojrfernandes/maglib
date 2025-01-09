#include "footprint.h"
#include "input_read.h"
#include "tcabr_collider.h"

int main() {
    // read params from input file
    std::string pathsFile = "params.txt";

    input_read input(pathsFile);
    bool readStatus = input.readInputFile();
    if (!readStatus) {
        std::cerr << "Error reading input file." << std::endl;
        return 1;
    }

    // initialize shape object
    tcabr_shape shape(input.shape_path);

    // create footprint object
    footprint footprint(input.plate, input.gridMin, input.gridMax, input.nGrid, input.nPhi);

    // set omp parametera
    omp_set_num_threads(input.num_theads);

    std::vector<maglit> tracer;
    tracer.reserve(input.num_theads);
    for (int i = 0; i < input.num_theads; ++i) {
        // Construct a maglit object directly in the vector
        tracer.emplace_back(input.source_path, FIO_M3DC1_SOURCE, input.timeslice);
        // Configure the newly created maglit object
        tracer.back().set_inside(shape, &tcabr_shape::tcabr_inside);
        tracer.back().configure(0.01, 1e-5, 0.2);
    }

#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        footprint.runGrid(tracer[tid]);
    }

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

    // print success message
    std::cout << "\nOutput file saved successfully.\n";
    std::cout << "Program finished successfully." << std::endl;

    return 0;
}