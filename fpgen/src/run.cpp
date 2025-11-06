#include "footprint.h"
#include "input_read.h"

int main() {
    // read params from input file
    std::string pathsFile = "fpgen_input.txt";

    std::cout << "\n-----------------------------------------------\n"
              << "FPGEN - Magnetic Footprint Generator\n"
              << "-----------------------------------------------\n"
              << std::endl;
    input_read input(pathsFile);
    bool       readStatus = input.readInputFile();
    if (!readStatus) {
        std::cerr << "Error reading input file." << std::endl;
        return 1;
    }

    std::cout << "--------------- I/O FILES ---------------------\n\n"
              << "source_path: " << input.source_path << "\n"
              << "shape_path: " << input.shape_path << "\n"
              << "output_path: " << input.output_path << "\n\n"
              << "--------------- MAPPING PARAMETERS ------------\n\n"
              << "timeslice: " << input.timeslice << "\n"
              << "manifold: " << input.manifold << "\n"
              << "grid_R1: " << input.grid_R1 << "\n"
              << "grid_Z1: " << input.grid_Z1 << "\n"
              << "grid_R2: " << input.grid_R2 << "\n"
              << "grid_Z2: " << input.grid_Z2 << "\n"
              << "nRZ: " << input.nRZ << "\n"
              << "nPhi: " << input.nPhi << "\n\n"
              << "--------------- ADDITIONAL PARAMETERS ---------\n\n"
              << "num_threads: " << input.num_threads << "\n"
              << "max_turns: " << input.max_turns << "\n"
              << "h_init: " << input.h_init << "\n"
              << "h_min: " << input.h_min << "\n"
              << "h_max: " << input.h_max << "\n\n"
              << "-----------------------------------------------" << std::endl;

    // create footprint object
    std::cout << "\nCreating footprint object...\n\n"
              << "-----------------------------------------------"
              << std::endl;
    footprint footprint(input.manifold, input.grid_R1, input.grid_Z1, input.grid_R2, input.grid_Z2, input.nRZ, input.nPhi);

    // set omp parameters
    omp_set_num_threads(input.num_threads);

    if (input.num_threads > 1) {
        std::cout << "\nA maglit object must be created for each thread. \nYou shall see the same message printed a number of times."
                  << std::endl;
    } else {
        std::cout << "\nCreating maglit object..."
                  << std::endl;
    }

    std::vector<maglit> tracer;
    tracer.reserve(input.num_threads);
    for (int i = 0; i < input.num_threads; ++i) {
        // Construct a maglit object directly in the vector
        tracer.emplace_back(input.source_path.c_str(), FIO_M3DC1_SOURCE, input.timeslice);
        // Configure the newly created maglit object
        tracer.back().set_monitor(input.shape_path);
        tracer.back().configure(input.h_init, input.h_min, input.h_max);
    }

    std::cout << "\nMaglit object(s) created successfully."
              << std::endl;

    std::cout << "\nRunning grid...\n"
              << std::endl;

#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        footprint.runGrid(tracer[tid]);
    }

    // check if output data ends with .dat
    if (input.output_path.size() >= 4 &&
        input.output_path.substr(input.output_path.size() - 4) != ".dat") {
        input.output_path += ".dat";
    }
    // create output file
    std::cout << "\nSaving output file at " << input.output_path << std::endl;
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
    std::cout << "\nProgram finished successfully.\n"
              << std::endl;

    return 0;
}