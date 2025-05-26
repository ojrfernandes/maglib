#include "input_read.h"
#include "manifold.h"

int main() {

    // read params from input file
    std::string pathsFile = "mf_input.txt";

    // read input file
    std::cout << "\nReading input file..."
              << std::endl;
    input_read input(pathsFile);
    bool readStatus = input.readInputFile();
    if (!readStatus) {
        std::cerr << "Error reading input file." << std::endl;
        return 1;
    }

    // print input parameters
    std::cout << "\nInput parameters:" << std::endl;
    std::cout << "source_path: " << input.source_path << std::endl;
    std::cout << "output_path: " << input.output_path << std::endl;
    std::cout << "stability: " << input.stability << std::endl;
    std::cout << "timeslice: " << input.timeslice << std::endl;
    std::cout << "Phi: " << input.Phi << std::endl;
    std::cout << "epsilon: " << input.epsilon << std::endl;
    std::cout << "nSeg: " << input.nSeg << std::endl;
    std::cout << "l_lim: " << input.l_lim << std::endl;
    std::cout << "theta_lim: " << input.theta_lim << std::endl;

    if (input.R_xPoint == 0 && input.Z_xPoint == 0) {
        // read xnull and znull from the HDF5 file
        std::cout << "\nReading xnull and znull from the HDF5 file..."
                  << std::endl;
        bool readHDF5 = input.readHDF5File();
        if (!readHDF5) {
            std::cerr << "Error reading HDF5 file." << std::endl;
            return 1;
        }

        // print xnull and znull
        std::cout << "\nXnull and Znull read from HDF5 file:" << std::endl;
        std::cout << "xnull: " << input.R_xPoint << std::endl;
        std::cout << "znull: " << input.Z_xPoint << std::endl;
    } else {
        std::cout << "\nUsing xnull and znull from input file:" << std::endl;
        std::cout << "xnull: " << input.R_xPoint << std::endl;
        std::cout << "znull: " << input.Z_xPoint << std::endl;
    }

    // create manifold object
    std::cout << "\nCreating manifold object...\n"
              << std::endl;
    manifold manifold(input.source_path.c_str(), input.timeslice, input.Phi, input.stability, input.epsilon);

    // find the x-point
    std::cout << "\nFinding x-point..."
              << std::endl;
    bool found_xp = manifold.find_xPoint(input.R_xPoint, input.Z_xPoint);
    if (!found_xp) {
        std::cerr << "Failed to find the x-point" << std::endl;
        return 1;
    }

    // print x-point
    std::cout << std::fixed << std::setprecision(16) << "\n"
              << "x-point found at: " << "\n"
              << "R: " << manifold.xPoint.R << " Z: " << manifold.xPoint.Z << "\n"
              << std::endl;

    size_t num_points = 10;             // number of points in the primary segment
    std::vector<point> primary_segment; // vector of points to store the primary segment

    // compute primary segment
    std::cout << "\nComputing the first primary segment...\n"
              << std::endl;
    manifold.primarySegment(primary_segment, num_points);

    // save the primary segment to file input.output_path
    std::ofstream output_file(input.output_path);
    if (!output_file) {
        std::cerr << "Error opening output file: " << input.output_path << std::endl;
        return 1;
    }
    output_file << std::fixed << std::setprecision(16);
    for (const auto &pt : primary_segment) {
        output_file << pt.R << " " << pt.Z << "\n";
    }
    output_file.close();

    std::vector<point> new_segment; // vector of points to store the new segment

    // loop to create new segments
    for (int i = 1; i < input.nSeg; ++i) {
        // print progress bar
        manifold.progressBar(i, input.nSeg);
        std::this_thread::sleep_for(std::chrono::milliseconds(50));

        // create new segment from the primary segment
        manifold.newSegment(primary_segment, new_segment, input.Phi, i, input.l_lim, input.theta_lim);

        // append the new segment to file input.output_path
        std::ofstream output_file(input.output_path, std::ios::app);
        if (!output_file) {
            std::cerr << "Error opening output file: " << input.output_path << std::endl;
            return 1;
        }
        output_file << std::fixed << std::setprecision(16);
        for (const auto &pt : new_segment) {
            output_file << pt.R << " " << pt.Z << "\n";
        }
        output_file.close();

        primary_segment = new_segment;

        // empty new segment vector
        new_segment.clear();
    }

    // Done!
    std::cout << "\nDone! \n";
    std::cout << "Program finished successfully." << std::endl;

    return 0;
}