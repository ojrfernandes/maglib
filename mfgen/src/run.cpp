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
    std::cout << "method: " << input.method << std::endl;
    std::cout << "slices: " << input.slices << std::endl;
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

    // create vector of input.slices points from input.phi_0 to input.phi_1
    std::vector<double> phi_slices = {input.phi_0};
    if (input.slices > 1) {
        std::cout << "\nCreating vector of points from phi_0 to phi_1..." << std::endl;
        double dPhi = (input.phi_1 - input.phi_0) / (input.slices - 1);
        std::cout << "Phi = " << input.phi_0;
        for (int i = 1; i < input.slices; ++i) {
            double phi = input.phi_0 + i * dPhi;
            phi_slices.push_back(phi);
            std::cout << ", " << phi;
        }
    } else {
        std::cout << "\nUsing single point with Phi = " << input.Phi << std::endl;
    }

    // loop over phi slices
    for (const auto &phi : phi_slices) {

        std::cout << "\n\n"
                  << std::string(50, '#') << "\n"
                  << std::endl;

        std::cout << "Processing phi = " << static_cast<int>(phi) << std::endl;
        // create manifold object
        std::cout << "\nCreating manifold object...\n"
                  << std::endl;
        manifold manifold(input.source_path.c_str(), input.timeslice, phi, input.stability, input.epsilon);

        // find the x-point
        std::cout << "\nFinding x-point..."
                  << std::endl;
        bool found_xp = manifold.find_xPoint(input.R_xPoint, input.Z_xPoint);
        if (!found_xp) {
            std::cerr << "Failed to find the x-point at Phi = " << phi << std::endl;
            continue;
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

        // concatenate input.output_path with the slice
        std::string slice_file = input.output_path + "_" + std::to_string(static_cast<int>(phi)) + ".dat";

        // save the primary segment to file slice_file
        std::ofstream output_file(slice_file);
        if (!output_file) {
            std::cerr << "Error opening output file: " << slice_file << std::endl;
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

            if (input.method == 0) {
                manifold.newSegment(primary_segment, new_segment, phi, i, input.l_lim, input.theta_lim);
            } else if (input.method == 1) {
                manifold.newSegment(primary_segment, new_segment, phi, input.l_lim, input.theta_lim);
                primary_segment = new_segment;
            } else {
                std::cerr << "Invalid method selected. Please choose 0 or 1." << std::endl;
                return 1;
            }

            // append the new segment to file slice_file
            std::ofstream output_file(slice_file, std::ios::app);
            if (!output_file) {
                std::cerr << "Error opening output file: " << slice_file << std::endl;
                return 1;
            }
            output_file << std::fixed << std::setprecision(16);
            for (const auto &pt : new_segment) {
                output_file << pt.R << " " << pt.Z << "\n";
            }
            output_file.close();

            // empty new segment vector
            new_segment.clear();
        }

        std::cout << "\nManifold successfully computed for slice Phi = " << static_cast<int>(phi) << "\n"
                  << "Data saved to " << slice_file << std::endl;
    }

    // Done!
    std::cout << "\nDone! \n";

    return 0;
}