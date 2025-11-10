#include "input_read.h"
#include "manifold.h"

int main() {

    // read params from input file
    std::string pathsFile = "mfgen_input.txt";

    std::cout << "\n-----------------------------------------------\n"
              << "MFGEN - Invariant Manifold Generator\n"
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
              << "output_path: " << input.output_path << "\n\n"
              << "--------------- TRACING PARAMETERS ------------\n\n"
              << "timeslice: " << input.timeslice << "\n"
              << "manifold: " << input.manifold << "\n"
              << "method: " << input.method << "\n"
              << "Phi: " << input.Phi << "\n\n"
              << "--------------- MULTIPLE POINCARE SECTIONS ----\n\n"
              << "nSections: " << input.nSections << "\n"
              << "phi_0: " << input.phi_0 << "\n"
              << "phi_1: " << input.phi_1 << "\n\n"
              << "--------------- ADDITIONAL PARAMETERS ---------\n\n"
              << "epsilon: " << input.epsilon << "\n"
              << "nSegments: " << input.nSegments << "\n"
              << "l_lim: " << input.l_lim << "\n"
              << "theta_lim: " << input.theta_lim << "\n"
              << "h_init: " << input.h_init << "\n"
              << "h_min: " << input.h_min << "\n"
              << "h_max: " << input.h_max << "\n"
              << "h_deriv: " << input.h_deriv << "\n"
              << "n_tol: " << input.n_tol << "\n"
              << "max_iter: " << input.max_iter << "\n"
              << "precision: " << input.precision << "\n"
              << "max_insertions: " << input.max_insertions << "\n\n"
              << "-----------------------------------------------" << std::endl;

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

    // create vector of input.nSections angular coordinates from input.phi_0 to input.phi_1 (deg)
    std::vector<double> poincare_sections;
    if (input.nSections > 1) {
        poincare_sections.push_back(input.phi_0);
        double dPhi = (input.phi_1 - input.phi_0) / (input.nSections - 1);
        std::cout << "\nPoincaré sections defined at Phi (deg) = " << input.phi_0;
        for (int i = 1; i < input.nSections; ++i) {
            double phi = input.phi_0 + i * dPhi;
            poincare_sections.push_back(phi);
            std::cout << ", " << phi;
        }
    } else {
        poincare_sections.push_back(input.Phi);
        std::cout << "\nPoincaré section defined at Phi (deg) = " << input.Phi << std::endl;
    }

    // loop over Poincaré sections
    for (const auto &phi : poincare_sections) {
        double phi_rad = phi * M_PI / 180.0; // convert degrees to radians

        std::cout << "\n\n"
                  << std::string(50, '-') << "\n"
                  << std::endl;

        std::cout << "Tracing section Phi = " << static_cast<int>(phi) << std::endl;

        // create maglit object
        std::cout << "\nCreating maglit object...\n"
                  << std::endl;
        maglit tracer(input.source_path.c_str(), FIO_M3DC1_SOURCE, input.timeslice);
        tracer.configure(input.h_init, input.h_min, input.h_max);

        // create manifold object
        std::cout << "\nCreating manifold object...\n"
                  << std::endl;
        manifold manifold(tracer, phi_rad, input.manifold);
        manifold.configure(input.epsilon, input.h_deriv, input.n_tol, input.max_iter, input.precision, input.max_insertions);

        // find the x-point
        std::cout << "\nFinding X-Point..."
                  << std::endl;
        bool found_xp = manifold.find_xPoint(input.R_xPoint, input.Z_xPoint);
        if (!found_xp) {
            std::cerr << "Failed to find the X-Point at Phi = " << phi << std::endl;
            continue;
        }

        // print x-point
        std::cout << std::fixed << std::setprecision(16) << "\n"
                  << "X-Point found at: " << "\n"
                  << "R: " << manifold.xPoint.R << " Z: " << manifold.xPoint.Z << "\n"
                  << std::endl;

        size_t             num_points = 10;       // number of points in the primary segment
        std::vector<point> first_primary_segment; // vector of points to store the first primary segment

        // compute first primary segment
        std::cout << "\nComputing the first primary segment...\n"
                  << std::endl;
        manifold.primarySegment(first_primary_segment, num_points);
        // if input.output_path ends with .dat, remove it
        if (input.output_path.size() >= 4 &&
            input.output_path.substr(input.output_path.size() - 4) == ".dat") {
            input.output_path = input.output_path.substr(0, input.output_path.size() - 4);
        }
        // concatenate input.output_path with the poincare sections
        std::string section_file = input.output_path + "_" + std::to_string(static_cast<int>(phi)) + ".dat";

        // save the primary segment to file section_file
        std::ofstream output_file(section_file);
        if (!output_file) {
            std::cerr << "Error opening output file: " << section_file << std::endl;
            return 1;
        }
        output_file << std::fixed << std::setprecision(16);
        for (const auto &pt : first_primary_segment) {
            output_file << pt.R << " " << pt.Z << "\n";
        }
        output_file.close();

        std::vector<point> new_segment; // vector of points to store the new segment

        // loop to create new segments
        for (int i = 1; i < input.nSegments; ++i) {
            // print progress bar
            manifold.progressBar(i, input.nSegments);
            std::this_thread::sleep_for(std::chrono::milliseconds(50));

            if (input.method == 0) {
                manifold.newSegment(first_primary_segment, new_segment, phi_rad, i, input.l_lim, input.theta_lim);
            } else if (input.method == 1) {
                manifold.newSegment(first_primary_segment, new_segment, phi_rad, input.l_lim, input.theta_lim);
                first_primary_segment = new_segment;
            } else {
                std::cerr << "Invalid method selected. Please choose 0 or 1." << std::endl;
                return 1;
            }

            // append the new segment to file section_file
            std::ofstream output_file(section_file, std::ios::app);
            if (!output_file) {
                std::cerr << "Error opening output file: " << section_file << std::endl;
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

        std::cout << "\nManifold successfully computed for Poincaré section phi = " << static_cast<int>(phi) << "\n"
                  << "Data saved to " << section_file << std::endl;
    }

    // Done!
    std::cout << "\nDone! \n";

    return 0;
}