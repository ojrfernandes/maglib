#include "input_read.h"

// Constructor
input_read::input_read(const std::string &readingPath) : reading_path(readingPath) {}

bool input_read::readInputFile() {
    std::ifstream file(this->reading_path);
    std::string line;

    // Open the file and check for errors
    if (!file) {
        std::cerr << "Error reading input: Unable to open file: " << reading_path << std::endl;
        return false;
    }

    // check if the file is empty
    if (file.peek() == std::ifstream::traits_type::eof()) {
        std::cerr << "Error reading input: Input file is empty: " << reading_path << std::endl;
        return false;
    }

    // Lambda function to trim whitespace from both ends of a string
    auto trim = [](std::string &s) {
        size_t start = s.find_first_not_of(" \t");
        size_t end = s.find_last_not_of(" \t");
        if (start == std::string::npos) {
            s.clear();
        } else {
            s = s.substr(start, end - start + 1);
        }
    };

    while (std::getline(file, line)) {
        // Trim the line first
        trim(line);

        // Ignore empty lines and lines that start with '#'
        if (line.empty() || line[0] == '#') {
            continue;
        }

        // Find '=' position
        auto eq_pos = line.find('=');
        if (eq_pos == std::string::npos) {
            std::cerr << "Error reading input: Invalid line (missing '='): " << line << std::endl;
            return false;
        }

        // Split into key and value
        std::string key = line.substr(0, eq_pos);
        std::string value = line.substr(eq_pos + 1);

        // Handle inline comments in value
        auto comment_pos = value.find('#');
        if (comment_pos != std::string::npos) {
            value = value.substr(0, comment_pos);
        }

        // Trim both key and value
        trim(key);
        trim(value);

        if (key.empty()) {
            std::cerr << "Error reading input: Empty key in line: " << line << std::endl;
            return false;
        }

        // Assign values to corresponding variables
        if (key == "source_path") {
            this->source_path = value;
        } else if (key == "output_path") {
            this->output_path = value;
            // check if there is a saved file with the same name
            std::ifstream f0(this->output_path);
            if (f0.is_open()) {
                std::cerr << "Warning: File " << this->output_path << " already exists. Please, change your output file name to avoid overwriting your data." << std::endl;
                return false;
            }
            f0.close();
        } else if (key == "timeslice") {
            this->timeslice = std::stoi(value);
            // check if timeslice is -1 or a non-negative integer
            if (this->timeslice < -1) {
                std::cerr << "Error: Timeslice must be -1 (equilibrium) or a non-negative integer." << std::endl;
                return false;
            }
        } else if (key == "manifold") {
            this->manifold = std::stoi(value);
            // check if manifold is 0 or 1
            if (this->manifold != 0 && this->manifold != 1) {
                std::cerr << "Error: Manifold must be 0 (unstable) or 1 (stable)." << std::endl;
                return false;
            }
        } else if (key == "method") {
            this->method = std::stoi(value);
            // check if method is 0 or 1
            if (this->method != 0 && this->method != 1) {
                std::cerr << "Error: Method must be 0 (exact) or 1 (interpolant)." << std::endl;
                return false;
            }
        } else if (key == "Phi") {
            this->Phi = std::stod(value);
        } else if (key == "nSections") {
            this->nSections = std::stoi(value);
            // check if nSections is a positive integer
            if (this->nSections < 1) {
                std::cerr << "Error: nSections must be a positive integer." << std::endl;
                return false;
            }
        } else if (key == "phi_0") {
            this->phi_0 = std::stod(value);
        } else if (key == "phi_1") {
            this->phi_1 = std::stod(value);
        } else if (key == "epsilon") {
            this->epsilon = std::stod(value);
            // check if epsilon is positive
            if (this->epsilon <= 0) {
                std::cerr << "Error: epsilon must be a positive value." << std::endl;
                return false;
            }
        } else if (key == "nSegments") {
            this->nSegments = std::stoi(value);
            // check if nSegments is a positive integer
            if (this->nSegments < 1) {
                std::cerr << "Error: nSegments must be a positive integer." << std::endl;
                return false;
            }
        } else if (key == "l_lim") {
            this->l_lim = std::stod(value);
            // check if l_lim is positive
            if (this->l_lim <= 0) {
                std::cerr << "Error: l_lim must be a positive number." << std::endl;
                return false;
            }
        } else if (key == "theta_lim") {
            this->theta_lim = std::stod(value);
            // check if theta_lim is positive
            if (this->theta_lim <= 0) {
                std::cerr << "Error: theta_lim must be a positive number." << std::endl;
                return false;
            }
        } else if (key == "h_init") {
            this->h_init = std::stod(value);
            // check if h_init is positive
            if (this->h_init <= 0) {
                std::cerr << "Error: h_init must be a positive number." << std::endl;
                return false;
            }
        } else if (key == "h_min") {
            this->h_min = std::stod(value);
            // check if h_min is positive
            if (this->h_min <= 0) {
                std::cerr << "Error: h_min must be a positive number." << std::endl;
                return false;
            }
        } else if (key == "h_max") {
            this->h_max = std::stod(value);
            if (this->h_max <= 0) {
                std::cerr << "Error: h_max must be a positive number." << std::endl;
                return false;
            }
        } else if (key == "h_deriv") {
            this->h_deriv = std::stod(value);
            // check if h_deriv is positive
            if (this->h_deriv <= 0) {
                std::cerr << "Error: h_deriv must be a positive number." << std::endl;
                return false;
            }
        } else if (key == "n_tol") {
            this->n_tol = std::stod(value);
            // check if n_tol is positive
            if (this->n_tol <= 0) {
                std::cerr << "Error: n_tol must be a positive number." << std::endl;
                return false;
            }
        } else if (key == "max_iter") {
            this->max_iter = std::stoi(value);
            // check if max_iter is a positive integer
            if (this->max_iter < 1) {
                std::cerr << "Error: max_iter must be a positive integer." << std::endl;
                return false;
            }
        } else if (key == "precision") {
            this->precision = std::stod(value);
            // check if precision is positive
            if (this->precision <= 0) {
                std::cerr << "Error: precision must be a positive number." << std::endl;
                return false;
            }
        } else if (key == "max_insertions") {
            this->max_insertions = std::stoi(value);
            // check if max_insertions is a positive integer
            if (this->max_insertions < 1) {
                std::cerr << "Error: max_insertions must be a positive integer." << std::endl;
                return false;
            }
        } else if (key == "verbose") {
            this->verbose = std::stoi(value);
            // check if verbose is 0 or 1
            if (this->verbose != 0 && this->verbose != 1) {
                std::cerr << "Error: verbose must be 0 (no) or 1 (yes)." << std::endl;
                return false;
            }
        } else if (key == "R_xPoint") {
            this->R_xPoint  = std::stod(value);
            this->xpoint_set = true;
        } else if (key == "Z_xPoint") {
            this->Z_xPoint  = std::stod(value);
            this->xpoint_set = true;
        } else if (key == "xpoint_index") {
            this->xpoint_index = std::stoi(value);
            if (this->xpoint_index < -1) {
                std::cerr << "Error: xpoint_index must be -1 (auto) or a non-negative integer." << std::endl;
                return false;
            }
        }

        else {
            std::cerr << "Error: Invalid key: " << key << std::endl;
            return false;
        }
    }

    // Check that required parameters are set
    if (source_path.empty() || output_path.empty()) {
        std::cerr << "Error reading input: Missing required parameters." << std::endl;
        return false;
    }

    // Validate step-size relationships after all keys are parsed
    if (this->h_max < this->h_min) {
        std::cerr << "Error: h_max must be greater than or equal to h_min." << std::endl;
        return false;
    }

    return true;
}

// read xnull and znull from the HDF5 file
bool input_read::readHDF5File() {
    // paths to the xnull and znull datasets in the HDF5 file
    const std::string xnullPath = "scalars/xnull";
    const std::string znullPath = "scalars/znull";

    try {
        hid_t file = H5Fopen(source_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file < 0) {
            std::cerr << "Error: could not open HDF5 file: " << source_path << std::endl;
            return false;
        }

        hid_t xnullDataset = H5Dopen(file, "scalars/xnull", H5P_DEFAULT);
        hid_t znullDataset = H5Dopen(file, "scalars/znull", H5P_DEFAULT);
        if (xnullDataset < 0 || znullDataset < 0) {
            std::cerr << "Error: could not open scalars/xnull or scalars/znull in HDF5 file." << std::endl;
            H5Fclose(file);
            return false;
        }

        hid_t xnullSpace = H5Dget_space(xnullDataset);
        hid_t znullSpace = H5Dget_space(znullDataset);

        hsize_t xnullSize;
        hsize_t znullSize;
        H5Sget_simple_extent_dims(xnullSpace, &xnullSize, nullptr);
        H5Sget_simple_extent_dims(znullSpace, &znullSize, nullptr);

        if (xnullSize == 0 || znullSize == 0) {
            std::cerr << "Error: xnull or znull dataset is empty in HDF5 file." << std::endl;
            H5Sclose(xnullSpace);
            H5Sclose(znullSpace);
            H5Dclose(xnullDataset);
            H5Dclose(znullDataset);
            H5Fclose(file);
            return false;
        }

        std::vector<double> xnullData(xnullSize);
        std::vector<double> znullData(znullSize);

        H5Dread(xnullDataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xnullData.data());
        H5Dread(znullDataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, znullData.data());

        // Deduplicate: collect indices of distinct X-points (exact equality).
        std::vector<hsize_t> unique_idx;
        for (hsize_t k = 0; k < xnullSize; ++k) {
            bool found = false;
            for (hsize_t u : unique_idx)
                if (xnullData[k] == xnullData[u] && znullData[k] == znullData[u]) { found = true; break; }
            if (!found) unique_idx.push_back(k);
        }

        if (unique_idx.size() == 1) {
            // One distinct X-point (possibly stored multiple times) — use it directly.
            R_xPoint = xnullData[unique_idx[0]];
            Z_xPoint = znullData[unique_idx[0]];
        } else if (xpoint_index == -1) {
            // Multiple distinct X-points and no index specified — list them and abort.
            std::cerr << "Error: HDF5 file contains " << unique_idx.size()
                      << " distinct X-points. Set xpoint_index in the input file to select one:\n";
            for (size_t k = 0; k < unique_idx.size(); ++k)
                std::cerr << "  xpoint_index = " << k
                          << "   R = " << xnullData[unique_idx[k]]
                          << "   Z = " << znullData[unique_idx[k]] << "\n";
            H5Sclose(xnullSpace);
            H5Sclose(znullSpace);
            H5Dclose(xnullDataset);
            H5Dclose(znullDataset);
            H5Fclose(file);
            return false;
        } else if (static_cast<size_t>(xpoint_index) >= unique_idx.size()) {
            std::cerr << "Error: xpoint_index=" << xpoint_index
                      << " is out of range (HDF5 file has " << unique_idx.size()
                      << " distinct X-points)." << std::endl;
            H5Sclose(xnullSpace);
            H5Sclose(znullSpace);
            H5Dclose(xnullDataset);
            H5Dclose(znullDataset);
            H5Fclose(file);
            return false;
        } else {
            R_xPoint = xnullData[unique_idx[xpoint_index]];
            Z_xPoint = znullData[unique_idx[xpoint_index]];
        }

        H5Sclose(xnullSpace);
        H5Sclose(znullSpace);
        H5Dclose(xnullDataset);
        H5Dclose(znullDataset);
        H5Fclose(file);

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return false;
    }

    return true;
}