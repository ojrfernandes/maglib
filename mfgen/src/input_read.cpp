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
        } else if (key == "xpoint_null") {
            this->xpoint_null = std::stoi(value);
            if (this->xpoint_null < 0 || this->xpoint_null > 2) {
                std::cerr << "Error: xpoint_null must be 0 (auto), 1 (xnull/znull), or 2 (xnull2/znull2)." << std::endl;
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

// Helper: read first element of a 1-D float32 HDF5 dataset as double.
// Returns false if the dataset cannot be opened.
static bool readFirstElement(hid_t file, const char *path, double &out) {
    hid_t ds = H5Dopen(file, path, H5P_DEFAULT);
    if (ds < 0) return false;
    float val;
    H5Dread(ds, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &val);
    H5Dclose(ds);
    out = static_cast<double>(val);
    return true;
}

// Read X-point initial guess from the M3DC1 C1.h5 file.
// M3DC1 stores coordinates in separate datasets (float32, shape (N,) with repeated values):
//   scalars/xnull  + scalars/znull  — primary null
//   scalars/xnull2 + scalars/znull2 — secondary null (zero-sentinel in single-null cases)
// Selection is controlled by xpoint_null: 0=auto, 1=primary, 2=secondary.
bool input_read::readHDF5File() {
    hid_t file = H5Fopen(source_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        std::cerr << "Error: could not open HDF5 file: " << source_path << std::endl;
        return false;
    }

    double R1 = 0.0, Z1 = 0.0;
    if (!readFirstElement(file, "scalars/xnull", R1) ||
        !readFirstElement(file, "scalars/znull", Z1)) {
        std::cerr << "Error: could not read scalars/xnull or scalars/znull from HDF5 file." << std::endl;
        H5Fclose(file);
        return false;
    }

    double R2 = 0.0, Z2 = 0.0;
    bool has_secondary = readFirstElement(file, "scalars/xnull2", R2) &&
                         readFirstElement(file, "scalars/znull2", Z2) &&
                         !(R2 == 0.0 && Z2 == 0.0);

    H5Fclose(file);

    if (xpoint_null == 1) {
        R_xPoint = R1;
        Z_xPoint = Z1;
    } else if (xpoint_null == 2) {
        if (!has_secondary) {
            std::cerr << "Error: xpoint_null=2 requested but no secondary X-point found "
                      << "(single-null configuration or xnull2/znull2 absent)." << std::endl;
            return false;
        }
        R_xPoint = R2;
        Z_xPoint = Z2;
    } else {
        // auto: valid only for single-null
        if (has_secondary) {
            std::cerr << "Error: double-null configuration detected. "
                      << "Set xpoint_null in the input file to select one:\n"
                      << "  xpoint_null = 1   R = " << R1 << "   Z = " << Z1 << "\n"
                      << "  xpoint_null = 2   R = " << R2 << "   Z = " << Z2 << "\n";
            return false;
        }
        R_xPoint = R1;
        Z_xPoint = Z1;
    }

    return true;
}