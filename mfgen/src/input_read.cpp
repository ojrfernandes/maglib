#include "input_read.h"

// Constructor
input_read::input_read(const std::string &readingPath) : reading_path(readingPath) {}

bool input_read::readInputFile() {
    std::ifstream file(this->reading_path);
    std::string line;

    if (!file) {
        std::cerr << "Error: Unable to open file: " << reading_path << std::endl;
        return false;
    }

    while (std::getline(file, line)) {
        // Ignore empty lines and comments
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::istringstream iss(line);
        std::string key, equals, value;

        if (!(iss >> key >> equals)) {
            std::cerr << "Error: Invalid line (missing key or '=' symbol): " << line << std::endl;
            continue;
        }

        if (equals != "=") {
            std::cerr << "Error: Invalid format (expected '=' after key): " << line << std::endl;
            continue;
        }

        // Read the rest of the line as value, ignoring anything after '#'
        std::getline(iss, value, '#');

        // Trim whitespace from value
        size_t start = value.find_first_not_of(" \t");
        size_t end = value.find_last_not_of(" \t");
        if (start != std::string::npos) {
            value = value.substr(start, end - start + 1);
        } else {
            value.clear();
        }

        // Assign values to corresponding variables
        if (key == "source_path") {
            this->source_path = value;
        } else if (key == "output_path") {
            this->output_path = value;
        } else if (key == "stability") {
            this->stability = std::stoi(value);
        } else if (key == "timeslice") {
            this->timeslice = std::stoi(value);
        } else if (key == "Phi") {
            this->Phi = std::stod(value);
        } else if (key == "epsilon") {
            this->epsilon = std::stod(value);
        } else if (key == "l_lim") {
            this->l_lim = std::stod(value);
        } else if (key == "theta_lim") {
            this->theta_lim = std::stod(value);
        } else if (key == "nSeg") {
            this->nSeg = std::stoi(value);
        } else {
            std::cerr << "Error: Invalid key: " << key << std::endl;
        }
    }

    return true;
}

// read xnull and znull from the HDF5 file
bool input_read::readHDF5File() {
    // paths to the xnull and znull datasets in the HDF5 file
    const std::string xnullPath = "scalars/xnull";
    const std::string znullPath = "scalars/znull";

    try {
        // Open the HDF5 file
        hid_t file = H5Fopen(source_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

        // Open the datasets
        hid_t xnullDataset = H5Dopen(file, "scalars/xnull", H5P_DEFAULT);
        hid_t znullDataset = H5Dopen(file, "scalars/znull", H5P_DEFAULT);

        // Get the dataspace of the datasets
        hid_t xnullSpace = H5Dget_space(xnullDataset);
        hid_t znullSpace = H5Dget_space(znullDataset);

        // Get the number of elements in the datasets
        hsize_t xnullSize;
        hsize_t znullSize;
        H5Sget_simple_extent_dims(xnullSpace, &xnullSize, nullptr);
        H5Sget_simple_extent_dims(znullSpace, &znullSize, nullptr);

        // Read the data into vectors
        std::vector<double> xnullData(xnullSize);
        std::vector<double> znullData(znullSize);

        H5Dread(xnullDataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xnullData.data());
        H5Dread(znullDataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, znullData.data());

        // save the xnull and znull values to R_xpoint and Z_xpoint
        R_xPoint = xnullData[0];
        Z_xPoint = znullData[0];

        // Close datasets and file
        H5Dclose(xnullDataset);
        H5Dclose(znullDataset);
        H5Fclose(file);

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return false;
    }

    return true;
}