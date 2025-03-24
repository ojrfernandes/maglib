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
        if (key == "hdf5File") {
            this->hdf5File = value;
        } else if (key == "equilibriumFile") {
            this->equilibriumFile = value;
        } else if (key == "perturbedFile") {
            this->perturbedFile = value;
        } else if (key == "intersectionFile") {
            this->intersectionFile = value;
        } else if (key == "lobeFile") {
            this->lobeFile = value;
        } else {
            std::cerr << "Error: Unknown key: " << key << std::endl;
        }
    }

    if (verbose == true) {
        std::cout << "hdf5File: " << this->hdf5File << "\n";
        std::cout << "equilibriumFile: " << this->equilibriumFile << "\n";
        std::cout << "perturbedFile: " << this->perturbedFile << "\n";
        std::cout << "intersectionFile: " << this->intersectionFile << "\n";
        std::cout << "lobeFile: " << this->lobeFile << std::endl;
    }

    return true;
}

// read coordinates of xmag and zmag from hdf5 file
bool input_read::readHDF5File() {
    // paths to the xmag and zmag datasets in the HDF5 file
    const std::string xmagPath = "scalars/xmag";
    const std::string zmagPath = "scalars/zmag";

    try {
        // Open the HDF5 file
        hid_t file = H5Fopen(this->hdf5File.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

        // Open the datasets
        hid_t xmagDataset = H5Dopen(file, xmagPath.c_str(), H5P_DEFAULT);
        hid_t zmagDataset = H5Dopen(file, zmagPath.c_str(), H5P_DEFAULT);

        // Get the dataspace of the datasets
        hid_t xmagSpace = H5Dget_space(xmagDataset);
        hid_t zmagSpace = H5Dget_space(zmagDataset);

        // Get the number of elements in the datasets
        hsize_t xmagSize;
        hsize_t zmagSize;
        H5Sget_simple_extent_dims(xmagSpace, &xmagSize, nullptr);
        H5Sget_simple_extent_dims(zmagSpace, &zmagSize, nullptr);

        // Read the data into vectors
        std::vector<double> xmagData(xmagSize);
        std::vector<double> zmagData(zmagSize);

        H5Dread(xmagDataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xmagData.data());
        H5Dread(zmagDataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, zmagData.data());

        // save the xmag and zmag values to R_xmag and Z_xmag
        this->xmag = xmagData[0];
        this->zmag = zmagData[0];

        // Close datasets and file
        H5Dclose(xmagDataset);
        H5Dclose(zmagDataset);
        H5Fclose(file);
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return false;
    }

    return true;
}