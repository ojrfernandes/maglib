#include "input_read.h"

// Constructor
input_read::input_read(const std::string &readingPath) : reading_path(readingPath) {}

bool input_read::readInputFile() {
    std::ifstream file(this->reading_path);
    std::string   line;

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

        // Ignore empty lines and comments
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
            std::cerr << "Error reading input: Unknown key: " << key << std::endl;
        }
    }

    // Check that required string parameters are not empty
    if (this->hdf5File.empty() || this->equilibriumFile.empty() || this->perturbedFile.empty() || this->intersectionFile.empty() || this->lobeFile.empty()) {
        std::cerr << "Error on input: One or more required file paths are missing." << std::endl;
        return false;
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