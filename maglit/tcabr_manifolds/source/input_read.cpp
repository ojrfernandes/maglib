#include "input_read.h"

// Constructor
input_read::input_read(const std::string &readingPath) : source_path(nullptr), reading_path(readingPath) {}

// Destructor to free allocated memory
input_read::~input_read() {
    delete[] source_path;
}

// read paths to the source, shape and output from a text file along with the initial grid parameters
bool input_read::readInputFile() {
    std::ifstream pathFile(this->reading_path);

    if (!pathFile.is_open()) {
        std::cerr << "Unable to open file: " << this->reading_path << std::endl;
        return false;
    }

    std::string line;
    int line_index = 0;

    while (std::getline(pathFile, line)) {
        if (!line.empty() && line[0] != '#') { // Ignore empty lines and comments starting with #
            switch (line_index) {
            case 0:
                // Allocate memory for source_path and copy the string
                source_path = new char[line.length() + 1];
                std::strcpy(source_path, line.c_str());
                break;
            case 1:
                // Allocate memory for output_path and copy the string
                output_path = line;
                // Check if the output file already exists
                if (access(output_path.c_str(), F_OK) != -1) {
                    std::cerr << "There is a file with the same name saved to the chosen directory. "
                              << "Please change the output file name to avoid overwriting your data."
                              << std::endl;
                    return false; // Exit the function to avoid overwriting
                }
                break;
            case 2:
                stability = std::stoi(line);
                break;
            case 3:
                timeslice = std::stoi(line);
                break;
            case 4:
                Phi = std::stod(line);
                break;
            case 5:
                epsilon = std::stod(line);
                break;
            case 6:
                nSeg = std::stoi(line);
                break;
            case 7:
                l_lim = std::stod(line);
                break;
            case 8:
                theta_lim = std::stod(line);
                break;
            default:
                std::cerr << "Unexpected line index: " << line_index << std::endl;
                return false;
            }
            line_index++;
        }
    }

    if (line_index < 9) {
        std::cerr << "File does not contain enough data. Expected 10 values, found " << line_index << "." << std::endl;
        return false;
    }

    pathFile.close();

    // call readHDF5File to read xnull and znull from the HDF5 file
    if (!readHDF5File()) {
        std::cerr << "Error reading xnull and znull from the HDF5 file." << std::endl;
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
        // Open the HDF5 file
        hid_t file = H5Fopen(source_path, H5F_ACC_RDONLY, H5P_DEFAULT);

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