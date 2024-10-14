#include "input_values.h"

// Constructor
input_values::input_values(const std::string &readingPath) : source_path(nullptr), shape_path(nullptr), output_path(nullptr), reading_path(readingPath) {}

// Destructor to free allocated memory
input_values::~input_values() {
    delete[] source_path;
    delete[] shape_path;
    delete[] output_path;
}

// read paths to the source, shape and output from a text file along with the initial grid parameters
bool input_values::readInputFile() {
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
                // Allocate memory for shape_path and copy the string
                shape_path = new char[line.length() + 1];
                std::strcpy(shape_path, line.c_str());
                break;
            case 2:
                // Allocate memory for output_path and copy the string
                output_path = new char[line.length() + 1];
                std::strcpy(output_path, line.c_str());
                // Check if the output file already exists
                if (access(output_path, F_OK) != -1) {
                    std::cerr << "There is a file with the same name saved to the chosen directory. "
                              << "Please change the output file name to avoid overwriting your data."
                              << std::endl;
                    return false; // Exit the function to avoid overwriting
                }
                break;
            case 3:
                num_theads = std::stoi(line);
                break;
            case 4:
                plate = std::stoi(line);
                break;
            case 5:
                timeslice = std::stoi(line);
                break;
            case 6:
                gridMin = std::stod(line);
                break;
            case 7:
                gridMax = std::stod(line);
                break;
            case 8:
                nGrid = std::stoi(line);
                break;
            case 9:
                nPhi = std::stoi(line);
                break;
            default:
                std::cerr << "Unexpected line index: " << line_index << std::endl;
                return false;
            }
            line_index++;
        }
    }

    if (line_index < 10) {
        std::cerr << "File does not contain enough data. Expected 10 values, found " << line_index << "." << std::endl;
        return false;
    }

    pathFile.close();
    return true;
}