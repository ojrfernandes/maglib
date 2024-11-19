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
                // Copy the string to shape_path
                shape_path = line;
                break;
            case 2:
                // Copy the string to output_path
                output_path = line;
                // Check if the output file already exists
                if (access(output_path.c_str(), F_OK) != -1) {
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