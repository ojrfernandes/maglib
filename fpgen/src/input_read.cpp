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
        } else if (key == "shape_path") {
            this->shape_path = value;
        } else if (key == "output_path") {
            this->output_path = value;
        } else if (key == "num_threads") {
            this->num_threads = std::stoi(value);
        } else if (key == "plate") {
            this->plate == std::stoi(value);
        } else if (key == "timeslice") {
            this->timeslice = std::stoi(value);
        } else if (key == "nGrid") {
            this->nGrid = std::stoi(value);
        } else if (key == "nPhi") {
            this->nPhi = std::stoi(value);
        } else if (key == "gridMin") {
            this->gridMin = std::stod(value);
        } else if (key == "gridMax") {
            this->gridMax = std::stod(value);
        } else {
            std::cerr << "Error: Unknown key: " << key << std::endl;
        }
    }

    return true;
}