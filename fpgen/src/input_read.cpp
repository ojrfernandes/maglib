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
        } else if (key == "shape_path") {
            this->shape_path = value;
        } else if (key == "output_path") {
            this->output_path = value;
            // check if there is a saved file with the same name
            std::ifstream f0(this->output_path);
            if (f0.is_open()) {
                std::cerr << "Warning: File " << this->output_path
                          << " already exists. Please, change your output file name to avoid overwriting your data."
                          << std::endl;
                return false;
            }
        } else if (key == "num_threads") {
            this->num_threads = std::stoi(value);
        } else if (key == "plate") {
            this->plate = std::stoi(value);
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
            std::cerr << "Error reading input: Unknown key: " << key << std::endl;
        }
    }

    return true;
}
