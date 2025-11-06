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
        } else if (key == "grid_R1") {
            this->grid_R1 = std::stod(value);
        } else if (key == "grid_Z1") {
            this->grid_Z1 = std::stod(value);
        } else if (key == "grid_R2") {
            this->grid_R2 = std::stod(value);
        } else if (key == "grid_Z2") {
            this->grid_Z2 = std::stod(value);
        } else if (key == "nRZ") {
            this->nRZ = std::stoi(value);
        } else if (key == "nPhi") {
            this->nPhi = std::stoi(value);
        } else if (key == "num_threads") {
            this->num_threads = std::stoi(value);
            // check if num_threads is a positive integer
            if (this->num_threads <= 0) {
                std::cerr << "Error: num_threads must be a positive integer." << std::endl;
                return false;
            }
        } else if (key == "max_turns") {
            this->max_turns = std::stoi(value);
            // check if max_turns is a positive integer
            if (this->max_turns <= 0) {
                std::cerr << "Error: max_turns must be a positive integer." << std::endl;
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
            // check if h_max is greater than or equal to h_min
            if (this->h_max < this->h_min) {
                std::cerr << "Error: h_max must be greater than or equal to h_min." << std::endl;
                return false;
            }
        } else {
            std::cerr << "Error reading input: Unknown key: " << key << std::endl;
        }
    }

    return true;
}
