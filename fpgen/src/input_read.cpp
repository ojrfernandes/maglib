#include "input_read.h"
#include <cstdio>

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

    std::string current_section = "";

    while (std::getline(file, line)) {
        // Trim the line first
        trim(line);

        // Ignore empty lines and lines that start with '#'
        if (line.empty() || line[0] == '#') {
            continue;
        }

        // Section header: [SECTION NAME]
        if (line.front() == '[' && line.back() == ']') {
            current_section = line.substr(1, line.size() - 2);
            trim(current_section);
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
        if (key == "first_wall_path") {
            this->first_wall_path = value;
        } else if (key == "output_path") {
            this->output_path = value;
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
            // check if nRZ is a positive integer
            if (this->nRZ <= 0) {
                std::cerr << "Error: nRZ must be a positive integer." << std::endl;
                return false;
            }
        } else if (key == "nPhi") {
            this->nPhi = std::stoi(value);
            // check if nPhi is a positive integer
            if (this->nPhi <= 0) {
                std::cerr << "Error: nPhi must be a positive integer." << std::endl;
                return false;
            }
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
            if (this->h_max <= 0) {
                std::cerr << "Error: h_max must be a positive number." << std::endl;
                return false;
            }
        } else if (key == "nsources") {
            this->nsources = std::stoi(value);
            if (this->nsources < 1) {
                std::cerr << "Error: nsources must be a positive integer." << std::endl;
                return false;
            }
        } else {
            // Pattern-match indexed superposition keys: source_N, timeslice_N, phase_N, amplitude_N
            int idx = -1;
            char prefix[32];
            if (sscanf(key.c_str(), "%31[a-z]_%d", prefix, &idx) == 2 && idx >= 0) {
                std::string pfx(prefix);
                if ((int)this->components.size() <= idx)
                    this->components.resize(idx + 1);
                if (pfx == "source") {
                    this->components[idx].path = value;
                } else if (pfx == "timeslice") {
                    this->components[idx].timeslice = std::stoi(value);
                    if (this->components[idx].timeslice < -1) {
                        std::cerr << "Error: timeslice_" << idx << " must be -1 or a non-negative integer." << std::endl;
                        return false;
                    }
                } else if (pfx == "phase") {
                    this->components[idx].phase = std::stod(value);
                } else if (pfx == "amplitude") {
                    this->components[idx].amplitude = std::stod(value);
                } else {
                    std::cerr << "Error reading input: Unknown key: " << key << std::endl;
                    return false;
                }
            } else {
                std::cerr << "Error reading input: Unknown key: " << key << std::endl;
                return false;
            }
        }
    }

    // Check that required string parameters are not empty
    if (first_wall_path.empty() || output_path.empty()) {
        std::cerr << "Error on input: One or more required file paths are missing." << std::endl;
        return false;
    }
    if (this->nsources < 1) {
        std::cerr << "Error: nsources must be set to a positive integer in [M3DC1 SOURCE]." << std::endl;
        return false;
    }
    if ((int)this->components.size() != this->nsources) {
        std::cerr << "Error: nsources = " << this->nsources
                  << " but " << this->components.size()
                  << " component(s) were defined." << std::endl;
        return false;
    }
    for (int i = 0; i < this->nsources; ++i) {
        if (this->components[i].path.empty()) {
            std::cerr << "Error: source_" << i << " path is missing." << std::endl;
            return false;
        }
    }

    // Validate step-size relationships after all keys are parsed
    if (this->h_max < this->h_min) {
        std::cerr << "Error: h_max must be greater than or equal to h_min." << std::endl;
        return false;
    }

    return true;
}
