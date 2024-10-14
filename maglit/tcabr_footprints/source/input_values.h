#ifndef INPUT_VALUES_H
#define INPUT_VALUES_H

#include <cstring>
#include <fstream>
#include <iostream>
#include <unistd.h>

class input_values {
  public:
    // Constructor to initialize reading_path
    input_values(const std::string &readingPath);
    // Destructor to clean up dynamically allocated memory
    ~input_values();
    // read paths to the source, shape and output from a text file along with the initial grid parameters
    bool readInputFile();

    // Variables to store the paths and parameters
    char *source_path;
    char *shape_path;
    char *output_path;
    int num_theads;
    int plate;
    int timeslice;
    double gridMin;
    double gridMax;
    int nGrid;
    int nPhi;

  private:
    std::string reading_path; // File path
};

#endif // INPUT_VALUES_H