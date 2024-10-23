#ifndef INPUT_READ_H
#define INPUT_READ_H

#include <cstring>
#include <fstream>
#include <iostream>
#include <unistd.h>

class input_read {
  public:
    // Constructor to initialize reading_path
    input_read(const std::string &readingPath);
    // Destructor to clean up dynamically allocated memory
    ~input_read();
    // read paths to the source, shape and output from a text file along with the initial grid parameters
    bool readInputFile();

    // Variables to store the paths and parameters
    char *source_path;
    char *output_path;
    int timeslice;
    int stability;
    double Phi;
    double R_xPoint;
    double Z_xPoint;
    double epsilon;
    double l_lim;
    double theta_lim;
    double nSeg;

  private:
    std::string reading_path; // File path
};

#endif // INPUT_READ_H