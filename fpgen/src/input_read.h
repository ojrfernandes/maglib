#ifndef INPUT_READ_H
#define INPUT_READ_H
#define INPUT_READE_V 241014 // version (yy.mm.dd)

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
    char *source_path;       // path to the source file
    std::string shape_path;  // path to the first wall shape file
    std::string output_path; // path to the output file
    int num_theads;          // number of threads to be used
    int plate;               // mapped divertor plate (floor=0; wall=1)
    int timeslice;           // timeslice to be read from the source file
    double gridMin;          // minimum value for the y coordinate (R or Z) of the grid
    double gridMax;          // maximum value for the y coordinate (R or Z) of the grid
    int nGrid;               // number of grid points in the y coordinate (R or Z) of the grid
    int nPhi;                // number of grid points in the phi coordinate of the grid

  private:
    std::string reading_path; // File path
};

#endif // INPUT_READ_H