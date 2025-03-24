#ifndef INPUT_READ_H
#define INPUT_READ_H
#define INPUT_READ_V 250324 // version (yy.mm.dd)

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

class input_read {
  public:
    // Constructor to initialize reading_path
    input_read(const std::string &readingPath);
    // read paths to the source, shape and output from a text file along with the initial grid parameters
    bool readInputFile();

    // Variables to store the paths and parameters
    std::string source_path; // path to the source file
    std::string shape_path;  // path to the first wall shape file
    std::string output_path; // path to the output file
    int num_threads;         // number of threads to be used
    int plate;               // mapped divertor plate (floor=0; wall=1)
    int timeslice;           // timeslice to be read from the source file
    int nGrid;               // number of grid points in the y coordinate (R or Z) of the grid
    int nPhi;                // number of grid points in the phi coordinate of the grid
    double gridMin;          // minimum value for the y coordinate (R or Z) of the grid
    double gridMax;          // maximum value for the y coordinate (R or Z) of the grid

  private:
    std::string reading_path; // File path
};

#endif // INPUT_READ_H