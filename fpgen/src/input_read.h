#ifndef INPUT_READ_H
#define INPUT_READ_H
#define INPUT_READ_V 251106 // version (yy.mm.dd)

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

    // I/O files
    std::string source_path; // path to the source file
    std::string shape_path;  // path to the first wall shape file
    std::string output_path; // path to the output file

    // Mapping parameters
    int    timeslice; // M3DC1 timeslice
    int    manifold;  // unstable manifold=0; stable manifold=1
    double grid_R1;   // first point (R,Z) delimiting the target plate mapped surface
    double grid_Z1;   // first point (R,Z) delimiting the target plate mapped surface
    double grid_R2;   // second point (R,Z) delimiting the target plate mapped surface
    double grid_Z2;   // second point (R,Z) delimiting the target plate mapped surface
    int    nRZ;       // grid dimension along the (R,Z) plane
    int    nPhi;      // grid dimension along toroidal direction (Phi)

    // Additional parameters
    int    num_threads; // number of threads for OpenMP parallelization
    int    max_turns;   // maximum toroidal turns for field line integration
    double h_init;      // initial step size for field line integration
    double h_min;       // minimum step size for field line integration
    double h_max;       // maximum step size for field line integration

  private:
    std::string reading_path; // File path
};

#endif // INPUT_READ_H