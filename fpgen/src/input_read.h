#ifndef INPUT_READ_H
#define INPUT_READ_H
// Last modified: 26.06.10

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

class input_read {
  public:
    // Constructor to initialize reading_path
    input_read(const std::string &readingPath);
    // read paths to the source, first wall and output from a text file along with the initial grid parameters
    bool readInputFile();

    // I/O files
    std::string source_path;     // path to the source file
    std::string first_wall_path; // path to the first wall boundary file
    std::string output_path; // path to the output file

    // Mapping parameters
    int    timeslice   = -1;  // M3DC1 timeslice
    int    manifold    = 0;   // 0 = unstable (against-B from target);  1 = stable (following-B from target)
    double grid_R1     = 0.0; // first point (R,Z) delimiting the target plate mapped surface
    double grid_Z1     = 0.0; // first point (R,Z) delimiting the target plate mapped surface
    double grid_R2     = 0.0; // second point (R,Z) delimiting the target plate mapped surface
    double grid_Z2     = 0.0; // second point (R,Z) delimiting the target plate mapped surface
    int    nRZ         = 0;   // grid dimension along the (R,Z) plane
    int    nPhi        = 0;   // grid dimension along toroidal direction (Phi)

    // Additional parameters
    int    num_threads = 1;    // number of threads for OpenMP parallelization
    int    max_turns   = 1000; // maximum toroidal turns for field line integration
    double h_init      = 1e-2; // initial step size for field line integration
    double h_min       = 1e-6; // minimum step size for field line integration
    double h_max       = 1e-1; // maximum step size for field line integration

  private:
    std::string reading_path; // File path
};

#endif // INPUT_READ_H
