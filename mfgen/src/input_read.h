#ifndef INPUT_READ_H
#define INPUT_READ_H
#define INPUT_READ_V 251110 // version (yy.mm.dd)

#include <fstream>
#include <hdf5.h>
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
    // read xnull and znull from the HDF5 file
    bool readHDF5File();

    // I/O files
    std::string source_path; // path to the source file
    std::string output_path; // path to the output file

    // Tracing parameters
    int    timeslice; // timeslice to be read from the source file
    int    manifold;  // stability of the manifold (stable or unstable)
    int    method;    // method to be used for the manifold generation
    double Phi;       // toroidal angle coordinate for of the poincaré section
    int    nSegments; // number of primary segments to be mapped

    // Multiple Poincare sections
    int    nSections; // number of poincare sections to be generated
    double phi_0;     // initial toroidal angle coordinate
    double phi_1;     // final toroidal angle coordinate

    // Additional parameters
    double epsilon;        // first primary segment distance to the x-point
    double l_lim;          // distance treshold for the refinement process
    double theta_lim;      // angle treshold for the refinement process
    double h_init;         // initial step-size for integration
    double h_min;          // minimum step-size for integration
    double h_max;          // maximum step-size for integration
    double h_deriv;        // step-size for numerical derivatives
    double n_tol;          // tolerance for Newton's method
    int    max_iter;       // maximum iterations for Newton's method
    double precision;      // precision for floating point comparison
    int    max_insertions; // maximum number of points per segment

    // X-point coordinates for HDF5 reading
    double R_xPoint = 0; // R coordinate of the x-point
    double Z_xPoint = 0; // Z coordinate of the x-point

  private:
    std::string reading_path; // File path
};

#endif // INPUT_READ_H