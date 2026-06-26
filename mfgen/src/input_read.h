#ifndef INPUT_READ_H
#define INPUT_READ_H
// Last modified: 26.06.26

#include <fstream>
#include <hdf5.h>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
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
    std::string output_path; // path to the output file

    // [M3DC1 SOURCE] section — always populated, even for nsources=1.
    struct SourceComponent {
        std::string path;
        int    timeslice  = 1;
        double phase      = 0.0; // degrees
        double amplitude  = 1.0;
    };
    int                       nsources   = 0;  // must be set explicitly in [M3DC1 SOURCE]
    std::vector<SourceComponent> components;   // indexed 0..nsources-1

    // Tracing parameters
    int    manifold    = 0;   // 0 = unstable (forward map);  1 = stable (inverse map)
    int    method      = 0;   // method: exact-map=0, interpolant=1
    double Phi         = 0.0; // toroidal angle of the Poincaré section (rad)
    int    nSegments   = 0;   // number of primary segments to be mapped

    // Multiple Poincare sections
    int    nSections   = 0;   // number of Poincaré sections to generate
    double phi_0       = 0.0; // initial toroidal angle coordinate
    double phi_1       = 0.0; // final toroidal angle coordinate

    // Additional parameters
    double epsilon       = 0.0;  // first primary segment distance to the x-point
    double l_lim         = 0.0;  // distance threshold for the refinement process
    double theta_lim     = 0.0;  // angle threshold for the refinement process
    double h_init        = 1e-2; // initial step-size for integration
    double h_min         = 1e-6; // minimum step-size for integration
    double h_max         = 1e-1; // maximum step-size for integration
    double h_deriv       = 0.0;  // step-size for numerical derivatives
    double n_tol         = 0.0;  // tolerance for Newton's method
    int    max_iter      = 0;    // maximum iterations for Newton's method
    double precision     = 0.0;  // precision for floating point comparison
    int    max_insertions = 0;   // maximum number of points per segment
    int    verbose       = 0;    // verbose output flag (0=off, 1=on)

    // Pivot direction branch signs (±1 each; default 1 = no flip)
    int    branch_R   = 1;
    int    branch_Z   = 1;

    // X-point coordinates (optional; if neither is set, falls back to HDF5 auto-read)
    double R_xPoint   = 0.0;
    double Z_xPoint   = 0.0;
    bool   xpoint_set = false; // true when R_xPoint or Z_xPoint is specified in input
    int    xpoint_null = 0;    // HDF5 null selector: 0=auto, 1=xnull/znull, 2=xnull2/znull2

  private:
    std::string reading_path; // File path
};

#endif // INPUT_READ_H
