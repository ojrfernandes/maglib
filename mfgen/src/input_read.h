#ifndef INPUT_READ_H
#define INPUT_READ_H
#define INPUT_READ_V 250528 // version (yy.mm.dd)

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

    // Variables to store the paths and parameters
    std::string source_path; // path to the source file
    std::string output_path; // path to the output file
    int timeslice;           // timeslice to be read from the source file
    int stability;           // stability of the manifold (stable or unstable)
    int method;              // method to be used for the manifold generation
    int slices;              // number of toroidal slices to be generated
    double Phi;              // toroidal angle coordinate for of the poincaré section
    double phi_0;            // initial toroidal angle coordinate
    double phi_1;            // final toroidal angle coordinate
    double R_xPoint = 0;     // R coordinate of the x-point
    double Z_xPoint = 0;     // Z coordinate of the x-point
    double epsilon;          // first primary segment distance to the x-point
    double l_lim;            // distance treshold for the refinement process
    double theta_lim;        // angle treshold for the refinement process
    int nSeg;                // number of primary segments to be mapped

  private:
    std::string reading_path; // File path
};

#endif // INPUT_READ_H