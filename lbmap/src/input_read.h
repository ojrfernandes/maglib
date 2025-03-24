#ifndef INPUT_READ_H
#define INPUT_READ_H
#define INPUT_READ_V 250324 // version (yy.mm.dd)

#include <fstream>
#include <hdf5.h>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

class input_read {
  public:
    // Constructor
    input_read(const std::string &readingPath);
    // read paths to the input and output paths
    bool readInputFile();
    // read xmag and zmag from the HDF5 file
    bool readHDF5File();

    // Variables to store paths and coordinates
    std::string hdf5File;
    std::string equilibriumFile;
    std::string perturbedFile;
    std::string intersectionFile;
    std::string lobeFile;
    double xmag; // R coordinate of the x-point
    double zmag; // Z coordinate of the x-point

    // Set verbose to debug
    bool verbose = false;

  private:
    std::string reading_path; // File path
};

#endif // INPUT_READ_H