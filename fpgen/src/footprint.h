#ifndef FOOTPRINT_H
#define FOOTPRINT_H
#define FOOTPRINT_V 251106 // version (yy.mm.dd)

#include "input_read.h"
#include <iomanip>
#include <maglit.h>
#include <omp.h>
#include <thread>

class footprint {
  public:
    // class constructor
    footprint(const int manifold, const double grid_R1, const double grid_Z1, const double grid_R2, const double grid_Z2, const int nRZ, const int nPhi, const int max_turns);
    void runGrid(maglit &tracer);

    std::vector<std::vector<double>> outputData; // two dimensional vector to store the output data

  private:
    typedef struct
    {
        double length;
        double psimin;
    } map_scalars; // structure to store the connection length and the minimum psi value

    // integrate the field line from a given initial condition
    void evolve_line(maglit &tracer, double R0, double Z0, double phi0, double &R1, double &Z1, double &phi1, map_scalars &scalars);
    // calculate the distance between two points
    double connection_length(double R0, double Z0, double phi0, double R1, double Z1, double phi1);
    // display a progress bar
    void progressBar(float progress);

    int    manifold;         // manifold type: unstable=0; stable=1
    double grid_R1;          // first point R delimiting the target plate mapped surface
    double grid_Z1;          // first point Z delimiting the target plate mapped surface
    double grid_R2;          // second point R delimiting the target plate mapped surface
    double grid_Z2;          // second point Z delimiting the target plate mapped surface
    int    nRZ;              // grid dimension along the (R,Z) plane
    int    nPhi;             // grid dimension along the phi direction
    int    max_turns = 1000; // maximum toroidal turns for field line integration
};

#endif // FOOTPRINT_H