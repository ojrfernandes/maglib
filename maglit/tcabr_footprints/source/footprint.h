#ifndef FOOTPRINT_H
#define FOOTPRINT_H
#define FOOTPRINT_V 241212 // version (yy.mm.dd)

#include "../../maglit.h"
#include <iomanip>
#include <omp.h>
#include <thread>

class footprint {
  public:
    // class constructor
    footprint(const double &plate, const double &gridMin, const double &gridMax, const int &nGrid, const int &nPhi);
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

    int plate;      // mapped divertor plate (floor=0; wall=1)
    double gridMin; // minimum value for the y coordinate (R or Z) of the grid
    double gridMax; // maximum value for the y coordinate (R or Z) of the grid
    int nGrid;      // number of grid points in the y coordinate (R or Z) of the grid
    int nPhi;       // number of grid points in the phi coordinate of the grid
};

#endif // FOOTPRINT_H