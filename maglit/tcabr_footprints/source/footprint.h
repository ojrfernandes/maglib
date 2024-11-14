#ifndef FOOTPRINT_H
#define FOOTPRINT_H
#define FOOTPRINT_V 241014 // version (yy.mm.dd)

#include "../../maglit.h"
#include <iomanip>
#include <omp.h>

class footprint {
  public:
    footprint(const double &plate, const double &gridMin, const double &gridMax, const int &nGrid, const int &nPhi);
    void runGrid(maglit &tracer);
    std::vector<std::vector<double>> outputData;

  private:
    typedef struct
    {
        double length;
        double psimin;
    } map_scalars;

    void evolve_line(maglit &tracer, double R0, double Z0, double phi0, double &R1, double &Z1, double &phi1, map_scalars &scalars);
    double connection_length(double R0, double Z0, double phi0, double R1, double Z1, double phi1);
    void progressBar(float progress);

    double gridMin;
    double gridMax;
    int nGrid;
    int nPhi;
    int plate;
};

#endif // FOOTPRINT_H