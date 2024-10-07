#ifndef TEST_H
#define TEST_H

#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <unistd.h>

#include "../maglit.h"
#include "tcabr_collider.h"

typedef struct
{
    double deltaPhi;
    double length;
    double psimin;
} map_scalars;

//
void floor_grid(double nPhi, double nR, double Rmin, double Rmax, double Zfloor, maglit &tracer, map_scalars &scalars, std::vector<std::string> &dataWrite);
//
void wall_grid(double nPhi, double nZ, double Zmin, double Zmax, double Rfloor, maglit &tracer, map_scalars &scalars, std::vector<std::string> &dataWrite);
// follow field lines and return scalar values when crossing the walls
void evolve_lines(maglit &tracer, double R0, double Z0, double phi0, double &R1, double &Z1, double &phi1, map_scalars &scalars);
// print sode integrator status
void status_printer(int status);
// compute the distance integrated over the magnetic line
double dist(double R0, double Z0, double phi0, double R1, double Z1, double phi1);
// read paths to the source, shape and output from a text file along with the initial variables grid parameters
void readParams(const std::string &readingPath, char *&source_path, char *&shape_path, char *&output_path, int &plate, int &timeslice, double &gridMin, double &gridMax, int &nGrid, int &nPhi);

#endif