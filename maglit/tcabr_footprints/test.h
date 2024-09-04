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

// follow field lines and return scalar values when crossing the walls
void map_wall(maglit &tracer, double R0, double Z0, double phi0, double &R1, double &Z1, double &phi1, map_scalars &scalars);
// print sode integrator status
void status_printer(int status);
// compute the distance integrated over the magnetic line
double dist(double R0, double Z0, double phi0, double R1, double Z1, double phi1);
// read paths to the source, shape and output from a text file along with the initial variables grid parameters
void readPaths(const std::string &readingPath, char *&source_path, char *&shape_path, char *&output_path, double &Rmin, double &Rmax, int &nR, int &nPhi);

#endif