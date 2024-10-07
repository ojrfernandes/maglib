#ifndef TCABR_COLLIDER
#define TCABR_COLLIDER

#include <cmath>
#include <stdio.h>

typedef struct {
    int np;     // number of nodes
    double *R;  // node coords
    double *Z;  //
    double Rc;  // center coords
    double Zc;  //
    double *th; // sector angles
    int idx;    // current index
} tcabr_shape;

// load vessel shape from file
bool load_shape(const char path[], tcabr_shape *shape);
// free allocated memory
void free_shape(tcabr_shape *shape);
// search the index of a given angle
bool search_index(double ang, tcabr_shape *shape);
// calculate AB x AC
double cross(double RA, double ZA, double RB, double ZB, double RC, double ZC);
// check if coordinate lies inside the vessel
bool tcabr_inside(double R, double Z, double phi, void *aux);

#endif