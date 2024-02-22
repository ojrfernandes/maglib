#ifndef TCABR_COLLIDER
#define TCABR_COLLIDER

#include <cmath>
#include <stdio.h>

typedef struct {
    int     np;  // number of nodes
    double *R;   // node coords
    double *Z;   //
    double  Rc;  // center coords
    double  Zc;  //
    double *th;  // sector angles
    int     idx; // current index
} tcabr_shape;

//
bool load_shape(const char path[], tcabr_shape *shape);
//
void free_shape(tcabr_shape *shape);
//
bool search_index(double ang, tcabr_shape *shape);
//
double cross(double RA, double ZA, double RB, double ZB, double RC, double ZC);
//
bool tcabr_inside(double R, double Z, double phi, void *aux);

#endif