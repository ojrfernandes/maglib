#ifndef TEST_H
#define TEST_H

#include <stdio.h>
#include <unistd.h>

#include "../../sode/sode.cpp"
#include "../auxfields.cpp"
#include "../maglit.cpp"
#include "tcabr_collider.cpp"

typedef struct
{
    double deltaPhi;
    double length;
    double psimin;
} map_scalars;

// map coordinates in vessel and get scalar values
void map_wall(maglit &tracer, auxfields &aux_field, double R0, double Z0, double phi0, double &R1, double &Z1, double &phi1, map_scalars &scalars);
// print sode integrator status
void status_printer(int status);
// compute the distance integrated over the magnetic line
double dist(double R0, double Z0, double phi0, double R1, double Z1, double phi1);

#endif