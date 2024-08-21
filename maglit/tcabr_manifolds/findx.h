#ifndef FINDX_H
#define FINDX_H

#include "../maglit.h"
#include <iomanip>
#include <utility>

// Apply the map to a given point
std::pair<double, double> apply_map(maglit &source, double R, double Z, double Phi);
// Function to evaluate the planar map's jacobian
void eval_jacobian(maglit &source, double R, double Z, double Phi, double h, double jacobian[2][2]);
// Iteratively find the closest 1 period fixed point from the initial guess
bool x_point(maglit &source, double &R, double &Z, double &Phi, double tol, int max_iter);

#endif