#ifndef FINDX_H
#define FINDX_H

#include "../maglit.h"

double eval_psi(maglit &source, double R, double Z, double Phi);
void   eval_field_derivative(maglit &source, double R, double Z, double Phi, double *gradient, double hessian[2][2], double h = 1e-6);
bool   newton_raphson(maglit &source, double &R, double &Z, double &Phi, double tol = 1e-6, int max_iter = 1000);

#endif