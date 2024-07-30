#include "findx.h"

// Function to evaluate the field at a given point
double eval_psi(maglit &source, double R, double Z, double Phi) {
    double psi;
    source.psi_eval(R, Phi, Z, &psi);
    return psi;
}

// Function to evaluate the field's gradient and Hessian
void eval_field_derivative(maglit &source, double R, double Z, double Phi, double *gradient, double hessian[2][2], double h) {

    // psi field at the point (R, Z)
    double f_xy = eval_psi(source, R, Z, Phi);

    // psi values at the four points around the point (R, Z)
    double f_x1 = eval_psi(source, R + h, Z, Phi);
    double f_x2 = eval_psi(source, R - h, Z, Phi);
    double f_y1 = eval_psi(source, R, Z + h, Phi);
    double f_y2 = eval_psi(source, R, Z - h, Phi);

    // gradient (fisr order aproximation using finite differences)
    gradient[0] = (f_x1 - f_x2) / (2 * h);
    gradient[1] = (f_y1 - f_y2) / (2 * h);

    // psi values at the four points around the point (R, Z)
    double f_xy1 = eval_psi(source, R + h, Z + h, Phi);
    double f_xy2 = eval_psi(source, R - h, Z + h, Phi);
    double f_xy3 = eval_psi(source, R + h, Z - h, Phi);
    double f_xy4 = eval_psi(source, R - h, Z - h, Phi);

    // Hessian (first order aproximation using finite differences)
    double f_xx = (f_x1 - 2 * f_xy + f_x2) / (h * h);
    double f_yy = (f_y1 - 2 * f_xy + f_y2) / (h * h);
    double d2fdxdy = (f_xy1 - f_xy2 - f_xy3 + f_xy4) / (4 * h * h);

    hessian[0][0] = f_xx;
    hessian[0][1] = d2fdxdy;
    hessian[1][0] = d2fdxdy;
    hessian[1][1] = f_yy;
}

// Newton-Raphson method to find the closest saddle point from the initial guess
bool newton_raphson(maglit &source, double &R, double &Z, double &Phi, double tol, int max_iter) {
    source.alloc_hint();
    for (int iter = 0; iter < max_iter; ++iter) {
        // allocate memory for the gradient and hessian
        double gradient[2];
        double hessian[2][2];
        // eval the field's gradient and derivative
        eval_field_derivative(source, R, Z, Phi, gradient, hessian);

        // check if gradient is close enough to zero
        if (std::sqrt(gradient[0] * gradient[0] + gradient[1] * gradient[1]) < tol) {
            return true;
        }

        // check if hessian is small enough to be singular, hence non invertible
        double det = hessian[0][0] * hessian[1][1] - hessian[0][1] * hessian[1][0];
        if (std::abs(det) < 1e-12) {
            return false;
        }

        // allocate memory for the inverse hessian
        double inv_hessian[2][2];
        inv_hessian[0][0] = hessian[1][1] / det;
        inv_hessian[0][1] = -hessian[0][1] / det;
        inv_hessian[1][0] = -hessian[1][0] / det;
        inv_hessian[1][1] = hessian[0][0] / det;

        // update the position
        double delta_R = -(inv_hessian[0][0] * gradient[0] + inv_hessian[0][1] * gradient[1]);
        double delta_Z = -(inv_hessian[1][0] * gradient[0] + inv_hessian[1][1] * gradient[1]);
        R += delta_R;
        Z += delta_Z;

        // check if incremented length is small enough
        if (std::sqrt(delta_R * delta_R + delta_Z * delta_Z) < tol) {
            source.clear_hint();
            return true;
        }
    }
    source.clear_hint();
    return false;
}