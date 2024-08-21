#include "findx.h"

// Apply the map to a given point
std::pair<double, double> apply_map(maglit &source, double R, double Z, double Phi) {
    int status = SODE_CONTINUE_GOOD_STEP;
    source.reset();
    source.alloc_hint();
    double phi_max = 2 * M_PI;

    do {
        status = source.step(R, Z, Phi, phi_max, 0);
    } while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP);
    source.clear_hint();

    if (status == SODE_SUCCESS_TIME) {
        return std::make_pair(R, Z);
    } else {
        std::cerr << "Failed to apply map when evaluating the x-point" << std::endl;
        return std::make_pair(std::nan(""), std::nan(""));
    }
}

// Function to evaluate the planar map's jacobian
void eval_jacobian(maglit &source, double R, double Z, double Phi, double h, double jacobian[2][2]) {

    std::pair<double, double> point_1 = apply_map(source, R + h, Z, Phi);
    std::pair<double, double> point_2 = apply_map(source, R - h, Z, Phi);
    std::pair<double, double> point_3 = apply_map(source, R, Z + h, Phi);
    std::pair<double, double> point_4 = apply_map(source, R, Z - h, Phi);

    jacobian[0][0] = (point_1.first - point_2.first) / (2 * h);
    jacobian[0][1] = (point_3.first - point_4.first) / (2 * h);
    jacobian[1][0] = (point_1.second - point_2.second) / (2 * h);
    jacobian[1][1] = (point_3.second - point_4.second) / (2 * h);
}

// Iteratively find the closest 1 period fixed point from the initial guess
bool x_point(maglit &source, double &R, double &Z, double &Phi, double tol, int max_iter) {
    for (int iter = 0; iter < max_iter; ++iter) {
        std::cout << "iter: " << iter << std::endl;
        // initialize jacobian and identity
        double jacobian[2][2];
        eval_jacobian(source, R, Z, Phi, 1e-5, jacobian);
        double identity[2][2] = {{1, 0}, {0, 1}};

        // computer (identity - jacobian)
        double diff[2][2];
        diff[0][0] = identity[0][0] - jacobian[0][0];
        diff[0][1] = identity[0][1] - jacobian[0][1];
        diff[1][0] = identity[1][0] - jacobian[1][0];
        diff[1][1] = identity[1][1] - jacobian[1][1];

        // check if (identity - jacobian) is singular
        double det = diff[0][0] * diff[1][1] - diff[0][1] * diff[1][0];
        if (std::abs(det) < 1e-12) {
            std::cout << "(I - J) is a singular matrix" << std::endl;
            return false;
        }

        // compute (identity - jacobian)^(-1)
        double inv_diff[2][2];
        inv_diff[0][0] = diff[1][1] / det;
        inv_diff[0][1] = -diff[0][1] / det;
        inv_diff[1][0] = -diff[1][0] / det;
        inv_diff[1][1] = diff[0][0] / det;

        // compute (R1,RZ)
        std::pair<double, double> R1_Z1 = apply_map(source, R, Z, Phi);
        double                    R1 = R1_Z1.first;
        double                    Z1 = R1_Z1.second;

        // update the position
        double delta_R = R1 - (jacobian[0][0] * R + jacobian[0][1] * Z);
        double delta_Z = Z1 - (jacobian[1][0] * R + jacobian[1][1] * Z);
        R = inv_diff[0][0] * delta_R + inv_diff[0][1] * delta_Z;
        Z = inv_diff[1][0] * delta_R + inv_diff[1][1] * delta_Z;

        // check if the increment is small enough
        double dR = R - R1;
        double dZ = Z - Z1;
        if (std::sqrt(dR * dR + dZ * dZ) < tol) {
            return true;
        }

        // check if iteration limit is reached
        if (iter == max_iter - 1) {
            std::cout << "Newton method reached maximum number of iterations" << std::endl;
            return false;
        }
    }
    return false;
}

// test the x_point function
// int main() {
//     char   source_path[] = "/home/jfernandes/Software/maglib_local/maglit/our_m3dc1_data/i_coils/n03/C1.h5";
//     maglit tracer(source_path, FIO_M3DC1_SOURCE);
//     double R = 0.46;
//     double Z = -0.21;
//     double Phi = 0;

//     bool status = x_point(tracer, R, Z, Phi, 1e-12, 100);
//     if (!status) {
//         std::cerr << "Failed to find the x-point" << std::endl;
//         return 1;
//     }

//     std::cout << std::fixed << std::setprecision(15);
//     std::cout << "R: " << R << " Z: " << Z << std::endl;
//     return 0;
// }