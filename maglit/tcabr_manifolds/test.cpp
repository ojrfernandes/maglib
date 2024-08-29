// #include "findx.h"

#include "../maglit.h"
#include <armadillo>
#include <iomanip>
#include <utility>

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
            std::cout << "Newton's method reached maximum number of iterations" << std::endl;
            return false;
        }
    }
    return false;
}

// test the x_point function
int main() {
    char   source_path[] = "/home/jfernandes/Software/maglib_local/maglit/our_m3dc1_data/i_coils/n03/C1.h5";
    maglit tracer(source_path, FIO_M3DC1_SOURCE);
    tracer.inverse_map(false);
    double R_xpoint = 0.48;
    double Z_xpoint = -0.2;
    double Phi = 0;

    bool found_xp = x_point(tracer, R_xpoint, Z_xpoint, Phi, 1e-12, 100);
    if (!found_xp) {
        std::cerr << "Failed to find the x-point" << std::endl;
        return 1;
    }

    std::cout << std::fixed << std::setprecision(15);
    std::cout << "x-point found at: " << std::endl;
    std::cout << "R: " << R_xpoint << " Z: " << Z_xpoint << std::endl;
    // std::cout << std::defaultfloat << std::setprecision(6);

    // Evaluate the jacobian at the x-point
    double jacobian[2][2];
    eval_jacobian(tracer, R_xpoint, Z_xpoint, Phi, 1e-5, jacobian);
    std::cout << "Jacobian at the x-point: " << std::endl;
    std::cout << jacobian[0][0] << " " << jacobian[0][1] << std::endl;
    std::cout << jacobian[1][0] << " " << jacobian[1][1] << std::endl;

    // Convert the jacobian to an armadillo matrix
    arma::mat jacobian_arma = {{jacobian[0][0], jacobian[0][1]}, {jacobian[1][0], jacobian[1][1]}};

    // Compute the eigenvalues and eigenvectors of the jacobian
    arma::cx_vec eigval;
    arma::cx_mat eigvec;

    arma::eig_gen(eigval, eigvec, jacobian_arma);

    std::cout << "Eigenvalues: " << std::endl;
    std::cout << eigval << std::endl;

    std::cout << "Eigenvectors: " << std::endl;
    std::cout << eigvec << std::endl;

    // Find the eigenvector corresponding to the eigenvalue with |lambda| > 1
    arma::cx_vec selected_eigvec;
    for (size_t i = 0; i < eigval.n_elem; ++i) {
        if (std::abs(eigval[i]) > 1) {
            selected_eigvec = eigvec.col(i);
            break;
        }
    }

    // Normalize the selected eigenvector to get the unit vector
    selected_eigvec = selected_eigvec / arma::norm(selected_eigvec);

    // Extract real parts of the eigenvector (assuming they are real)
    double v_R = selected_eigvec[0].real();
    double v_Z = selected_eigvec[1].real();
    std::cout << "Selected eigenvector (|lambda| > 1): " << std::endl;
    std::cout << v_R << std::endl;
    std::cout << v_Z << std::endl;
    std::cout << sqrt(v_R * v_R + v_Z * v_Z) << std::endl;

    // Define the small distance epsilon
    double epsilon = 1e-8;

    // Calculate the new point (R_new, Z_new)
    double R_new = R_xpoint + epsilon * v_R;
    double Z_new = Z_xpoint + epsilon * v_Z;

    std::cout << "New point in the vicinity x-point: " << std::endl;
    std::cout << "R_new: " << R_new << " Z_new: " << Z_new << std::endl;

    // Apply the map M to the new point: M(R_new, Z_new) = (MR, MZ)
    std::pair<double, double> M_rz = apply_map(tracer, R_new, Z_new, Phi);
    double                    MR = M_rz.first;
    double                    MZ = M_rz.second;

    std::cout << "Mapped point: " << std::endl;
    std::cout << "MR: " << MR << " MZ: " << MZ << std::endl;

    // Interpolate points between (R_new, Z_new) and (MR, MZ)
    int    num_points = 3000;
    double step_R = (MR - R_new) / (num_points - 1);
    double step_Z = (MZ - Z_new) / (num_points - 1);

    std::vector<double> R_vec(num_points);
    std::vector<double> Z_vec(num_points);

    for (int i = 0; i < num_points; ++i) {
        R_vec[i] = R_new + i * step_R;
        Z_vec[i] = Z_new + i * step_Z;
    }

    // Compute manifold
    int    crossings = 12;
    double phi_max = 2 * M_PI;
    double Phi_aux = 0;
    FILE  *f1 = fopen("crossings/crossings_icoil_n3_1_.dat", "a");
    int    status = SODE_CONTINUE_GOOD_STEP;
    for (int j = 0; j < num_points; j++) {
        for (int i = 0; i < crossings; i++) {
            tracer.reset();
            tracer.alloc_hint();
            Phi_aux = Phi;
            do {
                status = tracer.step(R_vec[j], Z_vec[j], Phi_aux, phi_max, 0);
            } while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP);
            Phi_aux = Phi;
            if (status == SODE_SUCCESS_TIME) {
                fprintf(f1, "%f  %f  %f\n", R_vec[j], Z_vec[j], Phi_aux);
            } else {
                printf("SODE_ERROR\n");
            }
            tracer.clear_hint();
        }
        // print j
        printf("j: %d\n", j);
    }
    fclose(f1);

    return 0;
}