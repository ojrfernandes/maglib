// #include "findx.h"

#include "../maglit.h"
#include <armadillo>
#include <chrono>
#include <iomanip>
#include <thread>
#include <utility>

// Apply map to a given point returning a pair (R, Z)
std::pair<double, double> apply_map(maglit &source, double R, double Z, double Phi, int nTurns) {
    // Reset the source and allocate hint
    int status = SODE_CONTINUE_GOOD_STEP;
    source.reset();
    source.alloc_hint();
    double phi_max = nTurns * 2 * M_PI;

    // Apply the map
    do {
        status = source.step(R, Z, Phi, phi_max, 0);
    } while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP);
    source.clear_hint();

    if (status == SODE_SUCCESS_TIME) {
        return std::make_pair(R, Z);
    } else {
        std::cerr << "Failed to apply map" << std::endl;
        return std::make_pair(std::nan(""), std::nan(""));
    }
}

// Function to evaluate the planar map's jacobian
void eval_jacobian(maglit &source, double R, double Z, double Phi, double h, double jacobian[2][2]) {

    // Apply maps to the points (R + h, Z), (R - h, Z), (R, Z + h), (R, Z - h)
    std::pair<double, double> point_1 = apply_map(source, R + h, Z, Phi, 1);
    std::pair<double, double> point_2 = apply_map(source, R - h, Z, Phi, 1);
    std::pair<double, double> point_3 = apply_map(source, R, Z + h, Phi, 1);
    std::pair<double, double> point_4 = apply_map(source, R, Z - h, Phi, 1);

    // Compute the jacobian
    jacobian[0][0] = (point_1.first - point_2.first) / (2 * h);
    jacobian[0][1] = (point_3.first - point_4.first) / (2 * h);
    jacobian[1][0] = (point_1.second - point_2.second) / (2 * h);
    jacobian[1][1] = (point_3.second - point_4.second) / (2 * h);
}

// Iteratively find the closest 1 period fixed point from the initial guess
bool x_point(maglit &source, double &R, double &Z, double &Phi, double tol, int max_iter) {
    for (int iter = 0; iter < max_iter; ++iter) {
        // std::cout << "iter: " << iter << std::endl;
        //  initialize jacobian and identity
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
        std::pair<double, double> R1_Z1 = apply_map(source, R, Z, Phi, 1);
        double R1 = R1_Z1.first;
        double Z1 = R1_Z1.second;

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

        std::cout << "R: " << R << " Z: " << Z << std::endl;
        // check if iteration limit is reached
        if (iter == max_iter - 1) {
            std::cout << "Newton's method reached maximum number of iterations" << std::endl;
            return false;
        }
    }
    return false;
}

void primarySegment(maglit &source, double R_xpoint, double Z_xpoint, double Phi, double epsilon, double num_points, std::vector<double> &R_vec, std::vector<double> &Z_vec) {
    // Evaluate the jacobian at the x-point
    double jacobian[2][2];
    eval_jacobian(source, R_xpoint, Z_xpoint, Phi, 1e-5, jacobian);
    std::cout << "Jacobian at the x-point: " << std::endl;
    std::cout << "[" << jacobian[0][0] << " " << jacobian[0][1] << "]" << std::endl;
    std::cout << "[" << jacobian[1][0] << " " << jacobian[1][1] << "]\n"
              << std::endl;

    // Convert the jacobian to an armadillo matrix
    arma::mat jacobian_arma = {{jacobian[0][0], jacobian[0][1]}, {jacobian[1][0], jacobian[1][1]}};

    // Compute the eigenvalues and eigenvectors of the jacobian
    arma::cx_vec eigval;
    arma::cx_mat eigvec;
    arma::eig_gen(eigval, eigvec, jacobian_arma);

    // Print the eigenvalues and eigenvectors
    std::cout << "Eigenvalues: " << std::endl;
    std::cout << eigval[0].real() << " " << eigval[1].real() << "\n"
              << std::endl;

    std::cout << "Eigenvectors: " << std::endl;
    std::cout << "[" << eigvec.col(0)[0].real() << " " << eigvec.col(0)[1].real() << "]" << std::endl;
    std::cout << "[" << eigvec.col(1)[0].real() << " " << eigvec.col(1)[1].real() << "]\n"
              << std::endl;

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
    std::cout << "[" << v_R << " " << v_Z << "]\n"
              << std::endl;

    // Calculate the new point (R_new, Z_new)
    double R_new = R_xpoint + epsilon * v_R;
    double Z_new = Z_xpoint + epsilon * v_Z;

    std::cout << "Primary segment first point: " << std::endl;
    std::cout << "R_new: " << R_new << " Z_new: " << Z_new << "\n"
              << std::endl;

    // Apply the map M to the new point: M(R_new, Z_new) = (MR, MZ)
    std::pair<double, double> M_rz = apply_map(source, R_new, Z_new, Phi, 1);
    double MR = M_rz.first;
    double MZ = M_rz.second;

    std::cout << "Primary segment last point: " << std::endl;
    std::cout << "MR: " << MR << " MZ: " << MZ << "\n"
              << std::endl;

    // Interpolate points between (R_new, Z_new) and (MR, MZ)
    double step_R = (MR - R_new) / (num_points - 1);
    double step_Z = (MZ - Z_new) / (num_points - 1);

    for (int i = 0; i < num_points; ++i) {
        R_vec[i] = R_new + i * step_R;
        Z_vec[i] = Z_new + i * step_Z;
    }
};

// Compute distance between two points
double computeDistance(double R1, double Z1, double R2, double Z2) {
    return std::sqrt((R1 - R2) * (R1 - R2) + (Z1 - Z2) * (Z1 - Z2));
}

// Compute angle between two vectors
double computeAngle(double R0, double Z0, double R1, double Z1, double R2, double Z2) {
    // compute vector from (R0, Z0) to (R1, Z1)
    double v1_R = R1 - R0;
    double v1_Z = Z1 - Z0;

    // compute vector from (R1, Z1) to (R2, Z2)
    double v2_R = R2 - R1;
    double v2_Z = Z2 - Z1;

    // compute angle between v1 and v2
    double dot_product = v1_R * v2_R + v1_Z * v2_Z;
    double norm_v1 = std::sqrt(v1_R * v1_R + v1_Z * v1_Z);
    double norm_v2 = std::sqrt(v2_R * v2_R + v2_Z * v2_Z);

    return std::acos(dot_product / (norm_v1 * norm_v2));
}

// Function to insert a new point in the vector
void insertPoint(std::vector<double> &R_vec, std::vector<double> &Z_vec, int i) {
    double R_new = (R_vec[i - 1] + R_vec[i]) / 2.0;
    double Z_new = (Z_vec[i - 1] + Z_vec[i]) / 2.0;
    R_vec.insert(R_vec.begin() + i, R_new);
    Z_vec.insert(Z_vec.begin() + i, Z_new);
}

// Print a progress bar
void progressBar(float progress) {
    int barWidth = 50;

    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

void newSegment(maglit &tracer, std::vector<double> &R_PrevSeg, std::vector<double> &Z_PrevSeg, std::vector<double> &R_NewSeg, std::vector<double> &Z_NewSeg, double Phi, int nSeg, double l_lim, double theta_lim) {
    size_t j = 1;
    while (j < R_PrevSeg.size() - 1) {
        std::pair<double, double> x_i = apply_map(tracer, R_PrevSeg[j - 1], Z_PrevSeg[j - 1], Phi, nSeg);
        std::pair<double, double> x_j = apply_map(tracer, R_PrevSeg[j], Z_PrevSeg[j], Phi, nSeg);
        std::pair<double, double> x_k = apply_map(tracer, R_PrevSeg[j + 1], Z_PrevSeg[j + 1], Phi, nSeg);

        // double l_p = computeDistance(R_PrevSeg[j - 1], Z_PrevSeg[j - 1], R_PrevSeg[j], Z_PrevSeg[j]);
        double l_i = computeDistance(x_i.first, x_i.second, x_j.first, x_j.second);
        double l_ii = computeDistance(x_j.first, x_j.second, x_k.first, x_k.second);
        double l_theta = computeAngle(x_i.first, x_i.second, x_j.first, x_j.second, x_k.first, x_k.second);

        double L_lim;
        if (l_theta > theta_lim) {
            L_lim = l_lim * 0.1;
        } else {
            L_lim = l_lim;
        }

        if (l_i > L_lim) {
            insertPoint(R_PrevSeg, Z_PrevSeg, j);
            if (j != 1) {
                j = j - 1;
            }
        } else if (l_ii > L_lim) {
            insertPoint(R_PrevSeg, Z_PrevSeg, j + 1);
        } else {
            // std::cout << "passed with:" << l_theta << std::endl;
            R_NewSeg.push_back(x_i.first);
            Z_NewSeg.push_back(x_i.second);
            // std::cout << x_i.first << " " << x_i.second << std::endl;
            if (j == R_PrevSeg.size() - 2) {
                R_NewSeg.push_back(x_j.first);
                Z_NewSeg.push_back(x_j.second);
                R_NewSeg.push_back(x_k.first);
                Z_NewSeg.push_back(x_k.second);
            }
            j = j + 1;
        }
        // std::cout << x_i.first << " " << x_i.second << std::endl;
    }
}

// test the x_point function
int main() {
    char source_path[] = "/home/jfernandes/Software/m3dc1_data/cp_coils/n06/C1.h5";
    maglit tracer(source_path, FIO_M3DC1_SOURCE);
    tracer.inverse_map(false);
    double R_xpoint = 0.49800405;
    double Z_xpoint = -0.21860696;
    double Phi = 0;

    bool found_xp = x_point(tracer, R_xpoint, Z_xpoint, Phi, 1e-12, 100);
    if (!found_xp) {
        std::cerr << "Failed to find the x-point" << std::endl;
        return 1;
    }

    std::cout << std::fixed << std::setprecision(16);
    std::cout << std::endl;
    std::cout << "x-point found at: " << std::endl;
    std::cout << "R: " << R_xpoint << " Z: " << Z_xpoint << "\n"
              << std::endl;
    // std::cout << std::defaultfloat << std::setprecision(6);

    // Define the small distance epsilon
    double epsilon = 1e-4;
    int num_points = 10;
    std::vector<double> R_primeseg(num_points);
    std::vector<double> Z_primeseg(num_points);
    primarySegment(tracer, R_xpoint, Z_xpoint, Phi, epsilon, num_points, R_primeseg, Z_primeseg);

    std::vector<double> R_newseg, Z_newseg;

    double l_lim = 0.01;
    double theta_lim = 0.17;
    double nSeg = 6;

    for (int i = 1; i < nSeg; ++i) {
        newSegment(tracer, R_primeseg, Z_primeseg, R_newseg, Z_newseg, Phi, i, l_lim, theta_lim);
        std::ofstream file2("cp_n06_S.dat");
        file2 << std::fixed << std::setprecision(16);
        if (file2.is_open()) {
            for (size_t i = 0; i < R_newseg.size(); ++i) {
                file2 << R_newseg[i] << " " << Z_newseg[i] << std::endl;
            }
            file2.close();
        }

        progressBar((float)i / nSeg);
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }

    std::cout << "\nDone! \n"
              << std::endl;

    return 0;
}