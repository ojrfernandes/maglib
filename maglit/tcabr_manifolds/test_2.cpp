#include "../maglit.h"
#include <armadillo>
#include <chrono>
#include <iomanip>
#include <thread>
#include <utility>

struct point {
    long double R;
    long double Z;
};

class manifold {
  public:
    manifold(const char *source_path, const int timeslice, double phi, int stability);
    bool find_xPoint(double rGuess, double zGuess);
    void primarySegment(std::vector<point> &segment, size_t num_points);
    void newSegment(std::vector<point> &prev_seg, std::vector<point> &new_seg, double Phi, int nSeg, double l_lim, double theta_lim);
    void progressBar(float progress);

    point xPoint; // x-point

  private:
    void eval_jacobian(double R, double Z, double Phi, double h, double jacobian[2][2]);
    void insertPoint(std::vector<point> &segment, size_t index);
    double computeDistance(double R1, double Z1, double R2, double Z2);
    double computeAngle(double R0, double Z0, double R1, double Z1, double R2, double Z2);
    point apply_map(double R, double Z, double Phi, int nTurns);
    point pivot();

    int stability;    // 0: forward, 1: backward
    int s_factor = 1; // sign factor for the manifold stability
    int max_iter;     // maximum number of iterations for Newton's method
    double phi;       // toroidal angle
    double tol;       // tolerance for Newton's method
    double epsilon;   // distance from the x-point
    double h;         // step size for numerical differentiation

    maglit tracer;
};

manifold::manifold(const char *source_path, const int timeslice, double phi, int stability)
    : stability(stability), phi(phi), tracer(source_path, FIO_M3DC1_SOURCE, timeslice) {
    this->max_iter = 100;
    this->tol = 1e-14;
    this->h = 1e-5;
    this->epsilon = 1e-6;

    switch (stability) {
    case 0:
        tracer.inverse_map(false);
        break;
    case 1:
        tracer.inverse_map(true);
        this->s_factor = -1;
        break;
    }
}

// Apply map to a given point returning a Point (R, Z)
point manifold::apply_map(double R, double Z, double Phi, int nTurns) {

    // Reset the source and allocate hint
    int status = SODE_CONTINUE_GOOD_STEP;
    tracer.reset();
    tracer.alloc_hint();
    double phi_max = nTurns * 2 * M_PI;

    // Apply the map
    do {
        status = tracer.step(R, Z, Phi, phi_max, 0);
    } while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP);
    tracer.clear_hint();

    // Check if the map was successful
    if (status == SODE_SUCCESS_TIME) {
        return {R, Z};
    } else {
        std::cerr << "Failed to apply map" << std::endl;
        return {std::nan(""), std::nan("")};
    }
}

void manifold::eval_jacobian(double R, double Z, double Phi, double h, double jacobian[2][2]) {

    // Apply maps to the points (R + h, Z), (R - h, Z), (R, Z + h), (R, Z - h)
    point P1 = this->apply_map(R + h, Z, Phi, 1);
    point P2 = this->apply_map(R - h, Z, Phi, 1);
    point P3 = this->apply_map(R, Z + h, Phi, 1);
    point P4 = this->apply_map(R, Z - h, Phi, 1);

    // Compute the jacobian
    jacobian[0][0] = (P1.R - P2.R) / (2 * h);
    jacobian[0][1] = (P3.R - P4.R) / (2 * h);
    jacobian[1][0] = (P1.Z - P2.Z) / (2 * h);
    jacobian[1][1] = (P3.Z - P4.Z) / (2 * h);
}

// Iteratively find the closest 1 period fixed point from the initial guess
bool manifold::find_xPoint(double R, double Z) {

    // Iterate to find the x-point
    for (int iter = 0; iter < this->max_iter; ++iter) {

        //  Initialize jacobian and identity
        double identity[2][2] = {{1, 0}, {0, 1}};
        double jacobian[2][2];
        this->eval_jacobian(R, Z, this->phi, this->h, jacobian);

        // Computer (identity - jacobian)
        double diff[2][2];
        diff[0][0] = identity[0][0] - jacobian[0][0];
        diff[0][1] = identity[0][1] - jacobian[0][1];
        diff[1][0] = identity[1][0] - jacobian[1][0];
        diff[1][1] = identity[1][1] - jacobian[1][1];

        // Check if (identity - jacobian) is singular
        double det = diff[0][0] * diff[1][1] - diff[0][1] * diff[1][0];
        if (std::abs(det) < 1e-12) {
            std::cout << "(I - J) is a singular matrix" << std::endl;
            return false;
        }

        // Compute (identity - jacobian)^(-1)
        double inv_diff[2][2];
        inv_diff[0][0] = diff[1][1] / det;
        inv_diff[0][1] = -diff[0][1] / det;
        inv_diff[1][0] = -diff[1][0] / det;
        inv_diff[1][1] = diff[0][0] / det;

        // Compute T(Rn, Zn)
        point T_Pn = this->apply_map(R, Z, this->phi, 1);

        // Update the position
        double delta_R = T_Pn.R - (jacobian[0][0] * R + jacobian[0][1] * Z);
        double delta_Z = T_Pn.Z - (jacobian[1][0] * R + jacobian[1][1] * Z);
        R = inv_diff[0][0] * delta_R + inv_diff[0][1] * delta_Z;
        Z = inv_diff[1][0] * delta_R + inv_diff[1][1] * delta_Z;

        // Check if the increment is small enough
        double dR = R - T_Pn.R;
        double dZ = Z - T_Pn.Z;
        if (std::sqrt(dR * dR + dZ * dZ) < this->tol) {
            this->xPoint = {R, Z};
            return true;
        }

        std::cout << std::setprecision(14) << "R: " << R << " Z: " << Z << std::endl;
    }

    // Maximum number of iterations reached
    std::cout << "Newton's method reached maximum number of iterations" << std::endl;
    return false;
}

point manifold::pivot() {
    double R_xPoint = this->xPoint.R;
    double Z_xPoint = this->xPoint.Z;
    // Evaluate the jacobian at the x-point
    double jacobian[2][2];
    this->eval_jacobian(R_xPoint, Z_xPoint, this->phi, this->h, jacobian);
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
    double v_Z = selected_eigvec[1].real() * this->s_factor;
    std::cout << "Selected eigenvector (|lambda| > 1): " << std::endl;
    std::cout << "[" << v_R << " " << v_Z << "]\n"
              << std::endl;

    // Calculate the new point (R_new, Z_new)
    double R_pivot = R_xPoint + this->epsilon * v_R;
    double Z_pivot = Z_xPoint + this->epsilon * v_Z;

    return {R_pivot, Z_pivot};
}

void manifold::primarySegment(std::vector<point> &segment, size_t num_points) {
    // Find the primary segment's first point
    point pivot = this->pivot();
    std::cout << "Primary segment's first point: " << "\n"
              << "R: " << pivot.R << " Z: " << pivot.Z << "\n"
              << std::endl;

    // Find the primary segment's last point
    point T_pivot = this->apply_map(pivot.R, pivot.Z, this->phi, 1);
    std::cout << "Primary segment's last point: " << std::endl;
    std::cout << "R: " << T_pivot.R << " Z: " << T_pivot.Z << "\n"
              << std::endl;

    for (size_t i = 0; i <= num_points; i++) {
        double t = (i) / static_cast<double>(num_points);
        double R = (1 - t) * pivot.R + t * T_pivot.R;
        double Z = (1 - t) * pivot.Z + t * T_pivot.Z;
        segment.push_back({R, Z});
        // std::cout << "R: " << segment[i].R << " Z: " << segment[i].Z << std::endl;
    }
}

// Insert a new point in the vector by linear interpolation
void manifold::insertPoint(std::vector<point> &segment, size_t index) {
    point new_point;
    new_point.R = (segment[index - 1].R + segment[index].R) / 2;
    new_point.Z = (segment[index - 1].Z + segment[index].Z) / 2;
    segment.insert(segment.begin() + index, new_point);
}

void manifold::newSegment(std::vector<point> &prev_seg, std::vector<point> &new_seg, double Phi, int nSeg, double l_lim, double theta_lim) {
    size_t j = 1;
    while (j < prev_seg.size() - 1) {

        // Apply map to the points
        point x_i = this->apply_map(prev_seg[j - 1].R, prev_seg[j - 1].Z, Phi, nSeg);
        point x_j = this->apply_map(prev_seg[j].R, prev_seg[j].Z, Phi, nSeg);
        point x_k = this->apply_map(prev_seg[j + 1].R, prev_seg[j + 1].Z, Phi, nSeg);

        // Compute distances and angle to check if the points are too far apart
        double l_i = this->computeDistance(x_i.R, x_i.Z, x_j.R, x_j.Z);
        double l_ii = this->computeDistance(x_j.R, x_j.Z, x_k.R, x_k.Z);
        double l_theta = this->computeAngle(x_i.R, x_i.Z, x_j.R, x_j.Z, x_k.R, x_k.Z);

        // Check if the angle is too large
        double L_lim;
        if (l_theta > theta_lim) {
            L_lim = l_lim * 0.1;
        } else {
            L_lim = l_lim;
        }

        // if (l_theta > theta_lim) {
        //     if (l_i > l_ii) {
        //         this->insertPoint(prev_seg, j);
        //         if (j != 1) {
        //             j = j - 1;
        //         }
        //     } else {
        //         this->insertPoint(prev_seg, j + 1);
        //     }
        // }

        // Check if the points are too far apart
        if (l_i > L_lim) {
            this->insertPoint(prev_seg, j);
            if (j != 1) {
                j = j - 1;
            }
        } else if (l_ii > L_lim) {
            this->insertPoint(prev_seg, j + 1);
        } else {
            new_seg.push_back(x_i);
            if (j == prev_seg.size() - 2) {
                new_seg.push_back(x_j);
                new_seg.push_back(x_k);
            }
            j = j + 1;
        }
        // std::cout << "  x_j: " << x_j.R << " " << x_j.Z << std::endl;
    }
}

// Compute distance between two points
double manifold::computeDistance(double R1, double Z1, double R2, double Z2) {
    return std::sqrt((R1 - R2) * (R1 - R2) + (Z1 - Z2) * (Z1 - Z2));
}

// Compute angle between two vectors
double manifold::computeAngle(double R0, double Z0, double R1, double Z1, double R2, double Z2) {
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

// Print a progress bar
void manifold::progressBar(float progress) {
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

int main() {
    char source_path[] = "/home/jfernandes/Software/maglib_local/maglit/our_m3dc1_data/i_coils/n03/C1.h5";
    int timeslice = 1;
    int stability = 1;
    double Phi = 0;

    manifold manifold(source_path, timeslice, Phi, stability);

    double R_xPoint = 0.46950445;
    double Z_xPoint = -0.21327174;
    bool found_xp = manifold.find_xPoint(R_xPoint, Z_xPoint);
    if (!found_xp) {
        std::cerr << "Failed to find the x-point" << std::endl;
        return 1;
    }

    std::cout << std::fixed << std::setprecision(16) << "\n"
              << "x-point found at: " << "\n"
              << "R: " << R_xPoint << " Z: " << Z_xPoint << "\n"
              << std::endl;

    size_t num_points = 10;
    std::vector<point> primary_segment;

    manifold.primarySegment(primary_segment, num_points);

    double l_lim = 0.02;
    double theta_lim = 0.17;
    double nSeg = 10;

    std::vector<point> new_segment;
    std::vector<point> prev_segment = primary_segment;

    for (int i = 1; i < nSeg; ++i) {
        // create new segment
        manifold.newSegment(prev_segment, new_segment, Phi, i, l_lim, theta_lim);

        // append new segment to the output file
        std::ofstream file2("i_n03_U.dat", std::ios::app);
        file2 << std::fixed << std::setprecision(16);
        if (file2.is_open()) {
            for (size_t i = 0; i < new_segment.size(); ++i) {
                file2 << new_segment[i].R << " " << new_segment[i].Z << std::endl;
            }
            file2.close();
        } else {
            std::cerr << "Unable to open file" << std::endl;
            return 1;
        }

        // empty new segment
        // prev_segment = new_segment;
        new_segment.clear();

        // progress bar
        manifold.progressBar((float)i / nSeg);
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }

    std::cout << "\nDone! \n"
              << std::endl;

    return 0;
}