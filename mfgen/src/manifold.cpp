#include "manifold.h"

// Constructor
manifold::manifold(const char *source_path, const int timeslice, double phi, int stability, double epsilon)
    : stability(stability), phi(phi), epsilon(epsilon), tracer(source_path, FIO_M3DC1_SOURCE, timeslice) {

    // Set the inverse map based on the stability
    switch (stability) {
    case 0:
        tracer.inverse_map(false);
        this->s_factor = 1;
        break;
    case 1:
        tracer.inverse_map(true);
        this->s_factor = -1;
        break;
    }
}

// Set warning flag
void manifold::setWarnings() {
    this->warnings = true;
    tracer.set_warnings();
}

// Apply map to a given point returning a Point (R, Z)
point manifold::apply_map(double R, double Z, double Phi, int nTurns) {

    // Reset the source and allocate hint
    int status = SODE_CONTINUE_GOOD_STEP;
    tracer.reset();
    tracer.alloc_hint();
    double phi_max = Phi + nTurns * 2 * M_PI;

    // Apply the map
    do {
        status = tracer.step(R, Z, Phi, phi_max, 0);
    } while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP);
    tracer.clear_hint();

    // Check if the map was successful
    if (status == SODE_SUCCESS_TIME) {
        return {R, Z};
    } else {
        if (this->warnings) {
            std::cerr << "Failed to apply map at (R, Z) = (" << R << ", " << Z << ")" << std::endl;
        }
        return {std::nan(""), std::nan("")};
    }
}

// Evaluate the jacobian of the map at a given point
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
        if (this->warnings) {
            std::cout << std::setprecision(16) << "R: " << R << " Z: " << Z << std::endl;
        }
    }

    // Maximum number of iterations reached
    std::cout << "Newton's method reached maximum number of iterations" << std::endl;
    return false;
}

// Find the pivot point for the primary segment
point manifold::pivot() {
    double R_xPoint = this->xPoint.R;
    double Z_xPoint = this->xPoint.Z;
    double jacobian[2][2];

    // Evaluate the jacobian at the x-point
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

// Compute the primary segment
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

    // Compute the primary segment
    for (size_t i = 0; i <= num_points; i++) {
        double t = (i) / static_cast<double>(num_points);
        double R = (1 - t) * pivot.R + t * T_pivot.R;
        double Z = (1 - t) * pivot.Z + t * T_pivot.Z;
        segment.push_back({R, Z});

        // std::cout << "R: " << segment[i].R << " Z: " << segment[i].Z << std::endl;
    }
}

// Compute distance between two points
double manifold::computeDistance(double R1, double Z1, double R2, double Z2) {
    // Compute vector from (R1, Z1) to (R2, Z2)
    double dR = R1 - R2;
    double dZ = Z1 - Z2;
    // Return the verctor's magnitude
    return std::sqrt(dR * dR + dZ * dZ);
}

// Compute angle between two vectors
double manifold::computeAngle(double R0, double Z0, double R1, double Z1, double R2, double Z2) {
    // Compute vector from (R0, Z0) to (R1, Z1)
    double v1_R = R1 - R0;
    double v1_Z = Z1 - Z0;

    // Compute vector from (R1, Z1) to (R2, Z2)
    double v2_R = R2 - R1;
    double v2_Z = Z2 - Z1;

    // Compute dot product of v1 and v2
    double dot_product = v1_R * v2_R + v1_Z * v2_Z;

    // Compute norms (magnitudes) of v1 and v2
    double norm_v1 = std::sqrt(v1_R * v1_R + v1_Z * v1_Z);
    double norm_v2 = std::sqrt(v2_R * v2_R + v2_Z * v2_Z);

    // Check for zero-length vectors to avoid division by zero
    if (norm_v1 == 0.0 || norm_v2 == 0.0) {
        return 0.0; // Angle is 0 if any vector is zero (e.g., points are the same)
    }

    // Compute cosine of the angle and clamp to avoid precision issues
    double cos_theta = dot_product / (norm_v1 * norm_v2);

    // Clamp value to the valid range [-1, 1] to avoid floating-point errors
    if (cos_theta > 1.0)
        cos_theta = 1.0;
    if (cos_theta < -1.0)
        cos_theta = -1.0;

    // Return the angle in radians
    return std::acos(cos_theta);
}

// Insert a new point in the vector by linear interpolation
void manifold::insertPoint(std::vector<point> &segment, size_t index) {
    // Check if the index is valid
    if (index == 0 || index >= segment.size()) {
        std::cerr << "Error: Invalid index for insertPoint." << std::endl;
        return;
    }

    // Compute the difference between the two points
    double dR = segment[index].R - segment[index - 1].R;
    double dZ = segment[index].Z - segment[index - 1].Z;

    // Only insert a point if the distance between the two points is larger than the precision limit
    if (std::sqrt(dR * dR + dZ * dZ) > this->precision_limit) {
        point new_point;
        new_point.R = (segment[index - 1].R + segment[index].R) / 2;
        new_point.Z = (segment[index - 1].Z + segment[index].Z) / 2;
        segment.insert(segment.begin() + index, new_point);
    } else {
        // Set the overlap flag to true if the precision limit is reached
        this->overlap = true;
        if (this->warnings) {
            std::cerr << "Warning: Skipping insertion due to floating point precision limits." << std::endl;
        }
    }
}

// Compute a refined new segment from a previous segment
void manifold::newSegment(std::vector<point> &prev_seg, std::vector<point> &new_seg, double Phi, int nSeg, double l_lim, double theta_lim) {

    size_t j = 1;                     // Start at the second point
    double theta_lim_aux = theta_lim; // Copy of the angle limit
    int insertion_count = 0;          // Counter for insertions in the current segment

    while (j < prev_seg.size() - 1) {
        // Apply map to the points
        point x_i = this->apply_map(prev_seg[j - 1].R, prev_seg[j - 1].Z, Phi, nSeg);
        point x_j = this->apply_map(prev_seg[j].R, prev_seg[j].Z, Phi, nSeg);
        point x_k = this->apply_map(prev_seg[j + 1].R, prev_seg[j + 1].Z, Phi, nSeg);

        // Compute distances and angle to check if the points are too far apart
        double l_i = this->computeDistance(x_i.R, x_i.Z, x_j.R, x_j.Z);
        double l_ii = this->computeDistance(x_j.R, x_j.Z, x_k.R, x_k.Z);
        double l_theta = this->computeAngle(x_i.R, x_i.Z, x_j.R, x_j.Z, x_k.R, x_k.Z);

        // if not refining angle, reset the angle limit
        if (!this->refining_angle) {
            theta_lim_aux = theta_lim;
        }

        // Check if angle exceeds threshold
        if (l_theta > theta_lim_aux) {
            if (l_i > l_ii) {
                this->insertPoint(prev_seg, j);
                insertion_count++;
                if (j != 1) {
                    j = j - 1;
                }
            } else {
                this->insertPoint(prev_seg, j + 1);
                insertion_count++;
            }
            this->refining_angle = true;
        }
        // Check if distance l_i exceed threshold
        else if (l_i > l_lim) {
            this->insertPoint(prev_seg, j);
            insertion_count++;
            if (j != 1) {
                j = j - 1;
            }
            this->refining_angle = false;
        }
        // Check if distance l_ii exceed threshold
        else if (l_ii > l_lim) {
            this->insertPoint(prev_seg, j + 1);
            insertion_count++;
            this->refining_angle = false;
        }
        // If conditions are met, add the point to the new segment
        else {
            new_seg.push_back(x_i);
            if (j == prev_seg.size() - 2) {
                new_seg.push_back(x_j);
                new_seg.push_back(x_k);
            }
            j = j + 1;
            insertion_count = 0;
        }
        // Relax the angle limit if the precision limit is reached
        if (this->overlap) {
            theta_lim_aux *= 1.5;
            this->overlap = false;
        }
        // Stop refinement if the maximum number of insertions is reached
        if (insertion_count >= this->max_insertions) {
            if (this->warnings) {
                std::cerr << "Warning: Maximum number of insertions reached. Stopping refinement at this segment." << std::endl;
            }
            break;
        }
    }
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