#include "manifold.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

manifold::manifold(maglit &tracer, double phi, int stability)
    : stability(stability), phi(phi), tracer(tracer) {
    switch (stability) {
    case 0:
        tracer.inverse_map(false);
        s_factor = 1;
        break;
    case 1:
        tracer.inverse_map(true);
        s_factor = -1;
        break;
    }
}

void manifold::setVerbose() {
    verbose = true;
    tracer.set_warnings();
}

void manifold::configure(double epsilon, double h, double tol, int max_iter,
                          double precision_limit, int max_insertions) {
    this->epsilon         = epsilon;
    this->h               = h;
    this->tol             = tol;
    this->max_iter        = max_iter;
    this->precision_limit = precision_limit;
    this->max_insertions  = max_insertions;
}

point manifold::apply_map(double R, double Z, double Phi, int nTurns) {
    int    status  = SODE_CONTINUE_GOOD_STEP;
    double phi_max = Phi + nTurns * 2 * M_PI;
    tracer.reset();
    do {
        status = tracer.step(R, Z, Phi, phi_max, 0);
    } while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP);

    if (status == SODE_SUCCESS_TIME)
        return {R, Z};
    if (verbose)
        std::cerr << "Failed to apply map at (R, Z) = (" << R << ", " << Z << ")\n";
    return {std::nan(""), std::nan("")};
}

void manifold::eval_jacobian(double R, double Z, double Phi, double h, double jacobian[2][2]) {
    point P1 = apply_map(R + h, Z,     Phi, 1);
    point P2 = apply_map(R - h, Z,     Phi, 1);
    point P3 = apply_map(R,     Z + h, Phi, 1);
    point P4 = apply_map(R,     Z - h, Phi, 1);

    jacobian[0][0] = (P1.R - P2.R) / (2 * h);
    jacobian[0][1] = (P3.R - P4.R) / (2 * h);
    jacobian[1][0] = (P1.Z - P2.Z) / (2 * h);
    jacobian[1][1] = (P3.Z - P4.Z) / (2 * h);
}

bool manifold::find_xPoint(double R, double Z) {
    for (int iter = 0; iter < max_iter; ++iter) {
        double jacobian[2][2];
        eval_jacobian(R, Z, phi, h, jacobian);

        double diff[2][2] = {
            {1.0 - jacobian[0][0], -jacobian[0][1]},
            {-jacobian[1][0],       1.0 - jacobian[1][1]}
        };

        double det = diff[0][0] * diff[1][1] - diff[0][1] * diff[1][0];
        if (std::abs(det) < 1e-12) {
            if (verbose) std::cout << "(I - J) is a singular matrix\n";
            return false;
        }

        double inv[2][2] = {
            { diff[1][1] / det, -diff[0][1] / det},
            {-diff[1][0] / det,  diff[0][0] / det}
        };

        point T_Pn = apply_map(R, Z, phi, 1);

        if (verbose) {
            std::cout << "Iteration " << iter + 1 << "\n"
                      << std::setprecision(16) << " R = " << R << ", Z = " << Z << "\n";
        }

        double dR = T_Pn.R - (jacobian[0][0] * R + jacobian[0][1] * Z);
        double dZ = T_Pn.Z - (jacobian[1][0] * R + jacobian[1][1] * Z);
        R = inv[0][0] * dR + inv[0][1] * dZ;
        Z = inv[1][0] * dR + inv[1][1] * dZ;

        double err = std::sqrt((R - T_Pn.R) * (R - T_Pn.R) + (Z - T_Pn.Z) * (Z - T_Pn.Z));
        if (verbose) std::cout << "error = " << err << "\n\n";

        if (err < tol) {
            xPoint = {R, Z};
            return true;
        }
    }

    if (verbose) std::cout << "Newton's method reached maximum number of iterations\n";
    return false;
}

point manifold::pivot() {
    double jacobian[2][2];
    eval_jacobian(xPoint.R, xPoint.Z, phi, h, jacobian);

    if (verbose) {
        std::cout << "Jacobian at the x-point:\n"
                  << "[" << jacobian[0][0] << " " << jacobian[0][1] << "]\n"
                  << "[" << jacobian[1][0] << " " << jacobian[1][1] << "]\n\n";
    }

    // Analytic 2x2 eigendecomposition. The Poincaré map Jacobian at a hyperbolic
    // X-point has real eigenvalues (discriminant > 0).
    double tr   = jacobian[0][0] + jacobian[1][1];
    double det  = jacobian[0][0] * jacobian[1][1] - jacobian[0][1] * jacobian[1][0];
    double sq   = std::sqrt(tr * tr - 4.0 * det);
    double lam1 = (tr + sq) / 2.0;
    double lam2 = (tr - sq) / 2.0;

    if (verbose) std::cout << "Eigenvalues: " << lam1 << " " << lam2 << "\n\n";

    double lam = (std::abs(lam1) > 1.0) ? lam1 : lam2;

    // Eigenvector from the null space of (J - λI). Use the row whose off-diagonal
    // entry is larger for numerical stability.
    double v_R, v_Z;
    if (std::abs(jacobian[0][1]) >= std::abs(jacobian[1][0])) {
        v_R = jacobian[0][1];
        v_Z = lam - jacobian[0][0];
    } else {
        v_R = lam - jacobian[1][1];
        v_Z = jacobian[1][0];
    }

    double norm = std::sqrt(v_R * v_R + v_Z * v_Z);
    v_R /= norm;
    v_Z  = (v_Z / norm) * s_factor;

    if (verbose)
        std::cout << "Selected eigenvector (|lambda| > 1): [" << v_R << " " << v_Z << "]\n\n";

    return {xPoint.R + epsilon * v_R, xPoint.Z + epsilon * v_Z};
}

std::vector<point> manifold::primarySegment(size_t n_intervals) {
    point pivot_pt = pivot();
    point T_pivot  = apply_map(pivot_pt.R, pivot_pt.Z, phi, 1);

    if (verbose) {
        std::cout << "Primary segment's first point: R=" << pivot_pt.R << " Z=" << pivot_pt.Z << "\n"
                  << "Primary segment's last point:  R=" << T_pivot.R  << " Z=" << T_pivot.Z  << "\n\n";
    }

    std::vector<point> segment;
    segment.reserve(n_intervals + 1);
    for (size_t i = 0; i <= n_intervals; ++i) {
        double t = static_cast<double>(i) / static_cast<double>(n_intervals);
        segment.push_back(pivot_pt * (1.0 - t) + T_pivot * t);
    }
    outputData.push_back(segment);
    return segment;
}

double manifold::computeDistance(point a, point b) {
    point d = b - a;
    return std::sqrt(d.R * d.R + d.Z * d.Z);
}

double manifold::computeAngle(point a, point b, point c) {
    point  v1 = b - a;
    point  v2 = c - b;
    double n1 = std::sqrt(v1.R * v1.R + v1.Z * v1.Z);
    double n2 = std::sqrt(v2.R * v2.R + v2.Z * v2.Z);
    if (n1 == 0.0 || n2 == 0.0) return 0.0;
    double cos_theta = std::clamp((v1.R * v2.R + v1.Z * v2.Z) / (n1 * n2), -1.0, 1.0);
    return std::acos(cos_theta);
}

void manifold::insertPoint(std::vector<point> &segment, size_t index, bool &overlap) {
    if (index == 0 || index >= segment.size()) {
        std::cerr << "Error: Invalid index for insertPoint.\n";
        return;
    }
    if (computeDistance(segment[index - 1], segment[index]) > precision_limit) {
        segment.insert(segment.begin() + index,
                       segment[index - 1] * 0.5 + segment[index] * 0.5);
    } else {
        overlap = true;
        if (verbose)
            std::cerr << "Warning: Skipping insertion due to floating point precision limits.\n";
    }
}

void manifold::insertPoint(std::vector<point> &segment, interpolantArc &arc, bool &overlap) {
    if (computeDistance(arc.x0, arc.x1) > precision_limit) {
        point new_point = arc.evalNewPoint(0.5);

        if (verbose) {
            std::cout << "x0: " << arc.x0.R << " " << arc.x0.Z << "\n"
                      << "x1: " << arc.x1.R << " " << arc.x1.Z << "\n"
                      << "New point: " << new_point.R << " " << new_point.Z << "\n";
            point img = apply_map(new_point.R, new_point.Z, phi, 1);
            if (std::isnan(img.R) || std::isnan(img.Z)) {
                std::cerr << "Warning: NaN in insertPoint R=" << new_point.R
                          << " Z=" << new_point.Z << "\n";
                return;
            }
        }

        if (static_cast<size_t>(arc.i0) < segment.size() &&
            static_cast<size_t>(arc.i1) < segment.size()) {
            segment.insert(segment.begin() + arc.i1, new_point);
            if (verbose)
                std::cout << "Inserted new point at index: " << arc.i1 << "\n\n";
        } else if (verbose) {
            std::cerr << "Warning: Invalid indices for insertion in insertPoint.\n";
        }
    } else {
        overlap = true;
        if (verbose)
            std::cerr << "Warning: Skipping insertion due to floating point precision limits.\n";
    }
}

std::vector<point> manifold::newSegment(std::vector<point> &prev_seg, double Phi,
                                         double l_lim, double theta_lim) {
    if (prev_seg.size() < 3) {
        std::cerr << "Error: Previous segment must have at least three points.\n";
        return {};
    }

    std::vector<point> new_seg = {
        prev_seg[prev_seg.size() - 3],
        prev_seg[prev_seg.size() - 2],
        prev_seg[prev_seg.size() - 1]
    };

    theta_lim *= M_PI / 180.0;
    double theta_lim_aux  = theta_lim;
    int    insertion_count = 0;
    bool   refining_angle  = false;
    bool   overlap         = false;

    std::vector<interpolantArc> arcs = buildInterpolants(prev_seg);
    size_t j = 1;

    while (j < arcs.size()) {
        point x_i = apply_map(arcs[j - 1].x0.R, arcs[j - 1].x0.Z, Phi, 1);
        point x_j = apply_map(arcs[j].x0.R,     arcs[j].x0.Z,     Phi, 1);
        point x_k = apply_map(arcs[j].x1.R,     arcs[j].x1.Z,     Phi, 1);

        double l_i     = computeDistance(x_i, x_j);
        double l_ii    = computeDistance(x_j, x_k);
        double l_theta = computeAngle(x_i, x_j, x_k);

        if (!refining_angle) theta_lim_aux = theta_lim;

        bool inserted = false;
        if (l_theta > theta_lim_aux) {
            if (l_i > l_ii) {
                insertPoint(prev_seg, arcs[j - 1], overlap);
                if (j != 1) --j;
            } else {
                insertPoint(prev_seg, arcs[j], overlap);
            }
            ++insertion_count;
            inserted       = true;
            refining_angle = true;
        } else if (l_i > l_lim) {
            insertPoint(prev_seg, arcs[j - 1], overlap);
            ++insertion_count;
            if (j != 1) --j;
            inserted       = true;
            refining_angle = false;
        } else if (l_ii > l_lim) {
            insertPoint(prev_seg, arcs[j], overlap);
            ++insertion_count;
            inserted       = true;
            refining_angle = false;
        } else {
            new_seg.push_back(x_i);
            if (j == arcs.size() - 1) {
                new_seg.push_back(x_j);
                new_seg.push_back(x_k);
            }
            ++j;
            insertion_count = 0;
        }

        if (overlap) { theta_lim_aux *= 1.5; overlap = false; }

        if (insertion_count >= max_insertions) {
            std::cerr << "Warning: Maximum number of insertions reached. "
                         "Stopping refinement at this segment.\n";
            break;
        }

        if (inserted) arcs = buildInterpolants(prev_seg);
    }

    outputData.push_back(new_seg);
    return new_seg;
}

std::vector<point> manifold::newSegment(std::vector<point> &prev_seg, double Phi,
                                         int nSeg, double l_lim, double theta_lim) {
    std::vector<point> new_seg;
    size_t j = 1;
    theta_lim *= M_PI / 180.0;
    double theta_lim_aux  = theta_lim;
    int    insertion_count = 0;
    bool   refining_angle  = false;
    bool   overlap         = false;

    while (j < prev_seg.size() - 1) {
        point x_i = apply_map(prev_seg[j - 1].R, prev_seg[j - 1].Z, Phi, nSeg);
        point x_j = apply_map(prev_seg[j].R,     prev_seg[j].Z,     Phi, nSeg);
        point x_k = apply_map(prev_seg[j + 1].R, prev_seg[j + 1].Z, Phi, nSeg);

        double l_i     = computeDistance(x_i, x_j);
        double l_ii    = computeDistance(x_j, x_k);
        double l_theta = computeAngle(x_i, x_j, x_k);

        if (!refining_angle) theta_lim_aux = theta_lim;

        if (l_theta > theta_lim_aux) {
            if (l_i > l_ii) {
                insertPoint(prev_seg, j, overlap);
                if (j != 1) --j;
            } else {
                insertPoint(prev_seg, j + 1, overlap);
            }
            ++insertion_count;
            refining_angle = true;
        } else if (l_i > l_lim) {
            insertPoint(prev_seg, j, overlap);
            ++insertion_count;
            if (j != 1) --j;
            refining_angle = false;
        } else if (l_ii > l_lim) {
            insertPoint(prev_seg, j + 1, overlap);
            ++insertion_count;
            refining_angle = false;
        } else {
            new_seg.push_back(x_i);
            if (j == prev_seg.size() - 2) {
                new_seg.push_back(x_j);
                new_seg.push_back(x_k);
            }
            ++j;
            insertion_count = 0;
        }

        if (overlap) { theta_lim_aux *= 1.5; overlap = false; }

        if (insertion_count >= max_insertions) {
            std::cerr << "Warning: Maximum number of insertions reached. "
                         "Stopping refinement at this segment.\n";
            break;
        }
    }

    outputData.push_back(new_seg);
    return new_seg;
}

bool manifold::save(const std::string &path) const {
    auto ends_with = [&](const std::string &ext) {
        return path.size() >= ext.size() &&
               path.substr(path.size() - ext.size()) == ext;
    };

    bool is_csv = ends_with(".csv");
    bool is_txt = ends_with(".dat") || ends_with(".txt");

    if (!is_csv && !is_txt) {
        std::cerr << "Error: unsupported extension in \"" << path
                  << "\". Use .dat, .txt, or .csv." << std::endl;
        return false;
    }

    std::ofstream f(path);
    if (!f.is_open()) {
        std::cerr << "Error: could not open output file: " << path << std::endl;
        return false;
    }

    f << std::fixed << std::setprecision(16);
    if (is_csv) {
        f << "seg,R,Z\n";
        for (size_t s = 0; s < outputData.size(); ++s)
            for (const auto &pt : outputData[s])
                f << s << "," << pt.R << "," << pt.Z << "\n";
    } else {
        f << "#seg R" << std::string(18, ' ') << "Z\n";
        for (size_t s = 0; s < outputData.size(); ++s)
            for (const auto &pt : outputData[s])
                f << s << " " << pt.R << " " << pt.Z << "\n";
    }

    return true;
}

void manifold::progressBar(int j, int nSeg) {
    const int barWidth = 50;
    float     progress = static_cast<float>(j) / nSeg;

    std::cout << "\nComputing primary segment " << j + 1 << " of " << nSeg << "...\n[";
    int pos = static_cast<int>(barWidth * progress);
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos)       std::cout << "=";
        else if (i == pos) std::cout << ">";
        else               std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r\n";
    std::cout.flush();
}

std::vector<interpolantArc> manifold::buildInterpolants(const std::vector<point> &segment) {
    std::vector<interpolantArc> arcs;
    for (size_t i = 1; i < segment.size() - 2; ++i) {
        const point &p1 = segment[i + 1];
        const point &p2 = segment[i + 2];

        if ((p1.R == 0 && p1.Z == 0) || (p2.R == 0 && p2.Z == 0)) {
            std::cerr << "Skipping arc between invalid points p1 or p2.\n";
            continue;
        }

        double theta   = computeAngle(segment[i - 1], segment[i],     segment[i + 1]);
        double theta_1 = computeAngle(segment[i],     segment[i + 1], segment[i + 2]);

        double m   = std::tan(theta);
        double m_1 = std::tan(theta_1);
        double a   = (std::sqrt(1 + m * m) - 1) / m;
        double b   = -(std::sqrt(1 + m_1 * m_1) - 1) / m_1;

        arcs.push_back({segment[i], segment[i + 1], a, b,
                        static_cast<int>(i), static_cast<int>(i + 1)});
    }
    return arcs;
}

point interpolantArc::evalNewPoint(double t) const {
    if (x0 == point{0, 0} && x1 == point{0, 0}) {
        std::cerr << "Error: Interpolating between two zero points!\n";
        std::exit(EXIT_FAILURE);
    }

    point  base   = x0 * (1.0 - t) + x1 * t;
    point  l      = x1 - x0;
    point  n      = {-l.Z, l.R};                                    // perpendicular
    double hshape = a * t * (1 - t) * (1 - t) - b * t * t * (1 - t); // shape function
    point  result = base + n * hshape;

    if (std::isnan(result.R) || std::isnan(result.Z))
        return base;
    return result;
}
