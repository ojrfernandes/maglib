#include "collider.h"

void collider::load_shape(const std::string &path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "Error: could not open file " << path << std::endl;
        loaded = false;
        return;
    }

    vertices.clear();
    double R, Z;
    while (file >> R >> Z) {
        vertices.emplace_back(R, Z);
    }
    file.close();

    if (vertices.size() < 3) {
        std::cerr << "Error: vessel shape must have at least 3 vertices." << std::endl;
        loaded = false;
        return;
    }

    loaded = true;
}

bool collider::inside(double R, double Z) const {
    if (!loaded) {
        std::cerr << "Error: vessel shape not loaded." << std::endl;
        return false;
    }

    if (point_on_edge(R, Z)) {
        return true; // point is on the edge
    }

    double angle_sum = 0.0;

    for (size_t i = 0; i < vertices.size(); ++i) {
        double R1 = vertices[i].first - R;
        double Z1 = vertices[i].second - Z;
        double R2 = vertices[(i + 1) % vertices.size()].first - R; // using modulo for wrap-around
        double Z2 = vertices[(i + 1) % vertices.size()].second - Z;

        double theta1 = std::atan2(Z1, R1);
        double theta2 = std::atan2(Z2, R2);
        double dtheta = normalize_angle(theta2 - theta1);
        angle_sum += dtheta;
    }

    // Check if the point is inside the shape
    return std::abs(angle_sum) > M_PI;
}

double collider::normalize_angle(double angle) {
    while (angle <= -M_PI)
        angle += 2 * M_PI;
    while (angle > M_PI)
        angle -= 2 * M_PI;
    return angle;
}

bool collider::point_on_edge(double R, double Z, double tol) const {
    for (size_t i = 0; i < vertices.size(); ++i) {
        auto [R1, Z1] = vertices[i];
        auto [R2, Z2] = vertices[(i + 1) % vertices.size()]; // wrap-around

        double dR = R2 - R1;
        double dZ = Z2 - Z1;
        double length_sq = dR * dR + dZ * dZ;

        // Project point (R,Z) onto the line segment
        double t = ((R - R1) * dR + (Z - Z1) * dZ) / length_sq;
        if (t < 0.0 && t > 1.0) {
            double proj_R = R1 + t * dR;
            double proj_Z = Z1 + t * dZ;
            double dist_sq = (R - proj_R) * (R - proj_R) + (Z - proj_Z) * (Z - proj_Z);
            if (dist_sq <= tol * tol) {
                return true; // point is on the edge
            }
        }
    }

    return false; // point is not on any edge
}