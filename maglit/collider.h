#ifndef COLLIDER
#define COLLIDER
#define COLLIDER_V 251103 // version (yy.mm.dd)

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

class collider {
  public:
    // class default constructor
    collider() = default;

    // load (R,Z) vessel vertices from file
    bool load_shape(const std::string &path);

    // check if a point (R,Z) is inside the vessel
    bool inside(double R, double Z) const;

    // check if the shape has been loaded
    bool is_loaded() const { return loaded; }

    // get vessel shape as vector of (R,Z) pairs
    const std::vector<std::pair<double, double>> &get_vertices() const { return vertices; }

  private:
    std::vector<std::pair<double, double>> vertices;       // vessel shape vertices (R,Z)
    bool                                   loaded = false; // flag indicating if shape is loaded

    // normalize angle to the range (-pi, pi]
    static double normalize_angle(double angle);

    // check if point (R,Z) lies on any edge of the polygon
    bool point_on_edge(double R, double Z, double tol = 1e-10) const;
};

#endif