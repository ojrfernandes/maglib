#ifndef MANIFOLD_H
#define MANIFOLD_H
#define MANIFOLD_V 241023 // version (yy.mm.dd)

#include "../../maglit.h"
#include <armadillo>
#include <chrono>
#include <iomanip>
#include <thread>
#include <utility>

// Define a point structure
struct point {
    double R;
    double Z;
};

class manifold {
  public:
    // Constructor
    manifold(const char *source_path, const int timeslice, double phi, int stability, double epsilon);
    // Iteratively find the closest 1 period fixed point from the initial guess
    bool find_xPoint(double rGuess, double zGuess);
    // Compute the primary segment
    void primarySegment(std::vector<point> &segment, size_t num_points);
    // Compute a refined new segment from a previous segment
    void newSegment(std::vector<point> &prev_seg, std::vector<point> &new_seg, double Phi, int nSeg, double l_lim, double theta_lim);
    // Print a progress bar
    void progressBar(float progress);
    // Set warning flag
    void setWarnings();

    point xPoint; // x-point

  private:
    // Evaluate the jacobian of the map at a given point
    void eval_jacobian(double R, double Z, double Phi, double h, double jacobian[2][2]);
    // Insert a new point in the vector by linear interpolation
    void insertPoint(std::vector<point> &segment, size_t index);
    // Compute distance between two points
    double computeDistance(double R1, double Z1, double R2, double Z2);
    // Compute angle between two vectors
    double computeAngle(double R0, double Z0, double R1, double Z1, double R2, double Z2);
    // Apply map to a given point returning a Point (R, Z)
    point apply_map(double R, double Z, double Phi, int nTurns);
    // Find the pivot point for the primary segment
    point pivot();

    // User defined parameters
    int stability;  // 0: forward, 1: backward
    double phi;     // toroidal angle
    double epsilon; // distance from the x-point

    // Default parameters
    bool warnings = false;          // verbose flag
    int s_factor = 1;               // sign factor for the manifold stability
    double h = 1e-8;                // step size for numerical differentiation
    double tol = 1e-14;             // tolerance for Newton's method
    int max_iter = 50;              // maximum number of iterations for Newton's method
    double precision_limit = 1e-14; // precision limit for floating point comparisons
    bool refining_angle = false;    // flag to indicate if the angle is being refined
    bool overlap = false;           // flag to indicate floating point precision limit overlap
    int max_insertions = 100;       // Max limit for new point insertions in a segment

    maglit tracer; // magnetic field lines tracer
};

#endif // MANIFOLD_H