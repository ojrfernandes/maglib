#ifndef MANIFOLD_H
#define MANIFOLD_H
#define MANIFOLD_V 260531 // version (yy.mm.dd)

#include <maglit.h>
#include <string>
#include <vector>

struct point {
    double R;
    double Z;

    point operator+(const point &o) const { return {R + o.R, Z + o.Z}; }
    point operator-(const point &o) const { return {R - o.R, Z - o.Z}; }
    point operator*(double s)        const { return {R * s,   Z * s};   }
    friend point operator*(double s, const point &p) { return p * s; }
    bool operator==(const point &o)  const { return R == o.R && Z == o.Z; }
};

struct interpolantArc {
    point  x0, x1;
    double a, b;
    int    i0, i1; // indices of the endpoints in the segment vector

    point evalNewPoint(double t) const;
};

class manifold {
  public:
    manifold(maglit &tracer, double phi, int stability);

    // Iteratively find the closest 1-period fixed point from the initial guess
    bool find_xPoint(double rGuess, double zGuess);
    // Compute the primary segment (n_intervals+1 points)
    std::vector<point> primarySegment(size_t n_intervals);
    // Compute a refined new segment from a previous segment (interpolant method)
    std::vector<point> newSegment(std::vector<point> &prev_seg,
                                  double l_lim, double theta_lim);
    // Compute a refined new segment by applying the map nSeg times (exact-map method)
    std::vector<point> newSegment(std::vector<point> &prev_seg, int nSeg,
                                  double l_lim, double theta_lim);
    // Compute all segments in one call: primarySegment + (n_segments-1) newSegment calls
    void run(size_t n_intervals, int n_segments, int method,
             double l_lim, double theta_lim);
    // Print a progress bar
    void progressBar(int j, int nSeg);
    // Enable verbose diagnostic output
    void setVerbose();
    // Configure numerical parameters
    void configure(double epsilon, double h, double tol, int max_iter,
                   double precision_limit, int max_insertions);
    // Save all accumulated segments to file. Format inferred from extension:
    //   .dat / .txt  — space-separated text with header (columns: seg, R, Z)
    //   .csv         — comma-separated text with header (columns: seg, R, Z)
    bool save(const std::string &path) const;

    std::vector<interpolantArc> buildInterpolants(const std::vector<point> &segment);

    point xPoint; // X-point coordinates (set by find_xPoint)
    std::vector<std::vector<point>> outputData; // segments accumulated by primarySegment/newSegment

  private:
    void eval_jacobian(double R, double Z, double Phi, double h, double jacobian[2][2]);
    void insertPoint(std::vector<point> &segment, size_t index, bool &overlap);
    void insertPoint(std::vector<point> &segment, interpolantArc &arc, bool &overlap);
    static double computeDistance(point a, point b);
    static double computeAngle(point a, point b, point c);
    point apply_map(double R, double Z, double Phi, int nTurns);
    point pivot();

    int    stability;
    double phi;

    double epsilon         = 1e-6;
    bool   verbose         = false;
    int    s_factor        = 1;
    double h               = 1e-8;
    double tol             = 1e-14;
    int    max_iter        = 50;
    double precision_limit = 1e-14;
    int    max_insertions  = 100;

    maglit &tracer;
};

#endif // MANIFOLD_H
