#ifndef FOOTPRINT_H
#define FOOTPRINT_H
// Last modified: 26.06.26

#include <atomic>
#include <iomanip>
#include <limits>
#include <maglit.h>
#include <omp.h>
#include <vector>

class footprint {
public:
  // Cache-line-aligned per-thread progress counter (prevents false sharing).
  struct alignas(64) ThreadProgress {
      std::atomic<int> value{0};
  };

  // class constructor
  footprint(const int manifold, const double grid_R1, const double grid_Z1,
            const double grid_R2, const double grid_Z2, const int nRZ,
            const int nPhi, const int max_turns, const int ntor = 1);

  // Magnetic axis (R, Z) [m] used as the origin of the poloidal angle. When not set,
  // run() reads it from tracers[0]->source() and throws std::runtime_error if unavailable.
  void set_magnetic_axis(double R, double Z);

  // Run the grid. Pass one tracer per thread for parallel execution,
  // or a single-element vector for serial execution.
  // When progress is non-null each thread increments its counter after
  // every field line; the display thread in run.cpp renders the bars.
  // When null the built-in single-bar progressBar() is used instead.
  void run(std::vector<maglit *> &tracers,
           std::vector<ThreadProgress> *progress = nullptr);

  // Columns of each outputData row: R0, Z0, phi0, length, psiMin, turn, polTurn
  static constexpr int ncols = 7;

  // Save outputData to file. Format is inferred from the extension:
  //   .dat / .txt  — space-separated with header
  //   .csv         — comma-separated with header
  // Returns false if the extension is unsupported or the file cannot be opened.
  bool save(const std::string &path) const;

  const std::vector<std::vector<double>> &get_output_data() const;

private:
  typedef struct {
    double length;
    double psimin;
    int turn;      // toroidal turn number
    int pol_turn;  // completed poloidal turns around the magnetic axis (net winding)
  } map_scalars; // structure to store the connection length and the minimum psi
                 // value

  // integrate the field line from a given initial condition
  void evolve_line(maglit &tracer, double R0, double Z0, double phi0,
                   double &R1, double &Z1, double &phi1, map_scalars &scalars);
  // calculate the distance between two points
  double connection_length(double R0, double Z0, double phi0, double R1,
                           double Z1, double phi1);
  // display a progress bar
  void progressBar(float progress);

  int manifold;   // 0 = unstable footprint (against-B from target);  1 = stable footprint (following-B from target)
  double grid_R1; // first point R delimiting the target plate mapped surface
  double grid_Z1; // first point Z delimiting the target plate mapped surface
  double grid_R2; // second point R delimiting the target plate mapped surface
  double grid_Z2; // second point Z delimiting the target plate mapped surface
  int nRZ;        // grid dimension along the (R,Z) plane
  int nPhi;       // grid dimension along the phi direction
  int ntor = 1;   // toroidal mode number; phi grid restricted to [0, 2*pi/ntor)
  int max_turns = 1000; // maximum toroidal turns for field line integration
  double R_axis = 0.0;  // magnetic axis R [m]
  double Z_axis = 0.0;  // magnetic axis Z [m]
  bool   axis_set = false;

  std::vector<std::vector<double>> outputData;
};

// Change of the poloidal angle about (R_ax, Z_ax) from (R0, Z0) to (R1, Z1), wrapped into (-pi, pi].
double poloidal_angle_step(double R_ax, double Z_ax, double R0, double Z0, double R1, double Z1);

#endif // FOOTPRINT_H
