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
            const int nPhi, const int max_turns);

  // Run the grid. Pass one tracer per thread for parallel execution,
  // or a single-element vector for serial execution.
  // When progress is non-null each thread increments its counter after
  // every field line; the display thread in run.cpp renders the bars.
  // When null the built-in single-bar progressBar() is used instead.
  void run(std::vector<maglit *> &tracers,
           std::vector<ThreadProgress> *progress = nullptr);

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
    int turn;    // toroidal turn number
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
  int max_turns = 1000; // maximum toroidal turns for field line integration

  std::vector<std::vector<double>> outputData;
};

#endif // FOOTPRINT_H
