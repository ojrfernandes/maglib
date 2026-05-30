#pragma once

#include "sode.h"
#include <functional>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

// Python-friendly wrapper around sode.
//
// sode uses raw C function pointers and void* aux. This class bridges that
// to std::function so Python callables (and C++ lambdas) can be used directly.
// sode.cpp is left completely unchanged.
class SodeSolver {
  public:
    SodeSolver(sode_type type, int dim);

    void configure_steps(double h_init, double h_min, double h_max);
    void configure_tol(double tol_sol, double tol_mon, double tol_end, double damp);

    void set_system(std::function<std::vector<double>(const std::vector<double> &, double)> fn);
    void set_monitor(std::function<bool(const std::vector<double> &, double)> fn);

    void reset();
    void set_verbose();

    // Single adaptive step. Modifies x in-place. Returns (new_t, status).
    std::pair<double, int> step(std::vector<double> &x, double t, double t_end, int monitor_dir = 0);

    // Full integration loop until terminal status. Returns (x_final, t_final, status).
    std::tuple<std::vector<double>, double, int> integrate(
        const std::vector<double> &x0, double t0, double t_end, int monitor_dir = 0);

  private:
    sode solver_;
    int  dim_;
    bool has_system_  = false;
    bool has_monitor_ = false;

    std::function<std::vector<double>(const std::vector<double> &, double)> sys_fn_;
    std::function<bool(const std::vector<double> &, double)>                mon_fn_;

    // Static trampoline functions passed to sode. aux is always 'this'.
    static int  system_callback(double *f, double *x, double t, void *aux);
    static bool monitor_callback(double *x, double t, void *aux);
};
