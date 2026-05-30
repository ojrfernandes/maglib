#include "sode_wrapper.h"

SodeSolver::SodeSolver(sode_type type, int dim) : solver_(type, dim), dim_(dim) {}

void SodeSolver::configure_steps(double h_init, double h_min, double h_max) {
    solver_.configure(h_init, h_min, h_max);
}

void SodeSolver::configure_tol(double tol_sol, double tol_mon, double tol_end, double damp) {
    solver_.configure(tol_sol, tol_mon, tol_end, damp);
}

void SodeSolver::set_system(std::function<std::vector<double>(const std::vector<double> &, double)> fn) {
    sys_fn_     = std::move(fn);
    has_system_ = true;
    solver_.set_system(system_callback);
}

void SodeSolver::set_monitor(std::function<bool(const std::vector<double> &, double)> fn) {
    mon_fn_      = std::move(fn);
    has_monitor_ = true;
    solver_.set_monitor(monitor_callback);
    solver_.set_aux(this);
}

void SodeSolver::reset() { solver_.reset(); }
void SodeSolver::set_verbose() { solver_.set_verb(); }

std::pair<double, int> SodeSolver::step(std::vector<double> &x, double t, double t_end, int monitor_dir) {
    if (!has_system_)
        throw std::runtime_error("System function not set — call set_system() first.");
    int status = solver_.evolve(x.data(), &t, t_end, monitor_dir, this);
    return {t, status};
}

std::tuple<std::vector<double>, double, int> SodeSolver::integrate(
    const std::vector<double> &x0, double t0, double t_end, int monitor_dir) {
    if (!has_system_)
        throw std::runtime_error("System function not set — call set_system() first.");
    std::vector<double> x = x0;
    double              t = t0;
    solver_.reset();
    int status;
    do {
        status = solver_.evolve(x.data(), &t, t_end, monitor_dir, this);
    } while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP);
    return {x, t, status};
}

int SodeSolver::system_callback(double *f, double *x, double t, void *aux) {
    auto *              self = static_cast<SodeSolver *>(aux);
    std::vector<double> xv(x, x + self->dim_);
    auto                result = self->sys_fn_(xv, t);
    for (int i = 0; i < self->dim_; i++)
        f[i] = result[i];
    return 0;
}

bool SodeSolver::monitor_callback(double *x, double t, void *aux) {
    auto *              self = static_cast<SodeSolver *>(aux);
    std::vector<double> xv(x, x + self->dim_);
    return self->mon_fn_(xv, t);
}
