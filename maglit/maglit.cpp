#include "maglit.h"

maglit::maglit(FieldSource &source)
    : source_(&source), solver_(SODE_RK56_CK, 2) {
    solver_.set_system(mag_system);
    solver_.configure(1e-8, 1e-12, 1e-15, 0.9);
}

void maglit::inverse_map(bool inverse) {
    inv_factor = inverse ? -1 : 1;
}

bool maglit::calc_mag_field(double *x, double *B) {
    bool ok = source_->eval_B(x[0], x[1], x[2], B);
    if (!ok && warnings)
        std::cerr << "maglit: eval_B failed at R=" << x[0]
                  << " phi=" << x[1] << " Z=" << x[2] << std::endl;
    return ok;
}

void maglit::configure(double dphi_init, double dphi_min, double dphi_max) {
    solver_.configure(dphi_init, dphi_min, dphi_max);
}

int maglit::step(double &R, double &Z, double &phi, double phi_max, int dir) {
    x[0]       = R;
    x[1]       = Z;
    int status = solver_.evolve(x, &phi, phi_max, dir, this);
    R          = x[0];
    Z          = x[1];
    return status;
}

void maglit::reset() {
    solver_.reset();
}

void maglit::set_verb() {
    verb = true;
    solver_.set_verb();
}

void maglit::set_warnings() {
    warnings = true;
}

void maglit::psin_eval(double &R, double &Phi, double &Z, double *psin) {
    if (!source_->eval_psin(R, Phi, Z, *psin) && warnings)
        std::cerr << "maglit: eval_psin failed" << std::endl;
}

void maglit::psi_eval(double &R, double &Phi, double &Z, double *psi) {
    if (!source_->eval_psi(R, Phi, Z, *psi) && warnings)
        std::cerr << "maglit: eval_psi failed" << std::endl;
}

void maglit::set_monitor(const std::string &collider_path) {
    if (!boundary.load_shape(collider_path)) {
        std::cerr << "maglit: error loading collider shape from "
                  << collider_path << std::endl;
        exit(EXIT_FAILURE);
    }
    solver_.set_monitor(&maglit::monitor_boundary);
    solver_.set_aux(this);
}

bool maglit::monitor_boundary(double *x, double t, void *mgl) {
    maglit *self = static_cast<maglit *>(mgl);
    return self->boundary.inside(x[0], x[1]);
}

int mag_system(double *f, double *x, double t, void *mgl) {
    maglit *tracer = static_cast<maglit *>(mgl);
    double  aux_x[3] = {x[0], t, x[1]}; // R, phi, Z
    double  aux_b[3];
    if (!tracer->calc_mag_field(aux_x, aux_b))
        return 1;
    f[0] = (aux_x[0] * aux_b[0] / aux_b[1]) * tracer->inv_factor;
    f[1] = (aux_x[0] * aux_b[2] / aux_b[1]) * tracer->inv_factor;
    return 0;
}
