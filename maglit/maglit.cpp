#include "maglit.h"

// maglit class constructor
maglit::maglit(const char *source_path, int source_type, const int timeslice) : solver(SODE_RK56_CK, 2) {
    // load source from file
    int result = fio_open_source(&src, source_type, source_path);
    if (result != FIO_SUCCESS) {
        std::cerr << "Error opening file: " << source_path << std::endl;
        delete (src);
        return;
    }

    // set options for fields obtained from this source
    src->get_field_options(&opt);
    opt.set_option(FIO_TIMESLICE, timeslice);
    opt.set_option(FIO_PART, FIO_TOTAL);

    // get magnetic field from source
    result = src->get_field(FIO_MAGNETIC_FIELD, &mag_field, &opt);
    if (result != FIO_SUCCESS) {
        std::cerr << "Error opening magnetic field" << std::endl;
        mag_field = 0;
        return;
    }

    // set dynamical system
    solver.set_system(mag_system);

    // configure precision of solution
    solver.configure(1e-8, 1e-12, 1e-15, 0.9);
}

// optional: sets the inverse map of the dynamical system
void maglit::inverse_map(bool inverse) {
    switch (inverse) {
    case true:
        inv_factor = -1;
        break;
    case false:
        inv_factor = 1;
        break;
    }
}

// optional: set inside region of interest
void maglit::set_inside(void *aux, bool (*inside)(double R, double Z, double phi, void *aux)) {
    this->aux = aux;
    this->inside = inside;
    solver.set_monitor(mag_monitor);
}

// evaluate magnetic field B at point x
bool maglit::calc_mag_field(double *x, double *B) {
    int result = mag_field->eval(x, B, hint);
    if (result != FIO_SUCCESS) {
        if (this->warnings) {
            std::cerr << "Fio mag field returned " << result << std::endl;
        }
        return false;
    } else
        return true;
}

// configure solvers step control parameter
void maglit::configure(double dphi_init, double dphi_min, double dphi_max) {
    solver.configure(dphi_init, dphi_min, dphi_max);
}

// performs step up to phi_max or if inside changes
// dir = -1 (true->false), 0 (dont monitor), 1 (false->true)
int maglit::step(double &R, double &Z, double &phi, double phi_max, int dir) {
    x[0] = R;
    x[1] = Z;
    int status = solver.evolve(x, &phi, phi_max, dir, this);
    R = x[0];
    Z = x[1];
    return status;
}

// reset solver
void maglit::reset() {
    solver.reset();
}

// set verbose mode
void maglit::set_verb() {
    this->verb = true;
    solver.set_verb();
}

// set warning mode
void maglit::set_warnings() {
    this->warnings = true;
}

// allocate memory for hint
void maglit::alloc_hint() {
    int result = src->allocate_search_hint(&hint);
    if (hint == nullptr || result != FIO_SUCCESS) {
        std::cerr << "Failed to allocate memory for hint." << std::endl;
        return;
    }
}

// // clear hint memory
void maglit::clear_hint() {
    src->deallocate_search_hint(&hint);
}

// evaluate normalized poloidal flux
void maglit::psin_eval(double &R, double &Phi, double &Z, double *psin) {
    int status = src->get_field(FIO_POLOIDAL_FLUX_NORM, &psin_field, &opt);
    if (status != FIO_SUCCESS) {
        std::cerr << "Error evaluating normalized poloidal flux" << std::endl;
        psin_field = 0;
        return;
    }
    double x[3] = {R, Phi, Z};
    bool result = psin_field->eval(x, psin, hint);
    if (result != FIO_SUCCESS && this->warnings) {
        std::cerr << "Fio psin field returned " << result << std::endl;
    }
}

// evaluate poloidal flux
void maglit::psi_eval(double &R, double &Phi, double &Z, double *psi) {
    int status = src->get_field(FIO_POLOIDAL_FLUX, &psi_field, &opt);
    if (status != FIO_SUCCESS) {
        std::cerr << "Error initializing poloidal flux" << std::endl;
        psi_field = 0;
        return;
    }
    double x[3] = {R, Phi, Z};
    bool result = psi_field->eval(x, psi, hint);
    if (result != FIO_SUCCESS && this->warnings) {
        std::cerr << "Fio psi field returned " << result << std::endl;
    }
}

// map of dynamical system x: (R,z); t: phi
int mag_system(double *f, double *x, double t, void *mgl) {
    maglit *tracer = (maglit *)mgl;
    double aux_x[3];
    double aux_b[3];
    aux_x[0] = x[0];                                     // R
    aux_x[1] = t;                                        // phi
    aux_x[2] = x[1];                                     // z
    bool success = tracer->calc_mag_field(aux_x, aux_b); // (B_R, B_phi, B_z)
    if (success) {
        f[0] = (aux_x[0] * aux_b[0] / aux_b[1]) * tracer->inv_factor;
        f[1] = (aux_x[0] * aux_b[2] / aux_b[1]) * tracer->inv_factor;
        return 0;
    } else {
        return 1;
    }
}

// monitor for inside region of interest
bool mag_monitor(double *x, double t, void *mgl) {
    maglit *tracer = (maglit *)mgl;
    // x: (R,z); t: phi
    return tracer->inside(x[0], x[1], t, tracer->aux);
}