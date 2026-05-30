#ifndef MAGLIT_H
#define MAGLIT_H
#define MAGLIT_V 260530

#include "collider.h"
#include "field_source.h"
#include "sode.h"
#include <iostream>
#include <string>

class maglit {
  public:
    // source must outlive this maglit instance.
    explicit maglit(FieldSource &source);

    ~maglit() = default;

    // set the inverse map of the dynamical system
    void inverse_map(bool inverse);
    // evaluate the magnetic field in cylindrical coordinates
    bool calc_mag_field(double *x, double *B); // x: (R, phi, Z)
    // configure solver step control parameters
    void configure(double dphi_init, double dphi_min, double dphi_max);
    // perform one step up to phi_max, or until the monitor changes
    // dir = -1 (true→false), 0 (don't monitor), 1 (false→true)
    int step(double &R, double &Z, double &phi, double phi_max, int dir);
    // reset integrator (call before starting a new field line)
    void reset();
    // enable verbose output
    void set_verb();
    // enable warning output
    void set_warnings();
    // evaluate the normalized poloidal flux
    void psin_eval(double &R, double &Phi, double &Z, double *psin);
    // evaluate the poloidal flux
    void psi_eval(double &R, double &Phi, double &Z, double *psi);
    // load vessel shape and activate boundary monitor
    void set_monitor(const std::string &path);

    int      inv_factor = 1;
    collider boundary;

  private:
    bool         verb     = false;
    bool         warnings = false;
    FieldSource *source_;  // non-owning
    sode         solver_;
    double       x[2];

    static bool monitor_boundary(double *x, double t, void *mgl);
};

// Right-hand side of the field-line ODE: x = (R, Z), t = phi
int mag_system(double *f, double *x, double t, void *mgl);

#endif // MAGLIT_H
