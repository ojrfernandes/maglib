#ifndef MAGLIT_H
#define MAGLIT_H
#define MAGLIT_V 230604

#include "../sode/sode.h"
#include <fusion_io.h>
#include <iostream>

class maglit {
  public:
    // opens source and enters additional parameters
    // source_type: FIO_M3DC1_SOURCE, FIO_GEQDSK_SOURCE, FIO_GPEC_SOURCE
    maglit(const char *source_path, int source_type, int timeslice);
    // optional: sets the inverse map of the dynamical system
    void inverse_map(bool inverse);
    // optional: defines the inside region of interest
    void set_inside(void *aux, bool (*inside)(double R, double Z, double phi, void *aux));
    // evaluates the magnetic field in cylindrical coordinates
    bool calc_mag_field(double *x, double *B); // x:(R, phi, Z)
    // configure numerical solver
    void configure(double dphi_init, double dphi_min, double dphi_max);
    // performs step up to phi_max or if inside changes
    // dir = -1 (true->false), 0 (dont monitor), 1 (false->true)
    int step(double &R, double &Z, double &phi, double phi_max, int dir);
    // run before start new line
    void reset();
    // run to print detailed messages
    void set_verb();
    // allocate search hint
    void alloc_hint();
    // clear hint allocated memory
    void clear_hint();
    // evaluate the normalized poloidal flux
    void psin_eval(double &R, double &Phi, double &Z, double *psin);
    // evaluate the poloidal flux
    void psi_eval(double &R, double &Phi, double &Z, double *psi);

    // public variables
    bool (*inside)(double R, double Z, double phi, void *aux);
    void *aux; // additional variables, e.g. for region monitor

  private:
    bool verb;
    fio_source *src;
    fio_field *mag_field;
    fio_field *psin_field;
    fio_field *psi_field;
    fio_option_list opt;
    sode solver;
    fio_hint hint;
    double x[2]; // auxiliary orbit variable
};

// dynamical system x: (R,z); t: phi
extern double aux_x[3];
extern double aux_b[3];
// map of dynamical system x: (R,z); t: phi
int mag_system(double *f, double *x, double t, void *mgl);
// inverse map of dynamical system x: (R,z); t: phi
int inverse_mag_system(double *f, double *x, double t, void *mgl);
// monitor for inside region of interest
bool mag_monitor(double *x, double t, void *mgl);

#endif // MAGLIT_H