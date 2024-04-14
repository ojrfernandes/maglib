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
    maglit(const char *source_path, int source_type);
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
    // evaluate the normalizes poloidal flux
    void psin_eval(double &R, double &Phi, double &Z, double *psin);

    // public variables
    bool (*inside)(double R, double Z, double phi, void *aux);
    void *aux; // additional variables, e.g. for region monitor

  private:
    bool            verb;
    fio_source     *src;
    fio_field      *mag_field;
    fio_field      *psin_field;
    fio_option_list opt;
    sode            solver;
    fio_hint        hint;
    double          x[2]; // auxiliary orbit variable
};

// dynamical system x: (R,z); t: phi
double aux_x[3];
double aux_b[3];
int    mag_system(double *f, double *x, double t, void *mgl);
bool   mag_monitor(double *x, double t, void *mgl);

#endif // MAGLIT_H