#ifndef MAGLIT_H
#define MAGLIT_H
#define MAGLIT_V 241119 // version (yy.mm.dd)

#include "fusion_io.h"
#include "sode.h"
#include <functional>
#include <iostream>

class maglit {
  public:
    // maglit class constructor
    // opens source and enters additional parameters
    // source_type: FIO_M3DC1_SOURCE, FIO_GEQDSK_SOURCE, FIO_GPEC_SOURCE
    maglit(const char *source_path, int source_type, int timeslice);

    // class destructor
    ~maglit();

    // set the inverse map of the dynamical system
    void inverse_map(bool inverse);
    // evaluate the magnetic field in cylindrical coordinates
    bool calc_mag_field(double *x, double *B); // x:(R, phi, Z)
    // configure solver's step control parameter
    void configure(double dphi_init, double dphi_min, double dphi_max);
    // performs step up to phi_max or if inside changes
    // dir = -1 (true->false), 0 (dont monitor), 1 (false->true)
    int step(double &R, double &Z, double &phi, double phi_max, int dir);
    // reset integrator (run before start new field line)
    void reset();
    // run to print detailed messages
    void set_verb();
    // run to print warning messages
    void set_warnings();
    // allocate search hint
    void alloc_hint();
    // clear hint allocated memory
    void clear_hint();
    // evaluate the normalized poloidal flux
    void psin_eval(double &R, double &Phi, double &Z, double *psin);
    // evaluate the poloidal flux
    void psi_eval(double &R, double &Phi, double &Z, double *psi);
    // monitor for inside region of interest
    static bool mag_monitor(double *x, double t, void *mgl);

    // template function
    // defines the inside region of interest
    template <typename auxClass>
    void set_inside(auxClass &aux, bool (auxClass::*inside)(double, double, double, void *)) {
        this->aux = static_cast<void *>(&aux);
        this->inside = [inside](double R, double Z, double phi, void *aux) {
            return (static_cast<auxClass *>(aux)->*inside)(R, Z, phi, aux);
        };
        solver.set_monitor(mag_monitor);
    }

    // public variables
    std::function<bool(double, double, double, void *)> inside;         // generic function for "inside" logic
    void                                               *aux;            // void pointer to auxiliary class
    int                                                 inv_factor = 1; // factor for inverse map

  private:
    bool            verb = false;     // verbose mode
    bool            warnings = false; // warning mode
    fio_source     *src;              // fusion-io data source
    fio_field      *mag_field;        // fusion-io magnetic field
    fio_field      *psin_field;       // fusion-io normalized poloidal flux
    fio_field      *psi_field;        // fusion-io poloidal flux
    fio_option_list opt;              // fusion-io field options
    sode            solver;           // solver for the dynamical system
    fio_hint        hint;             // fusion-io hint for finite element search
    double          x[2];             // auxiliary orbit variable
};

// map of dynamical system x: (R,z); t: phi
int mag_system(double *f, double *x, double t, void *mgl);

#endif // MAGLIT_H