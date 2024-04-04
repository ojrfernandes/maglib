#ifndef AUXFIELDS_H
#define AUXFIELDS_H

#include <fusion_io.h>

class auxfields {
  public:
    // opens source and enters additional parameters
    // source_type: FIO_M3DC1_SOURCE, FIO_GEQDSK_SOURCE, FIO_GPEC_SOURCE
    auxfields(const char *source_path, int source_type);
    // eval poloidal flux 'psi' on point x:(R, phi, Z)
    void psi_eval(double &R, double &Phi, double &Z, double *psi);
    // eval normalized poloidal flux 'psin' on point x:(R, phi, Z)
    void psin_eval(double &R, double &Phi, double &Z, double *psin);

  private:
    fio_source     *src;
    fio_field      *psin_field;
    fio_field      *psi_field;
    fio_option_list opt;
};

#endif