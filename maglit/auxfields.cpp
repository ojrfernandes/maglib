#include "auxfields.h"

// auxfields(char *source_path, int source_type)
auxfields::auxfields(char *source_path, int source_type) {
    // load source from file
    int result = fio_open_source(&src, source_type, source_path);
    if (result != FIO_SUCCESS) {
        printf("Error opening file: %s\n", source_path);
        delete (src);
        return;
    }
    // set options for fields obtained from this source
    src->get_field_options(&opt);
    opt.set_option(FIO_TIMESLICE, 0);
    opt.set_option(FIO_PART, FIO_TOTAL);
}

// int psi_eval(double *x, double *psi)
void auxfields::psi_eval(double &R, double &Phi, double &Z, double *psi) {
    int status = src->get_field(FIO_POLOIDAL_FLUX, &psi_field, &opt);
    if (status != FIO_SUCCESS) {
        printf("Error evaluating poloidal flux\n");
        psi_field = 0;
        return;
    }
    double x[3] = {R, Phi, Z};
    psi_field->eval(x, psi);
    // printf("psi = %f\n", *psi);
}

// int psin_eval(double *x, double *psin)
void auxfields::psin_eval(double &R, double &Phi, double &Z, double *psin) {
    int status = src->get_field(FIO_POLOIDAL_FLUX_NORM, &psin_field, &opt);
    if (status != FIO_SUCCESS) {
        printf("Error evaluating normalized poloidal flux\n");
        psin_field = 0;
        return;
    }
    double x[3] = {R, Phi, Z};
    psin_field->eval(x, psin);
    // printf("psin = %f\n", *psin);
}