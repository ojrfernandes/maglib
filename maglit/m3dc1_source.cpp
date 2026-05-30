#include "m3dc1_source.h"
#include <iostream>

M3DC1Source::M3DC1Source(const char *path, int timeslice) {
    int result = fio_open_source(&src_, FIO_M3DC1_SOURCE, path);
    if (result != FIO_SUCCESS) {
        std::cerr << "M3DC1Source: error opening " << path << std::endl;
        src_ = nullptr;
        return;
    }

    src_->get_field_options(&opt_);
    opt_.set_option(FIO_TIMESLICE, timeslice);
    opt_.set_option(FIO_PART, FIO_TOTAL);

    result = src_->get_field(FIO_MAGNETIC_FIELD, &b_field_, &opt_);
    if (result != FIO_SUCCESS) {
        std::cerr << "M3DC1Source: error loading magnetic field" << std::endl;
        b_field_ = nullptr;
        return;
    }

    // Hint accelerates finite-element mesh searches; optional (nullptr falls back to full search).
    result = src_->allocate_search_hint(&hint_);
    if (result != FIO_SUCCESS || hint_ == nullptr)
        hint_ = nullptr;
}

M3DC1Source::~M3DC1Source() {
    if (hint_)       src_->deallocate_search_hint(&hint_);
    if (b_field_)    fio_close_field(&b_field_);
    if (psin_field_) fio_close_field(&psin_field_);
    if (psi_field_)  fio_close_field(&psi_field_);
    if (src_)        fio_close_source(&src_);
}

bool M3DC1Source::eval_B(double R, double phi, double Z, double B[3]) {
    double x[3] = {R, phi, Z};
    return b_field_->eval(x, B, hint_) == FIO_SUCCESS;
}

bool M3DC1Source::eval_psin(double R, double phi, double Z, double &psin) {
    if (!psin_field_) {
        int status = src_->get_field(FIO_POLOIDAL_FLUX_NORM, &psin_field_, &opt_);
        if (status != FIO_SUCCESS) {
            std::cerr << "M3DC1Source: error loading psin field" << std::endl;
            psin_field_ = nullptr;
            return false;
        }
    }
    double x[3] = {R, phi, Z};
    return psin_field_->eval(x, &psin, hint_) == FIO_SUCCESS;
}

bool M3DC1Source::eval_psi(double R, double phi, double Z, double &psi) {
    if (!psi_field_) {
        int status = src_->get_field(FIO_POLOIDAL_FLUX, &psi_field_, &opt_);
        if (status != FIO_SUCCESS) {
            std::cerr << "M3DC1Source: error loading psi field" << std::endl;
            psi_field_ = nullptr;
            return false;
        }
    }
    double x[3] = {R, phi, Z};
    return psi_field_->eval(x, &psi, hint_) == FIO_SUCCESS;
}
