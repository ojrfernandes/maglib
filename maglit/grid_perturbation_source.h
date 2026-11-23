#pragma once

#include "field_source.h"
#include <cstddef>
#include <string>

// Perturbation field on a uniform (phi, R, Z) grid, read from a restore_field() stream binary:
//
//   int32  Nf, Nr, Nz, Nc                    (Nf-1 stored planes, Nc components)
//   double Ra, Rb, za, zb, fa, fb            (m, m, m, m, rad, rad)
//   double ba(Nc, Nf-1, Nr, Nz)              Fortran order; components 0=B_phi, 1=B_R, 2=B_Z [T]
//
// Little-endian. Planes phi_i = fa + i (fb-fa)/(Nf-1), i = 0..Nf-2; the grid is periodic in phi with
// period P = fb - fa, which may be 2 pi / N (one field period). P must divide 2 pi.
// Frame: right-handed (R, phi, Z), phi counter-clockwise from above.
//
// Interpolation: tricubic Catmull-Rom (C1, exact at nodes) on B only.
// shortcut: the stored derivative components (Nc = 12) are unused; switch to Hermite interpolation
// with them if adaptive step control struggles with the C1 field.
//
// The file is memory-mapped read-only, so the per-thread instances share one copy in the page cache.
// eval_B returns false outside [Ra,Rb] x [za,zb]. eval_psi/eval_psin always return false (no flux).
class GridPerturbationSource : public FieldSource {
  public:
    explicit GridPerturbationSource(const std::string &path);
    ~GridPerturbationSource() override;

    GridPerturbationSource(const GridPerturbationSource &)            = delete;
    GridPerturbationSource &operator=(const GridPerturbationSource &) = delete;

    bool eval_B(double R, double phi, double Z, double B[3]) override;
    bool eval_psin(double, double, double, double &) override { return false; }
    bool eval_psi(double, double, double, double &) override { return false; }

    bool   is_valid() const { return data_ != nullptr; }
    double period() const { return P_; } // toroidal period of the stored data [rad]

  private:
    void         *map_      = nullptr;
    size_t        map_size_ = 0;
    const double *data_     = nullptr; // ba(Nc, Np, Nr, Nz), Fortran order

    int    Nc_ = 0, Np_ = 0, Nr_ = 0, Nz_ = 0; // Np_ = Nf - 1 stored planes
    double Ra_ = 0, Rb_ = 0, za_ = 0, zb_ = 0, fa_ = 0, P_ = 0;
    double dR_ = 0, dZ_ = 0, dphi_ = 0;
};
