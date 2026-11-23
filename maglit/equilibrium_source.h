#pragma once

#include "field_source.h"
#include <string>
#include <vector>

// Axisymmetric equilibrium read from a 'magnetics' XML file or a G-EQDSK file.
// The format is detected from content: a first non-blank '<' means XML.
//
// Field (right-handed cylindrical (R, phi, Z), phi counter-clockwise from above):
//   B_R   = -sigma * dpsi/dZ / (2 pi R)
//   B_Z   =  sigma * dpsi/dR / (2 pi R)
//   B_phi =  F / R,   F = F(psi_N) inside the separatrix, F_vac elsewhere
// psi is the total poloidal flux [Wb] in the file's own sign convention; sigma = +-1.
//
// XML:     psi = PFM [T m^2] on (zj, Ri), sigma = +1, F = mu0 Jpol(psi_N) / (2 pi),
//          psi_N boundary = flux of the first genuine saddle in singularities/xp.
// G-EQDSK: psi = 2 pi psirz (psirz in Wb/rad), sigma = sign(Ip) * sign(simag - sibry),
//          F = fpol(psi_N), psi_N boundary = sibry.
//
// psi is interpolated by a natural bicubic spline (C2 psi, C1 B), evaluated cell-wise as a
// bicubic Hermite patch from nodal psi, psi_R, psi_Z, psi_RZ.
//
// Immutable after construction: one instance per thread is cheap (the file is re-read).
class EquilibriumSource : public FieldSource {
  public:
    explicit EquilibriumSource(const std::string &path);

    EquilibriumSource(const EquilibriumSource &)            = delete;
    EquilibriumSource &operator=(const EquilibriumSource &) = delete;

    bool eval_B(double R, double phi, double Z, double B[3]) override;
    bool eval_psin(double R, double phi, double Z, double &psin) override;
    bool eval_psi(double R, double phi, double Z, double &psi) override;

    // Returns false if the file could not be read or is inconsistent.
    bool is_valid() const { return valid_; }

  private:
    bool load_xml(const std::string &text);
    bool load_geqdsk(const std::string &text);
    void build_derivatives();
    bool eval_flux(double R, double Z, double &psi, double &psi_R, double &psi_Z) const;

    // Uniform grid: R_i = R0_ + i dR_, Z_k = Z0_ + k dZ_ [m]
    int    nR_ = 0, nZ_ = 0;
    double R0_ = 0, dR_ = 0, Z0_ = 0, dZ_ = 0;

    // Nodal psi and derivatives, shape=(nZ,nR), row-major (index k*nR + i), units Wb, Wb/m, Wb/m^2
    std::vector<double> psi_, psiR_, psiZ_, psiRZ_;

    double psi_axis_ = 0, psi_bdry_ = 0; // Wb
    double sigma_    = 1;                // sign of the poloidal field relative to grad psi x grad phi

    // F(psi_N) [T m]: abscissae, values, natural-spline second derivatives
    std::vector<double> fx_, fy_, fy2_;
    double              F_vac_ = 0;

    // Z range of the confined region (between the X-points); excludes the private-flux region
    double z_lo_ = 0, z_hi_ = 0;

    bool valid_ = false;
};
