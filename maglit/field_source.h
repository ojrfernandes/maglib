#pragma once

// Abstract interface for magnetic field data sources.
//
// Concrete implementations provide B, psi_n, and psi evaluated at cylindrical
// coordinates (R, phi, Z). Each implementation is responsible for managing its
// own internal state (file handles, mesh search hints, lazy-loaded fields, etc.).
//
// One FieldSource instance per thread is required for multi-threaded use since
// implementations may hold mutable internal state (search hints, cached fields).
class FieldSource {
  public:
    virtual ~FieldSource() = default;

    // Evaluate magnetic field B = (B_R, B_phi, B_Z) at (R, phi, Z).
    // Returns true on success, false if the point is outside the domain.
    virtual bool eval_B(double R, double phi, double Z, double B[3]) = 0;

    // Evaluate normalized poloidal flux at (R, phi, Z).
    // Returns true on success.
    virtual bool eval_psin(double R, double phi, double Z, double &psin) = 0;

    // Evaluate poloidal flux at (R, phi, Z).
    // Returns true on success.
    virtual bool eval_psi(double R, double phi, double Z, double &psi) = 0;
};
