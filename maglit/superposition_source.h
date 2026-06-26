#pragma once
// Last modified: 26.06.25

#include "field_source.h"
#include "m3dc1_source.h"
#include <memory>
#include <string>
#include <vector>

// FieldSource that linearly superposes N M3DC1 perturbation components
// onto a shared axisymmetric equilibrium:
//
//   B(R, φ, Z) = B_eq(R, φ, Z) + Σ_i A_i · [B_i(R, φ − δ_i, Z) − B_eq(R, φ, Z)]
//
// The phase shift δ_i is exact: rotating a coil by δ in the toroidal direction
// and evaluating at (R, φ, Z) equals evaluating the unshifted field at (R, φ−δ, Z)
// because the cylindrical basis vectors transform as T_δ·R̂(φ−δ) = R̂(φ) and
// T_δ·φ̂(φ−δ) = φ̂(φ).
//
// The equilibrium is loaded at timeslice = -1 from the first component's file.
// All components must share the same equilibrium (typical for linear M3DC1 runs
// from a common base equilibrium).
//
// psin and psi are delegated to the equilibrium source — flux surfaces are
// defined by the equilibrium and are not materially perturbed by small RMP fields.
//
// Non-copyable: owns M3DC1Source instances with Fusion-IO file handles.
// One instance per thread is required for multi-threaded use.
class SuperpositionSource : public FieldSource {
  public:
    SuperpositionSource()  = default;
    ~SuperpositionSource() = default;

    SuperpositionSource(const SuperpositionSource &)            = delete;
    SuperpositionSource &operator=(const SuperpositionSource &) = delete;

    // Add a perturbation component. Must be called at least once before use.
    // The first call also opens the equilibrium (timeslice = -1) from `path`.
    // path:        M3DC1 HDF5 file for this component
    // timeslice:   0 = vacuum, 1 = full single-fluid response (not -1: that is the equilibrium)
    // phase_shift: δ_i in radians; applied as φ → φ − δ_i at evaluation
    // amplitude:   A_i, linear scale factor (can be negative for anti-phase contributions)
    void add_component(const std::string &path, int timeslice,
                       double phase_shift, double amplitude);

    bool eval_B(double R, double phi, double Z, double B[3]) override;
    bool eval_psin(double R, double phi, double Z, double &psin) override;
    bool eval_psi(double R, double phi, double Z, double &psi) override;

    // Returns false if no components have been added or any source failed to open.
    bool is_valid() const;

    int num_components() const { return static_cast<int>(components_.size()); }

  private:
    struct Component {
        std::unique_ptr<M3DC1Source> pert_source;
        double phase_shift; // radians
        double amplitude;
    };

    std::unique_ptr<M3DC1Source> eq_source_;
    std::vector<Component>       components_;
};
