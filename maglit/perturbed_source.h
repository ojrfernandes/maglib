#pragma once

#include "equilibrium_source.h"
#include "field_source.h"
#include "grid_perturbation_source.h"
#include <memory>
#include <stdexcept>
#include <string>

// Equilibrium plus a scaled, toroidally shifted perturbation:
//
//   B(R, phi, Z) = B_eq(R, phi, Z) + A * B_pert(R, phi - delta, Z)
//
// psi and psi_n come from the equilibrium. Owns both sources; delta is in radians.
// Fails wherever either source fails (e.g. outside the perturbation grid).
class PerturbedSource : public FieldSource {
  public:
    PerturbedSource(std::unique_ptr<FieldSource> eq, std::unique_ptr<FieldSource> pert,
                    double amplitude, double phase_shift)
        : eq_(std::move(eq)), pert_(std::move(pert)), amplitude_(amplitude), phase_shift_(phase_shift) {}

    bool eval_B(double R, double phi, double Z, double B[3]) override {
        double Bp[3];
        if (!eq_->eval_B(R, phi, Z, B) || !pert_->eval_B(R, phi - phase_shift_, Z, Bp)) return false;
        for (int c = 0; c < 3; ++c)
            B[c] += amplitude_ * Bp[c];
        return true;
    }
    bool eval_psin(double R, double phi, double Z, double &psin) override {
        return eq_->eval_psin(R, phi, Z, psin);
    }
    bool eval_psi(double R, double phi, double Z, double &psi) override {
        return eq_->eval_psi(R, phi, Z, psi);
    }

  private:
    std::unique_ptr<FieldSource> eq_, pert_;
    double                       amplitude_, phase_shift_;
};

// EquilibriumSource from `equilibrium_path`, plus amplitude * GridPerturbationSource(phi - phase_shift)
// when `perturbation_path` is non-empty. phase_shift in radians. Throws std::runtime_error if a file
// fails to load.
inline std::unique_ptr<FieldSource> make_equilibrium_source(const std::string &equilibrium_path,
                                                            const std::string &perturbation_path,
                                                            double amplitude, double phase_shift) {
    auto eq = std::make_unique<EquilibriumSource>(equilibrium_path);
    if (!eq->is_valid())
        throw std::runtime_error("failed to load equilibrium: " + equilibrium_path);
    if (perturbation_path.empty())
        return eq;
    auto pert = std::make_unique<GridPerturbationSource>(perturbation_path);
    if (!pert->is_valid())
        throw std::runtime_error("failed to load perturbation: " + perturbation_path);
    return std::make_unique<PerturbedSource>(std::move(eq), std::move(pert), amplitude, phase_shift);
}
