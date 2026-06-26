#include "superposition_source.h"
#include <iostream>

void SuperpositionSource::add_component(const std::string &path, int timeslice,
                                        double phase_shift, double amplitude) {
    // First component also provides the equilibrium (timeslice = -1).
    if (!eq_source_) {
        eq_source_ = std::make_unique<M3DC1Source>(path.c_str(), -1);
        if (!eq_source_->is_valid()) {
            std::cerr << "SuperpositionSource: failed to open equilibrium from "
                      << path << std::endl;
            eq_source_.reset();
            return;
        }
    }

    auto pert = std::make_unique<M3DC1Source>(path.c_str(), timeslice);
    if (!pert->is_valid()) {
        std::cerr << "SuperpositionSource: failed to open timeslice " << timeslice
                  << " from " << path << std::endl;
        return;
    }

    components_.push_back({std::move(pert), phase_shift, amplitude});
}

bool SuperpositionSource::eval_B(double R, double phi, double Z, double B[3]) {
    double B_eq[3];
    if (!eq_source_->eval_B(R, phi, Z, B_eq))
        return false;

    B[0] = B_eq[0];
    B[1] = B_eq[1];
    B[2] = B_eq[2];

    for (const auto &comp : components_) {
        double B_total[3];
        if (!comp.pert_source->eval_B(R, phi - comp.phase_shift, Z, B_total))
            return false;
        B[0] += comp.amplitude * (B_total[0] - B_eq[0]);
        B[1] += comp.amplitude * (B_total[1] - B_eq[1]);
        B[2] += comp.amplitude * (B_total[2] - B_eq[2]);
    }
    return true;
}

bool SuperpositionSource::eval_psin(double R, double phi, double Z, double &psin) {
    return eq_source_->eval_psin(R, phi, Z, psin);
}

bool SuperpositionSource::eval_psi(double R, double phi, double Z, double &psi) {
    return eq_source_->eval_psi(R, phi, Z, psi);
}

bool SuperpositionSource::is_valid() const {
    if (!eq_source_ || !eq_source_->is_valid() || components_.empty())
        return false;
    for (const auto &comp : components_)
        if (!comp.pert_source || !comp.pert_source->is_valid())
            return false;
    return true;
}
