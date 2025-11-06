#include "footprint.h"

footprint::footprint(const int manifold, const double grid_R1, const double grid_Z1, const double grid_R2, const double grid_Z2, const int nRZ, const int nPhi) : manifold(manifold), grid_R1(grid_R1), grid_Z1(grid_Z1), grid_R2(grid_R2), grid_Z2(grid_Z2), nRZ(nRZ), nPhi(nPhi) {
    this->outputData.resize(nRZ * nPhi, std::vector<double>(5));
}

void footprint::runGrid(maglit &tracer) {
    map_scalars scalars;

    if (manifold == 0) {
        // evaluating unstable manifold
        // since the field line is traced backwards, the map remains direct
        tracer.inverse_map(false);
    } else {
        // evaluating stable manifold
        // since the field line is traced backwards, the map becomes inverse
        tracer.inverse_map(true);
    }

// loop over the grid
#pragma omp for schedule(dynamic)
    for (int i = 0; i < nPhi; i++) {
        for (int j = 0; j < nRZ; j++) {

            // set initial conditions
            double R_init = grid_R1 + (grid_R2 - grid_R1) * j / (nRZ - 1);
            double Z_init = grid_Z1 + (grid_Z2 - grid_Z1) * j / (nRZ - 1);
            double phi_init = 2 * M_PI * i / nPhi;

            // variables to store final conditions
            double R_final, phi_final, Z_final;

            // evolve lines
            tracer.alloc_hint();
            this->evolve_line(tracer, R_init, Z_init, phi_init, R_final, Z_final, phi_final, scalars);
            tracer.clear_hint();

            // Calculate the index for dataWrite
            int index = i * nRZ + j;

            // Assign the values directly into the preallocated matrix
            outputData[index][0] = R_init;
            outputData[index][1] = Z_init;
            outputData[index][2] = phi_init;
            outputData[index][3] = scalars.length;
            outputData[index][4] = scalars.psimin;

            // Print progress bar only in the first thread (thread 0)
            if (omp_get_thread_num() == 0) {
                float progress = (float)(i * nRZ + j) / (nRZ * nPhi);
                this->progressBar(progress);
                std::this_thread::sleep_for(std::chrono::milliseconds(50));
            }
        }
    }
}

void footprint::evolve_line(maglit &tracer, double R0, double Z0, double phi0, double &R1, double &Z1, double &phi1, map_scalars &scalars) {
    R1 = R0;
    Z1 = Z0;
    phi1 = phi0;
    int     status;
    double  phi_max = 10000 * 2 * M_PI;
    double  arc = 0;
    double  psin0 = 5.0;
    double *psin1 = &psin0;
    scalars.psimin = *psin1;
    tracer.reset();
    do {
        R0 = R1;
        Z0 = Z1;
        phi0 = phi1;
        status = tracer.step(R1, Z1, phi1, phi_max, -1);
        if (status == SODE_CONTINUE_GOOD_STEP) {
            arc += this->connection_length(R0, Z0, phi0, R1, Z1, phi1);
            tracer.psin_eval(R1, phi1, Z1, psin1);
            if (*psin1 < scalars.psimin) {
                scalars.psimin = *psin1;
            }
        }
    } while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP);
    scalars.length = arc;
}

double footprint::connection_length(double R0, double Z0, double phi0, double R1, double Z1, double phi1) {
    double dx = R1 * cos(phi1) - R0 * cos(phi0);
    double dy = R1 * sin(phi1) - R0 * sin(phi0);
    double dz = Z1 - Z0;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

// Print a progress bar
void footprint::progressBar(float progress) {
    int barWidth = 50;

    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}