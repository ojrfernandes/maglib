#include "footprint.h"

footprint::footprint(const double &plate, const double &gridMin, const double &gridMax, const int &nGrid, const int &nPhi) {
    this->plate = plate;
    this->gridMin = gridMin;
    this->gridMax = gridMax;
    this->nGrid = nGrid;
    this->nPhi = nPhi;
    this->outputData.resize(nGrid * nPhi, std::vector<double>(9));
    ;
}

void footprint::runGrid(maglit &tracer) {
    map_scalars scalars;
    double R0, Z0;

    switch (plate) {
    case 0:
        // nGrid = nR; gridMin = Rmin; gridMax = Rmax
        Z0 = -0.24; // zFloor
        break;
    default:
        // nGrid = nZ; gridMin = Zmin; gridMax = Zmax
        R0 = 0.435;               // rFloor
        tracer.inverse_map(true); // mat must be inverted
        break;
    }

// loop over the grid
#pragma omp for schedule(dynamic)
    for (int i = 0; i < nPhi; i++) {
        for (int j = 0; j < nGrid; j++) {

            if (plate == 0) {
                R0 = gridMin + (gridMax - gridMin) * j / nGrid;
            } else {
                Z0 = gridMin + (gridMax - gridMin) * j / nGrid;
            }

            double phi0 = 2 * M_PI * i / nPhi;
            double R1, phi1, Z1;

            // evolve lines
            tracer.alloc_hint();
            evolve_line(tracer, R0, Z0, phi0, R1, Z1, phi1, scalars);
            tracer.clear_hint();

            // Calculate the index for dataWrite
            int index = i * nGrid + j;

            // Assign the values directly into the preallocated matrix
            outputData[index][0] = R0;
            outputData[index][1] = Z0;
            outputData[index][2] = phi0;
            outputData[index][3] = R1;
            outputData[index][4] = Z1;
            outputData[index][5] = phi1;
            outputData[index][6] = scalars.deltaPhi;
            outputData[index][7] = scalars.length;
            outputData[index][8] = scalars.psimin;

            // Optional: Progress indication
            printf("Thread %d: progress = %f\n", omp_get_thread_num(), 100.0 * (i * nGrid + j) / (nGrid * nPhi));
        }
    }
}

void footprint::evolve_line(maglit &tracer, double R0, double Z0, double phi0, double &R1, double &Z1, double &phi1, map_scalars &scalars) {
    R1 = R0;
    Z1 = Z0;
    phi1 = phi0;
    int status;
    double phi_max = 100 * 2 * M_PI;
    double arc = 0;
    double psin0 = 5.0;
    double *psin1 = &psin0;
    scalars.psimin = *psin1;
    tracer.reset();
    do {
        R0 = R1;
        Z0 = Z1;
        phi0 = phi1;
        status = tracer.step(R1, Z1, phi1, phi_max, -1);
        if (status == SODE_CONTINUE_GOOD_STEP) {
            arc += connection_length(R0, Z0, phi0, R1, Z1, phi1);
            tracer.psin_eval(R1, phi1, Z1, psin1);
            if (*psin1 < scalars.psimin) {
                scalars.psimin = *psin1;
            }
        }
    } while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP);
    scalars.deltaPhi = phi1 - phi0;
    scalars.length = arc;
}

double footprint::connection_length(double R0, double Z0, double phi0, double R1, double Z1, double phi1) {
    double dx = R1 * cos(phi1) - R0 * cos(phi0);
    double dy = R1 * sin(phi1) - R0 * sin(phi0);
    double dz = Z1 - Z0;
    return sqrt(dx * dx + dy * dy + dz * dz);
}