#include "footprint.h"

footprint::footprint(const double &plate, const double &gridMin, const double &gridMax, const int &nGrid, const int &nPhi) : plate(plate), nGrid(nGrid), nPhi(nPhi), gridMin(gridMin), gridMax(gridMax) {
}

void footprint::runGrid(maglit &tracer, std::ofstream &output_file) {
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

    // Thread-local buffer
    const size_t bufferSize = 500; // Flush every 500 entries
    std::vector<double> localBuffer;
    localBuffer.reserve(bufferSize * 5); // Reserve space for efficiency

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

            tracer.alloc_hint();
            this->evolve_line(tracer, R0, Z0, phi0, R1, Z1, phi1, scalars);
            tracer.clear_hint();

            // Append to the thread-local buffer
            localBuffer.push_back(R0);
            localBuffer.push_back(Z0);
            localBuffer.push_back(phi0);
            localBuffer.push_back(scalars.length);
            localBuffer.push_back(scalars.psimin);

            // Flush if buffer is full
            if (localBuffer.size() >= bufferSize * 5) {
                std::lock_guard<std::mutex> lock(file_mutex);
                for (size_t k = 0; k < localBuffer.size(); k += 5) {
                    output_file << std::fixed << std::setprecision(16)
                                << localBuffer[k] << " "
                                << localBuffer[k + 1] << " "
                                << localBuffer[k + 2] << " "
                                << localBuffer[k + 3] << " "
                                << localBuffer[k + 4] << "\n";
                }
                localBuffer.clear();
            }

            if (omp_get_thread_num() == 0) {
                float progress = (float)(i * nGrid + j) / (nGrid * nPhi);
                this->progressBar(progress);
                std::this_thread::sleep_for(std::chrono::milliseconds(50));
            }
        }
    }

    // Flush remaining buffer at the end of the loop
    std::lock_guard<std::mutex> lock(file_mutex);
    for (size_t k = 0; k < localBuffer.size(); k += 5) {
        output_file << std::fixed << std::setprecision(16)
                    << localBuffer[k] << " "
                    << localBuffer[k + 1] << " "
                    << localBuffer[k + 2] << " "
                    << localBuffer[k + 3] << " "
                    << localBuffer[k + 4] << "\n";
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