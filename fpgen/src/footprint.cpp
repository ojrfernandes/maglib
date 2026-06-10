#include "footprint.h"
#include <fstream>
#include <stdexcept>

footprint::footprint(const int manifold, const double grid_R1, const double grid_Z1, const double grid_R2, const double grid_Z2, const int nRZ, const int nPhi, const int max_turns) : manifold(manifold), grid_R1(grid_R1), grid_Z1(grid_Z1), grid_R2(grid_R2), grid_Z2(grid_Z2), nRZ(nRZ), nPhi(nPhi), max_turns(max_turns) {
    if (nRZ < 1 || nPhi < 1)
        throw std::invalid_argument("nRZ and nPhi must each be at least 1");
    this->outputData.resize(nRZ * nPhi, std::vector<double>(6));
}

const std::vector<std::vector<double>> &footprint::get_output_data() const {
    return outputData;
}

void footprint::run(std::vector<maglit*> &tracers) {
    int nthreads = static_cast<int>(tracers.size());

    // set manifold direction on all tracers
    for (auto *t : tracers) {
        if (manifold == 1) {
            t->inverse_map(false); // unstable: forward map
        } else {
            t->inverse_map(true);  // stable: inverse map
        }
    }

#pragma omp parallel num_threads(nthreads)
    {
        int     tid    = omp_get_thread_num();
        maglit &tracer = *tracers[tid];

#pragma omp for schedule(dynamic)
        for (int i = 0; i < nPhi; i++) {
            for (int j = 0; j < nRZ; j++) {

                double t        = (nRZ > 1) ? static_cast<double>(j) / (nRZ - 1) : 0.0;
                double R_init   = grid_R1 + (grid_R2 - grid_R1) * t;
                double Z_init   = grid_Z1 + (grid_Z2 - grid_Z1) * t;
                double phi_init = 2 * M_PI * i / nPhi;

                double      R_final, Z_final, phi_final;
                map_scalars scalars; // local to each iteration — thread-safe

                this->evolve_line(tracer, R_init, Z_init, phi_init,
                                  R_final, Z_final, phi_final, scalars);

                int index          = i * nRZ + j;
                outputData[index][0] = R_init;
                outputData[index][1] = Z_init;
                outputData[index][2] = phi_init;
                outputData[index][3] = scalars.length;
                outputData[index][4] = scalars.psimin;
                outputData[index][5] = scalars.turn;

                if (tid == 0) {
                    float progress = (float)(i * nRZ + j) / (nRZ * nPhi);
                    this->progressBar(progress);
                }
            }
        }
    }
}

void footprint::evolve_line(maglit &tracer, double R0, double Z0, double phi0, double &R1, double &Z1, double &phi1, map_scalars &scalars) {
    R1 = R0;
    Z1 = Z0;
    phi1 = phi0;

    // integration variables
    int    status;
    double phi_max = this->max_turns * 2 * M_PI;
    double arc = 0;

    double psin_current = std::numeric_limits<double>::max();
    scalars.psimin = psin_current;

    const double phi_start = phi0;
    scalars.turn = 0;
    tracer.reset();
    do {
        R0 = R1;
        Z0 = Z1;
        phi0 = phi1;
        status = tracer.step(R1, Z1, phi1, phi_max, -1);
        if (status == SODE_CONTINUE_GOOD_STEP) {
            arc += this->connection_length(R0, Z0, phi0, R1, Z1, phi1);
            if (tracer.psin_eval(R1, phi1, Z1, &psin_current)) {
                if (psin_current < scalars.psimin)
                    scalars.psimin = psin_current;
            }
        }
    } while (status == SODE_CONTINUE_GOOD_STEP || status == SODE_CONTINUE_BAD_STEP);

    scalars.turn = static_cast<int>(floor((phi1 - phi_start) / (2 * M_PI)));
    scalars.length = arc;
}

double footprint::connection_length(double R0, double Z0, double phi0, double R1, double Z1, double phi1) {
    double dx = R1 * cos(phi1) - R0 * cos(phi0);
    double dy = R1 * sin(phi1) - R0 * sin(phi0);
    double dz = Z1 - Z0;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

bool footprint::save(const std::string &path) const {
    auto ends_with = [&](const std::string &ext) {
        return path.size() >= ext.size() &&
               path.substr(path.size() - ext.size()) == ext;
    };

    bool is_csv = ends_with(".csv");
    bool is_txt = ends_with(".dat") || ends_with(".txt");

    if (!is_csv && !is_txt) {
        std::cerr << "Error: unsupported extension in \"" << path
                  << "\". Use .dat, .txt, or .csv." << std::endl;
        return false;
    }

    std::ofstream f(path);
    if (!f.is_open()) {
        std::cerr << "Error: could not open output file: " << path << std::endl;
        return false;
    }

    if (is_csv) {
        f << "R0,Z0,phi0,length,psiMin,turn\n";
        for (const auto &row : outputData) {
            f << std::fixed << std::setprecision(16)
              << row[0] << "," << row[1] << "," << row[2] << ","
              << row[3] << "," << row[4] << ","
              << std::setprecision(0) << row[5] << "\n";
        }
    } else {
        f << "#R0" << std::string(17, ' ') << "Z0" << std::string(17, ' ')
          << "phi0" << std::string(15, ' ') << "length" << std::string(13, ' ')
          << "psiMin" << std::string(13, ' ') << "turn\n";
        for (const auto &row : outputData) {
            f << std::fixed << std::setprecision(16)
              << row[0] << " " << row[1] << " " << row[2]
              << " " << row[3] << " " << row[4] << " "
              << std::setprecision(0) << row[5] << "\n";
        }
    }

    return true;
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