#include "grid_perturbation_source.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sys/mman.h>
#include <unistd.h>

namespace {

// Catmull-Rom weights for nodes -1, 0, 1, 2 at t in [0, 1].
void catmull_rom(double t, double w[4]) {
    double t2 = t * t, t3 = t2 * t;
    w[0] = 0.5 * (-t3 + 2 * t2 - t);
    w[1] = 0.5 * (3 * t3 - 5 * t2 + 2);
    w[2] = 0.5 * (-3 * t3 + 4 * t2 + t);
    w[3] = 0.5 * (t3 - t2);
}

} // namespace

GridPerturbationSource::GridPerturbationSource(const std::string &path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) {
        std::cerr << "GridPerturbationSource: cannot open " << path << std::endl;
        return;
    }
    int32_t hi[4];
    double  hd[6];
    f.read(reinterpret_cast<char *>(hi), sizeof hi);
    f.read(reinterpret_cast<char *>(hd), sizeof hd);
    if (!f) {
        std::cerr << "GridPerturbationSource: truncated header in " << path << std::endl;
        return;
    }
    f.close();

    const int32_t Nf = hi[0], Nr = hi[1], Nz = hi[2], Nc = hi[3];
    if (Nf < 2 || Nr < 2 || Nz < 2 || Nc < 3) {
        std::cerr << "GridPerturbationSource: bad header Nf=" << Nf << " Nr=" << Nr << " Nz=" << Nz
                  << " Ncomp=" << Nc << " in " << path << " (need Nf, Nr, Nz >= 2, Ncomp >= 3)" << std::endl;
        return;
    }
    const uint64_t expected = 64 + 8ull * Nc * (Nf - 1) * Nr * Nz;
    std::error_code ec;
    const uint64_t  size = std::filesystem::file_size(path, ec);
    if (ec || size != expected) {
        std::cerr << "GridPerturbationSource: file size " << size << " != expected " << expected
                  << " for header Nf=" << Nf << " Nr=" << Nr << " Nz=" << Nz << " Ncomp=" << Nc
                  << " in " << path << std::endl;
        return;
    }
    const double Ra = hd[0], Rb = hd[1], za = hd[2], zb = hd[3], fa = hd[4], fb = hd[5];
    const double P  = fb - fa;
    if (!(Rb > Ra) || !(zb > za) || !(P > 0)) {
        std::cerr << "GridPerturbationSource: bad grid extents in " << path << std::endl;
        return;
    }
    const double N = 2 * M_PI / P;
    if (std::fabs(N - std::round(N)) > 1e-6) {
        std::cerr << "GridPerturbationSource: toroidal period " << P << " rad does not divide 2 pi in "
                  << path << std::endl;
        return;
    }

    int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) {
        std::cerr << "GridPerturbationSource: cannot open " << path << std::endl;
        return;
    }
    void *m = ::mmap(nullptr, size, PROT_READ, MAP_SHARED, fd, 0);
    ::close(fd); // the mapping stays valid
    if (m == MAP_FAILED) {
        std::cerr << "GridPerturbationSource: mmap failed for " << path << std::endl;
        return;
    }

    map_      = m;
    map_size_ = size;
    data_     = reinterpret_cast<const double *>(static_cast<const char *>(m) + 64);
    Nc_ = Nc;
    Np_ = Nf - 1;
    Nr_ = Nr;
    Nz_ = Nz;
    Ra_ = Ra; Rb_ = Rb; za_ = za; zb_ = zb; fa_ = fa; P_ = P;
    dR_   = (Rb - Ra) / (Nr - 1);
    dZ_   = (zb - za) / (Nz - 1);
    dphi_ = P / Np_;
}

GridPerturbationSource::~GridPerturbationSource() {
    if (map_) ::munmap(map_, map_size_);
}

bool GridPerturbationSource::eval_B(double R, double phi, double Z, double B[3]) {
    if (!data_ || !(R >= Ra_ && R <= Rb_ && Z >= za_ && Z <= zb_)) return false; // also rejects NaN
    if (!std::isfinite(phi)) return false;

    double x = (R - Ra_) / dR_, y = (Z - za_) / dZ_;
    int    i = std::min(static_cast<int>(x), Nr_ - 2);
    int    k = std::min(static_cast<int>(y), Nz_ - 2);

    double u = std::fmod(phi - fa_, P_); // map into one stored period
    if (u < 0) u += P_;
    u /= dphi_;
    int ip = static_cast<int>(u);
    double tp = u - ip;
    if (ip >= Np_) ip -= Np_; // u == Np_ after rounding

    double wr[4], wz[4], wp[4];
    catmull_rom(x - i, wr);
    catmull_rom(y - k, wz);
    catmull_rom(tp, wp);

    int jr[4], jz[4], jp[4];
    for (int m = 0; m < 4; ++m) {
        jr[m] = std::clamp(i - 1 + m, 0, Nr_ - 1);
        jz[m] = std::clamp(k - 1 + m, 0, Nz_ - 1);
        jp[m] = (ip - 1 + m + Np_) % Np_;
    }

    double acc[3] = {0.0, 0.0, 0.0}; // (B_phi, B_R, B_Z) as stored
    for (int mz = 0; mz < 4; ++mz) {
        for (int mr = 0; mr < 4; ++mr) {
            const double  w    = wz[mz] * wr[mr];
            const double *base = data_ + size_t(Nc_) * Np_ * (jr[mr] + size_t(Nr_) * jz[mz]);
            for (int mp = 0; mp < 4; ++mp) {
                const double *v  = base + size_t(Nc_) * jp[mp];
                const double  ww = w * wp[mp];
                acc[0] += ww * v[0];
                acc[1] += ww * v[1];
                acc[2] += ww * v[2];
            }
        }
    }
    B[0] = acc[1];
    B[1] = acc[0];
    B[2] = acc[2];
    return true;
}
