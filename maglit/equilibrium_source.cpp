#include "equilibrium_source.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

namespace {

constexpr double kTwoPi = 2.0 * M_PI;

// Natural cubic spline second derivatives of y(x), x strictly increasing, shape=(n,).
std::vector<double> spline_y2(const std::vector<double> &x, const std::vector<double> &y) {
    const size_t n = x.size();
    std::vector<double> y2(n, 0.0), u(n, 0.0);
    for (size_t i = 1; i + 1 < n; ++i) {
        double sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
        double p   = sig * y2[i - 1] + 2.0;
        y2[i]      = (sig - 1.0) / p;
        double d   = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
        u[i]       = (6.0 * d / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
    }
    y2[n - 1] = 0.0;
    for (size_t k = n - 1; k-- > 0;)
        y2[k] = y2[k] * y2[k + 1] + u[k];
    return y2;
}

// Spline value at xv; clamped to the end values outside [x0, x_{n-1}].
double spline_eval(const std::vector<double> &x, const std::vector<double> &y,
                   const std::vector<double> &y2, double xv) {
    if (xv <= x.front()) return y.front();
    if (xv >= x.back()) return y.back();
    size_t k = std::upper_bound(x.begin(), x.end(), xv) - x.begin() - 1;
    double h = x[k + 1] - x[k];
    double a = (x[k + 1] - xv) / h, b = 1.0 - a;
    return a * y[k] + b * y[k + 1] + ((a * a * a - a) * y2[k] + (b * b * b - b) * y2[k + 1]) * h * h / 6.0;
}

// Derivative of the natural spline through n equally spaced samples y[0], y[s], y[2s], ...
// (spacing h) at every sample; written to dy[0], dy[s], ...
void spline_node_deriv(const double *y, double *dy, int n, size_t s, double h) {
    std::vector<double> x(n), v(n);
    for (int i = 0; i < n; ++i) {
        x[i] = i * h;
        v[i] = y[i * s];
    }
    std::vector<double> y2 = spline_y2(x, v);
    for (int i = 0; i + 1 < n; ++i)
        dy[i * s] = (v[i + 1] - v[i]) / h - h * (2.0 * y2[i] + y2[i + 1]) / 6.0;
    dy[(n - 1) * s] = (v[n - 1] - v[n - 2]) / h + h * (y2[n - 2] + 2.0 * y2[n - 1]) / 6.0;
}

// Cubic Hermite basis on t in [0,1]: h[0], h[1] value weights (left, right node),
// h[2], h[3] slope weights (left, right; multiply by the cell width). dh = d/dt.
void hermite(double t, double h[4], double dh[4]) {
    double t2 = t * t, t3 = t2 * t;
    h[0]  = 2 * t3 - 3 * t2 + 1;
    h[1]  = -2 * t3 + 3 * t2;
    h[2]  = t3 - 2 * t2 + t;
    h[3]  = t3 - t2;
    dh[0] = 6 * t2 - 6 * t;
    dh[1] = -6 * t2 + 6 * t;
    dh[2] = 3 * t2 - 4 * t + 1;
    dh[3] = 3 * t2 - 2 * t;
}

bool parse_numbers(const std::string &s, std::vector<double> &out) {
    out.clear();
    std::string t = s;
    std::replace(t.begin(), t.end(), ',', ' ');
    const char *p = t.c_str();
    char       *end;
    while (true) {
        while (*p && std::isspace(static_cast<unsigned char>(*p))) ++p;
        if (!*p) break;
        double v = std::strtod(p, &end);
        if (end == p) return false;
        out.push_back(v);
        p = end;
    }
    return true;
}

// Text content of the first element `name` in s (matches "<name>" or "<name ...>").
bool xml_text(const std::string &s, const std::string &name, std::string &out) {
    size_t pos = 0;
    while ((pos = s.find("<" + name, pos)) != std::string::npos) {
        char c = s[pos + 1 + name.size()];
        if (c == '>' || std::isspace(static_cast<unsigned char>(c))) {
            size_t open  = s.find('>', pos);
            size_t close = s.find("</" + name + ">", open);
            if (open == std::string::npos || close == std::string::npos) return false;
            out = s.substr(open + 1, close - open - 1);
            return true;
        }
        pos += 1;
    }
    return false;
}

bool xml_numbers(const std::string &s, const std::string &name, std::vector<double> &out) {
    std::string text;
    if (!xml_text(s, name, text) || !parse_numbers(text, out) || out.empty()) {
        std::cerr << "EquilibriumSource: missing or unreadable XML element <" << name << ">" << std::endl;
        return false;
    }
    return true;
}

// True if v is uniformly spaced to within a relative tolerance of the mean step.
bool is_uniform(const std::vector<double> &v) {
    if (v.size() < 4) return false;
    double h = (v.back() - v.front()) / (v.size() - 1);
    if (!(h > 0)) return false;
    for (size_t i = 1; i < v.size(); ++i)
        if (std::fabs(v[i] - v[i - 1] - h) > 1e-4 * h) return false;
    return true;
}

} // namespace

EquilibriumSource::EquilibriumSource(const std::string &path) {
    std::ifstream f(path);
    if (!f) {
        std::cerr << "EquilibriumSource: cannot open " << path << std::endl;
        return;
    }
    std::stringstream ss;
    ss << f.rdbuf();
    const std::string text = ss.str();

    size_t first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        std::cerr << "EquilibriumSource: empty file " << path << std::endl;
        return;
    }
    valid_ = (text[first] == '<') ? load_xml(text) : load_geqdsk(text);
    if (!valid_)
        std::cerr << "EquilibriumSource: failed to load " << path << std::endl;
}

bool EquilibriumSource::load_xml(const std::string &text) {
    std::string eq;
    if (!xml_text(text, "equilibrium", eq)) {
        std::cerr << "EquilibriumSource: no <equilibrium> element" << std::endl;
        return false;
    }
    std::vector<double> Ri, zj, pfm, op, xp, psin, jpol;
    if (!xml_numbers(eq, "Ri", Ri) || !xml_numbers(eq, "zj", zj) || !xml_numbers(eq, "PFM", pfm) ||
        !xml_numbers(eq, "op", op) || !xml_numbers(eq, "xp", xp) ||
        !xml_numbers(eq, "psin", psin) || !xml_numbers(eq, "Jpol", jpol))
        return false;

    if (!is_uniform(Ri) || !is_uniform(zj)) {
        std::cerr << "EquilibriumSource: Ri/zj grid is not uniform" << std::endl;
        return false;
    }
    if (pfm.size() != Ri.size() * zj.size()) {
        std::cerr << "EquilibriumSource: PFM has " << pfm.size() << " values, not "
                  << zj.size() << "x" << Ri.size() << std::endl;
        return false;
    }
    if (op.size() != 3 || xp.size() % 3 != 0 || psin.size() != jpol.size() || psin.size() < 3) {
        std::cerr << "EquilibriumSource: inconsistent op/xp/psin/Jpol sizes" << std::endl;
        return false;
    }

    nR_ = static_cast<int>(Ri.size());
    nZ_ = static_cast<int>(zj.size());
    R0_ = Ri.front();
    dR_ = (Ri.back() - Ri.front()) / (nR_ - 1);
    Z0_ = zj.front();
    dZ_ = (zj.back() - zj.front()) / (nZ_ - 1);
    psi_    = std::move(pfm); // PFM dim "nz,nR", C order == our (nZ,nR) row-major, already in Wb
    sigma_  = 1.0;
    build_derivatives();

    psi_axis_ = op[0];
    R_axis_   = op[1];
    Z_axis_   = op[2];
    const double Z_axis = op[2];

    // xp: rows (psi, R, Z), one column per entry; keep genuine saddles only
    const size_t nx  = xp.size() / 3;
    bool         got = false;
    z_lo_            = -std::numeric_limits<double>::infinity();
    z_hi_            = std::numeric_limits<double>::infinity();
    for (size_t j = 0; j < nx; ++j) {
        double px = xp[j], Rx = xp[nx + j], Zx = xp[2 * nx + j];
        double p, pR, pZ;
        if (!eval_flux(Rx, Zx, p, pR, pZ)) continue;
        if (std::hypot(pR, pZ) >= 1e-2 * std::fabs(psi_axis_ - px)) continue; // not a saddle
        if (!got) {
            psi_bdry_ = px;
            got       = true;
        }
        if (Zx < Z_axis) z_lo_ = std::max(z_lo_, Zx);
        else             z_hi_ = std::min(z_hi_, Zx);
    }
    if (!got) {
        std::cerr << "EquilibriumSource: no X-point (saddle of PFM) in singularities/xp; "
                     "separatrix flux unknown" << std::endl;
        return false;
    }

    for (size_t i = 1; i < psin.size(); ++i)
        if (!(psin[i] > psin[i - 1])) {
            std::cerr << "EquilibriumSource: profiles1D/psin is not increasing" << std::endl;
            return false;
        }
    fx_ = psin;
    fy_.resize(jpol.size());
    for (size_t i = 0; i < jpol.size(); ++i)
        fy_[i] = 2e-7 * jpol[i]; // mu0 Jpol / (2 pi) [T m]
    fy2_   = spline_y2(fx_, fy_);
    F_vac_ = fy_.back();
    return true;
}

bool EquilibriumSource::load_geqdsk(const std::string &text) {
    std::vector<std::string> lines;
    {
        std::istringstream is(text);
        std::string        l;
        while (std::getline(is, l)) {
            if (!l.empty() && l.back() == '\r') l.pop_back();
            lines.push_back(l);
        }
    }
    if (lines.empty()) return false;

    // Header: the last two integers are nw, nh
    std::vector<std::string> tok;
    {
        std::istringstream is(lines[0]);
        std::string        t;
        while (is >> t) tok.push_back(t);
    }
    if (tok.size() < 2) return false;
    const int nw = std::atoi(tok[tok.size() - 2].c_str());
    const int nh = std::atoi(tok[tok.size() - 1].c_str());
    if (nw < 4 || nh < 4) {
        std::cerr << "EquilibriumSource: bad G-EQDSK header (nw=" << nw << ", nh=" << nh << ")" << std::endl;
        return false;
    }

    // Fixed-width 16-character fields; each block starts on a new line.
    size_t li         = 1;
    auto   read_block = [&](size_t count, std::vector<double> &out) -> bool {
        out.clear();
        while (out.size() < count) {
            if (li >= lines.size()) return false;
            const std::string &l = lines[li++];
            for (size_t k = 0; k + 16 <= l.size() && out.size() < count; k += 16) {
                std::string chunk = l.substr(k, 16);
                if (chunk.find_first_not_of(' ') == std::string::npos) continue;
                char  *end;
                double v = std::strtod(chunk.c_str(), &end);
                if (end == chunk.c_str()) return false;
                out.push_back(v);
            }
        }
        return true;
    };

    std::vector<double> b1, b2, b3, b4, fpol, pres, ffp, pp, psirz, qpsi;
    if (!read_block(5, b1) || !read_block(5, b2) || !read_block(5, b3) || !read_block(5, b4) ||
        !read_block(nw, fpol) || !read_block(nw, pres) || !read_block(nw, ffp) || !read_block(nw, pp) ||
        !read_block(size_t(nw) * nh, psirz) || !read_block(nw, qpsi)) {
        std::cerr << "EquilibriumSource: truncated or unreadable G-EQDSK" << std::endl;
        return false;
    }
    const double rdim = b1[0], zdim = b1[1], rleft = b1[3], zmid = b1[4];
    const double simag = b2[2], sibry = b2[3], current = b3[0];

    // Boundary (optional): "nbbbs limitr" then rbbbs/zbbbs interleaved
    z_lo_ = -std::numeric_limits<double>::infinity();
    z_hi_ = std::numeric_limits<double>::infinity();
    int nbbbs = 0, limitr = 0;
    if (li < lines.size()) {
        std::istringstream is(lines[li++]);
        is >> nbbbs >> limitr;
    }
    std::vector<double> bbbs;
    if (nbbbs > 0 && read_block(2 * size_t(nbbbs), bbbs)) {
        z_lo_ = std::numeric_limits<double>::infinity();
        z_hi_ = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < nbbbs; ++i) {
            z_lo_ = std::min(z_lo_, bbbs[2 * i + 1]);
            z_hi_ = std::max(z_hi_, bbbs[2 * i + 1]);
        }
    } else {
        std::cerr << "EquilibriumSource: warning: no boundary in G-EQDSK; "
                     "private-flux region not excluded from F(psi_N)" << std::endl;
    }

    if (current == 0.0 || simag == sibry) {
        std::cerr << "EquilibriumSource: G-EQDSK has Ip = 0 or simag = sibry" << std::endl;
        return false;
    }

    nR_ = nw;
    nZ_ = nh;
    R0_ = rleft;
    dR_ = rdim / (nw - 1);
    Z0_ = zmid - 0.5 * zdim;
    dZ_ = zdim / (nh - 1);
    // shortcut: assumes psirz in Wb/rad and phi counter-clockwise (COCOS 1-8, odd frame);
    // add an explicit COCOS input if a COCOS 11+ or clockwise-phi file has to be read.
    psi_.resize(psirz.size());
    for (size_t i = 0; i < psirz.size(); ++i)
        psi_[i] = kTwoPi * psirz[i]; // psirz[iz*nw + ir] == our (nZ,nR) row-major
    psi_axis_ = kTwoPi * simag;
    R_axis_   = b2[0]; // rmaxis
    Z_axis_   = b2[1]; // zmaxis
    psi_bdry_ = kTwoPi * sibry;
    // B_Z at the outboard midplane must have the sign of -Ip (Ampere, phi counter-clockwise)
    sigma_ = ((current > 0) == (simag > sibry)) ? 1.0 : -1.0;
    build_derivatives();

    fx_.resize(nw);
    for (int i = 0; i < nw; ++i)
        fx_[i] = double(i) / (nw - 1);
    fy_    = fpol;
    fy2_   = spline_y2(fx_, fy_);
    F_vac_ = fpol.back();
    return true;
}

void EquilibriumSource::build_derivatives() {
    const size_t n = psi_.size();
    psiR_.assign(n, 0.0);
    psiZ_.assign(n, 0.0);
    psiRZ_.assign(n, 0.0);
    for (int k = 0; k < nZ_; ++k)
        spline_node_deriv(&psi_[k * nR_], &psiR_[k * nR_], nR_, 1, dR_);
    for (int i = 0; i < nR_; ++i) {
        spline_node_deriv(&psi_[i], &psiZ_[i], nZ_, nR_, dZ_);
        spline_node_deriv(&psiR_[i], &psiRZ_[i], nZ_, nR_, dZ_);
    }
}

bool EquilibriumSource::eval_flux(double R, double Z, double &psi, double &psi_R, double &psi_Z) const {
    double x = (R - R0_) / dR_, y = (Z - Z0_) / dZ_;
    if (!(x >= 0 && x <= nR_ - 1 && y >= 0 && y <= nZ_ - 1)) return false; // also rejects NaN
    int i = std::min(static_cast<int>(x), nR_ - 2);
    int k = std::min(static_cast<int>(y), nZ_ - 2);
    double hr[4], dhr[4], hz[4], dhz[4];
    hermite(x - i, hr, dhr);
    hermite(y - k, hz, dhz);

    psi = psi_R = psi_Z = 0.0;
    for (int b = 0; b < 2; ++b) {
        for (int a = 0; a < 2; ++a) {
            size_t idx = size_t(k + b) * nR_ + (i + a);
            double p = psi_[idx], pR = psiR_[idx] * dR_, pZ = psiZ_[idx] * dZ_, pRZ = psiRZ_[idx] * dR_ * dZ_;
            double cz = p * hz[b] + pZ * hz[2 + b], czr = pR * hz[b] + pRZ * hz[2 + b];
            double dz = p * dhz[b] + pZ * dhz[2 + b], dzr = pR * dhz[b] + pRZ * dhz[2 + b];
            psi   += hr[a] * cz + hr[2 + a] * czr;
            psi_R += dhr[a] * cz + dhr[2 + a] * czr;
            psi_Z += hr[a] * dz + hr[2 + a] * dzr;
        }
    }
    psi_R /= dR_;
    psi_Z /= dZ_;
    return true;
}

bool EquilibriumSource::eval_B(double R, double /*phi*/, double Z, double B[3]) {
    double p, pR, pZ;
    if (!valid_ || !eval_flux(R, Z, p, pR, pZ)) return false;
    double psin   = (p - psi_axis_) / (psi_bdry_ - psi_axis_);
    bool   inside = psin < 1.0 && Z > z_lo_ && Z < z_hi_;
    double F      = inside ? spline_eval(fx_, fy_, fy2_, psin) : F_vac_;
    B[0]          = -sigma_ * pZ / (kTwoPi * R);
    B[1]          = F / R;
    B[2]          = sigma_ * pR / (kTwoPi * R);
    return true;
}

bool EquilibriumSource::eval_psin(double R, double /*phi*/, double Z, double &psin) {
    double p, pR, pZ;
    if (!valid_ || !eval_flux(R, Z, p, pR, pZ)) return false;
    psin = (p - psi_axis_) / (psi_bdry_ - psi_axis_);
    return true;
}

bool EquilibriumSource::eval_psi(double R, double /*phi*/, double Z, double &psi) {
    double pR, pZ;
    return valid_ && eval_flux(R, Z, psi, pR, pZ);
}
