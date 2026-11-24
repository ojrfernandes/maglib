"""
Tests for EquilibriumSource (XML / G-EQDSK), GridPerturbationSource and PerturbedSource,
on synthetic files written to tmp_path (no external data needed).

Analytic equilibrium, flux per radian (COCOS-1 like, Ip > 0, phi counter-clockwise):
    psi(R, Z) = (R - R0)^2 + Z^2 + C Z^3        axis (R0, 0), saddle (R0, -2/(3C))
    B_R = psi_Z / R,  B_Z = -psi_R / R,  B_phi = F(psi_N) / R
    F = F0 (1 + EPS (1 - psi_N)) inside the separatrix, F0 elsewhere.
"""
import math

import numpy as np
import pytest

from maglib import EquilibriumSource, Footprint, GridPerturbationSource, Maglit, PerturbedSource

R0, C, F0, EPS = 1.6, 1.0, -3.0, 0.05
ZX = -2.0 / (3.0 * C)  # m, X-point Z
PSI_X = ZX**2 + C * ZX**3  # Wb/rad, separatrix flux (axis flux = 0)
# top of the separatrix on R = R0: Z^2 + C Z^3 = PSI_X, Z > 0
ZTOP = float(next(r.real for r in np.roots([C, 1.0, 0.0, -PSI_X]) if abs(r.imag) < 1e-12 and r.real > 0))
R1 = np.linspace(1.0, 2.2, 121)  # m, shape=(nR,)
Z1 = np.linspace(-1.0, 0.8, 181)  # m, shape=(nZ,)
MU0 = 4e-7 * math.pi


def psi_rad(R, Z):
    return (R - R0) ** 2 + Z**2 + C * Z**3


def b_exact(R, Z):
    """(B_R, B_phi, B_Z) [T] of the analytic equilibrium."""
    psiR, psiZ = 2.0 * (R - R0), 2.0 * Z + 3.0 * C * Z**2
    psin = psi_rad(R, Z) / PSI_X
    inside = psin < 1.0 and ZX < Z < ZTOP
    F = F0 * (1.0 + EPS * (1.0 - psin)) if inside else F0
    return np.array([psiZ / R, F / R, -psiR / R])


def _fmt5(vals):
    vals = list(vals)
    return "".join(
        "".join(f"{v:16.9E}" for v in vals[i : i + 5]) + "\n" for i in range(0, len(vals), 5)
    )


@pytest.fixture(scope="module")
def geqdsk_path(tmp_path_factory):
    nw, nh = R1.size, Z1.size
    RR, ZZ = np.meshgrid(R1, Z1)  # shape=(nZ,nR)
    psirz = psi_rad(RR, ZZ)
    psin = np.linspace(0.0, 1.0, nw)
    fpol = F0 * (1.0 + EPS * (1.0 - psin))
    zero = np.zeros(nw)
    txt = f"  SYNTHETIC TEST                                  0 {nw:3d} {nh:3d}\n"
    txt += _fmt5([R1[-1] - R1[0], Z1[-1] - Z1[0], R0, R1[0], 0.5 * (Z1[0] + Z1[-1])])
    txt += _fmt5([R0, 0.0, 0.0, PSI_X, F0 / R0])
    txt += _fmt5([1e6, 0.0, 0.0, R0, 0.0])
    txt += _fmt5([0.0, 0.0, PSI_X, 0.0, 0.0])
    for arr in (fpol, zero, zero, zero, psirz.ravel(), zero):
        txt += _fmt5(arr)
    txt += f"{2:5d}{0:5d}\n"
    txt += _fmt5([R0, ZX, R0, ZTOP])
    p = tmp_path_factory.mktemp("eq") / "g_eqdsk"
    p.write_text(txt)
    return str(p)


@pytest.fixture(scope="module")
def xml_path(tmp_path_factory):
    # XML convention: B_R = -(1/2piR) dPFM/dZ, so the same field needs PFM = -2 pi psi_rad
    RR, ZZ = np.meshgrid(R1, Z1)
    pfm = -2 * math.pi * psi_rad(RR, ZZ)  # shape=(nZ,nR), C order
    psin = np.append(np.linspace(0.0, 1.0, 51), 1.2)
    jpol = F0 * (1.0 + EPS * (1.0 - np.minimum(psin, 1.0))) / 2e-7  # A, F = mu0 Jpol / 2pi
    xp = np.array([[-2 * math.pi * PSI_X, -2 * math.pi * psi_rad(1.2, 0.0)],  # psi row
                   [R0, 1.2],  # R row; column 2 is not a saddle
                   [ZX, 0.0]])  # Z row
    j = lambda a: ", ".join(f"{v:.17e}" for v in np.ravel(a))
    txt = f"""<?xml version="1.0" ?>
<magnetics version="2.3">
  <equilibrium R0="{R0}" Ip="1e6" B0="{F0 / R0}">
    <singularities>
      <op unit="T m^2,m,m" dim="3" type="float64">0.0, {R0}, 0.0</op>
      <xp unit="T m^2,m,m" dim="3,2" type="float64">{j(xp)}</xp>
    </singularities>
    <profiles1D>
      <psin unit="1" dim="{psin.size}" type="float64">{j(psin)}</psin>
      <Jpol unit="A" dim="{psin.size}" type="float64">{j(jpol)}</Jpol>
      <Jpolprime unit="A" dim="{psin.size}" type="float64">{j(np.zeros_like(psin))}</Jpolprime>
    </profiles1D>
    <Ri unit="m" dim="{R1.size}" type="float64">{j(R1)}</Ri>
    <zj unit="m" dim="{Z1.size}" type="float64">{j(Z1)}</zj>
    <PFM unit="T m^2" dim="{Z1.size},{R1.size}" type="float64">{j(pfm)}</PFM>
  </equilibrium>
</magnetics>
"""
    p = tmp_path_factory.mktemp("eq") / "equilibrium.xml"
    p.write_text(txt)
    return str(p)


# interior points (away from the natural-spline end conditions), inside and outside the separatrix,
# plus the private-flux region below the X-point
PTS = [(1.45, 0.1), (1.75, -0.3), (1.9, 0.2), (1.3, -0.5), (1.6, -0.8), (1.55, -0.75), (2.0, 0.5)]


@pytest.mark.parametrize("which", ["geqdsk", "xml"])
def test_equilibrium_field_matches_analytic(which, geqdsk_path, xml_path):
    src = EquilibriumSource(geqdsk_path if which == "geqdsk" else xml_path)
    assert src.is_valid()
    for R, Z in PTS:
        ok, B = src.eval_B(R, 0.3, Z)
        assert ok
        ref = b_exact(R, Z)
        assert np.allclose(B, ref, rtol=1e-6, atol=1e-6 * np.abs(ref).max()), (R, Z, B, ref)


@pytest.mark.parametrize("which", ["geqdsk", "xml"])
def test_equilibrium_psin_and_domain(which, geqdsk_path, xml_path):
    src = EquilibriumSource(geqdsk_path if which == "geqdsk" else xml_path)
    ok, pa = src.eval_psin(R0, 0.0, 0.0)
    assert ok and abs(pa) < 1e-9
    ok, px = src.eval_psin(R0, 0.0, ZX)
    assert ok and abs(px - 1.0) < 1e-9
    ok, psi = src.eval_psi(1.8, 0.0, 0.1)
    sign = 1.0 if which == "geqdsk" else -1.0  # file's own sign convention, Wb
    assert ok and abs(psi - sign * 2 * math.pi * psi_rad(1.8, 0.1)) < 1e-9
    assert not src.eval_B(0.5, 0.0, 0.0)[0]
    assert not src.eval_B(1.5, 0.0, 1.5)[0]
    assert not src.eval_psin(float("nan"), 0.0, 0.0)[0]


@pytest.mark.parametrize("which", ["geqdsk", "xml"])
def test_equilibrium_magnetic_axis(which, geqdsk_path, xml_path):
    src = EquilibriumSource(geqdsk_path if which == "geqdsk" else xml_path)
    assert src.magnetic_axis() == pytest.approx((R0, 0.0), abs=1e-12)


def test_equilibrium_ampere_sign(geqdsk_path):
    # Ip > 0 along +phi (counter-clockwise from above) => B_Z < 0 on the outboard midplane
    ok, B = EquilibriumSource(geqdsk_path).eval_B(1.9, 0.0, 0.0)
    assert ok and B[2] < 0


def test_equilibrium_bad_files(tmp_path):
    assert not EquilibriumSource(str(tmp_path / "missing")).is_valid()
    p = tmp_path / "trunc"
    p.write_text("  HEADER 0 65 65\n 1.0 2.0\n")
    assert not EquilibriumSource(str(p)).is_valid()
    p = tmp_path / "noeq.xml"
    p.write_text("<?xml version='1.0'?><magnetics></magnetics>")
    assert not EquilibriumSource(str(p)).is_valid()


# ---------------------------------------------------------------------------------------------------
# gridded perturbation: period 2 pi / 3, analytic field
NTOR = 3
PER = 2 * math.pi / NTOR
NP, NRP, NZP = 36, 41, 41
RP = np.linspace(1.0, 2.0, NRP)  # m
ZP = np.linspace(-0.5, 0.5, NZP)  # m
FA = 0.1  # rad, phi of the first plane


def b_pert(R, phi, Z):
    """(B_phi, B_R, B_Z) [T], period 2 pi / NTOR."""
    a = (1 + 0.3 * R) * (1 - 0.2 * Z**2) * 1e-3
    return np.array([a * np.cos(NTOR * phi), a * np.sin(NTOR * phi + 0.3), a * R * Z * np.cos(NTOR * phi - 1.0)])


def write_dat(path, nc=3, fb=FA + PER, nbytes_cut=0):
    ph = FA + PER * np.arange(NP) / NP
    P, RR, ZZ = np.meshgrid(ph, RP, ZP, indexing="ij")  # shape=(Np,Nr,Nz)
    ba = np.zeros((nc, NP, NRP, NZP))
    ba[:3] = b_pert(RR, P, ZZ)
    raw = np.array([NP + 1, NRP, NZP, nc], "<i4").tobytes()
    raw += np.array([RP[0], RP[-1], ZP[0], ZP[-1], FA, fb], "<f8").tobytes()
    raw += np.asarray(ba, "<f8").tobytes(order="F")
    path.write_bytes(raw[: len(raw) - nbytes_cut])
    return str(path)


@pytest.mark.parametrize("nc", [3, 12])
def test_perturbation_nodes_periodicity_accuracy(tmp_path, nc):
    src = GridPerturbationSource(write_dat(tmp_path / "p.dat", nc))
    assert src.is_valid() and abs(src.period() - PER) < 1e-12
    # exact at a node (stored order B_phi, B_R, B_Z -> returned B_R, B_phi, B_Z)
    R, Z, phi = RP[7], ZP[30], FA + PER * 5 / NP
    ref = b_pert(R, phi, Z)[[1, 0, 2]]
    ok, B = src.eval_B(R, phi, Z)
    assert ok and np.allclose(B, ref, rtol=0, atol=1e-15)
    # off-node: periodic copies agree, and the interpolant is close to the analytic field
    R, Z, phi = 1.513, 0.137, 0.777
    ok, B = src.eval_B(R, phi, Z)
    assert ok
    for k in (1, -5, 7):
        okk, Bk = src.eval_B(R, phi + k * PER, Z)
        assert okk and np.allclose(Bk, B, rtol=0, atol=1e-14)
    ref = b_pert(R, phi, Z)[[1, 0, 2]]
    assert np.abs(B - ref).max() < 5e-3 * np.abs(ref).max()
    # outside the R-Z box
    assert not src.eval_B(0.99, 0.0, 0.0)[0]
    assert not src.eval_B(1.5, 0.0, 0.51)[0]
    assert not src.eval_psin(1.5, 0.0, 0.0)[0]


def test_perturbation_bad_files(tmp_path):
    assert not GridPerturbationSource(str(tmp_path / "missing.dat")).is_valid()
    assert not GridPerturbationSource(write_dat(tmp_path / "t.dat", nbytes_cut=8)).is_valid()
    assert not GridPerturbationSource(write_dat(tmp_path / "b.dat", fb=FA + 2.5)).is_valid()


def test_perturbed_source(tmp_path, geqdsk_path):
    pert_path = write_dat(tmp_path / "p.dat")
    amp, delta = 2.5, 0.4
    src = PerturbedSource(geqdsk_path, pert_path, amplitude=amp, phase_shift=delta)
    eq = EquilibriumSource(geqdsk_path)
    pert = GridPerturbationSource(pert_path)
    R, phi, Z = 1.7, 1.234, 0.05
    ok, B = src.eval_B(R, phi, Z)
    assert ok
    expected = eq.eval_B(R, phi, Z)[1] + amp * pert.eval_B(R, phi - delta, Z)[1]
    assert np.allclose(B, expected, rtol=1e-14, atol=0)
    assert src.eval_psin(R, phi, Z) == eq.eval_psin(R, phi, Z)
    # equilibrium grid is larger than the perturbation box: outside the box the sum fails
    assert eq.eval_B(2.1, 0.0, 0.0)[0] and not src.eval_B(2.1, 0.0, 0.0)[0]
    with pytest.raises(RuntimeError):
        PerturbedSource(geqdsk_path, str(tmp_path / "missing.dat"))


def test_footprint_needs_axis(tmp_path, geqdsk_path):
    """Perturbation-only source has no axis: Footprint.run raises unless mag_axis is given."""
    pert = GridPerturbationSource(write_dat(tmp_path / "p.dat"))
    assert pert.magnetic_axis() is None
    src = PerturbedSource(geqdsk_path, str(tmp_path / "p.dat"))
    assert src.magnetic_axis() == pytest.approx((R0, 0.0), abs=1e-12)

    t = Maglit(pert)
    t.configure(0.01, 1e-6, 0.1)
    with pytest.raises(RuntimeError, match="magnetic axis"):
        Footprint(0, 1.5, 0.0, 1.5, 0.0, 1, 1, 1).run([t])
