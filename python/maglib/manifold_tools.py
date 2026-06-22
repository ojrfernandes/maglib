import numpy as np

try:
    from shapely.geometry import LineString
    _HAS_SHAPELY = True
except ImportError:
    _HAS_SHAPELY = False

_SHAPELY_MSG = "shapely is required. Install with: pip install shapely"


def lobe_map(equilibrium, perturbed, mag_axis, x_point=None):
    """
    Compute lobe boundary map between equilibrium and perturbed manifold curves.

    Finds all intersection points between the two curves, then for each pair of
    consecutive intersections computes the enclosed lobe's geometric properties.
    Intersection points are ordered by their position along the perturbed curve.

    Parameters
    ----------
    equilibrium : array-like of shape (N, 2)
        Equilibrium manifold curve, columns [R, Z] in metres.
    perturbed : array-like of shape (M, 2)
        Perturbed manifold curve, columns [R, Z] in metres.
    mag_axis : array-like of shape (2,)
        [R, Z] of the magnetic axis.
    x_point : array-like of shape (2,), optional
        [R, Z] of the X-point. Defaults to perturbed[0] when omitted.

    Returns
    -------
    ndarray of shape (N_lobes, 10)
        Columns: Rmid, Zmid, angle [rad], perimeter [m], area [m²], h_param [m],
        R1 [m], Z1 [m], R2 [m], Z2 [m].
        (R1, Z1) and (R2, Z2) are the two intersection points bounding each lobe,
        ordered along the perturbed curve.
        Returns shape (0, 10) when fewer than two intersections are found.
    """
    if not _HAS_SHAPELY:
        raise ImportError(_SHAPELY_MSG)

    from shapely.geometry import Polygon
    from shapely.ops import substring

    eq = np.asarray(equilibrium, dtype=float)
    ptb = np.asarray(perturbed, dtype=float)
    mag_axis = np.asarray(mag_axis, dtype=float)
    x_point = ptb[0].copy() if x_point is None else np.asarray(x_point, dtype=float)

    eq_ls = LineString(eq)
    ptb_ls = LineString(ptb)

    inter_geom = ptb_ls.intersection(eq_ls)
    if inter_geom.is_empty:
        return np.empty((0, 10))

    gtype = inter_geom.geom_type
    if gtype == 'Point':
        pts = [inter_geom]
    elif gtype == 'MultiPoint':
        pts = list(inter_geom.geoms)
    else:
        pts = [g for g in inter_geom.geoms if g.geom_type == 'Point']

    if len(pts) < 2:
        return np.empty((0, 10))

    pts = sorted(pts, key=lambda p: ptb_ls.project(p))

    v_xp = x_point - mag_axis
    mag_xp = np.linalg.norm(v_xp)

    rows = []
    for i in range(len(pts) - 1):
        p1, p2 = pts[i], pts[i + 1]

        t1_ptb = ptb_ls.project(p1)
        t2_ptb = ptb_ls.project(p2)
        t1_eq = eq_ls.project(p1)
        t2_eq = eq_ls.project(p2)

        # Perturbed sub-curve (always forward — pts are sorted by ptb position)
        ptb_coords = np.array(substring(ptb_ls, t1_ptb, t2_ptb).coords)

        # Equilibrium sub-curve (may be forward or backward along the curve)
        if t1_eq <= t2_eq:
            eq_coords = np.array(substring(eq_ls, t1_eq, t2_eq).coords)
        else:
            eq_coords = np.array(substring(eq_ls, t2_eq, t1_eq).coords)[::-1]

        # Closed lobe polygon: eq forward + ptb reversed, skipping shared endpoints
        inner_ptb = ptb_coords[-2:0:-1]
        lobe_pts = np.vstack([eq_coords, inner_ptb]) if len(inner_ptb) else eq_coords
        poly = Polygon(lobe_pts)

        mid_pt = eq_ls.interpolate(0.5 * (t1_eq + t2_eq))
        mid = np.array([mid_pt.x, mid_pt.y])

        # Reference angle: arc from mag_axis between midpoint and x_point
        v_mid = mid - mag_axis
        mag_mid = np.linalg.norm(v_mid)
        if mag_mid > 0.0 and mag_xp > 0.0:
            cos_a = np.clip(np.dot(v_mid, v_xp) / (mag_mid * mag_xp), -1.0, 1.0)
            angle = float(np.arccos(cos_a))
        else:
            angle = 0.0

        perimeter = poly.length
        area = poly.area
        base = np.hypot(p2.x - p1.x, p2.y - p1.y)
        h_param = 2.0 * area / base if base > 0.0 else 0.0

        rows.append([mid[0], mid[1], angle, perimeter, area, h_param,
                     p1.x, p1.y, p2.x, p2.y])

    return np.array(rows)


def simplify(source, tolerance):
    """
    Remove redundant points from manifold segments using Ramer-Douglas-Peucker.

    Points are only removed — no interpolation occurs, so no artificial
    information is introduced. The shape of each segment is preserved to within
    'tolerance' meters.

    Parameters
    ----------
    source : Manifold or list of ndarray of shape (N, 2)
        Manifold object or list of segments, each a (N, 2) array of [R, Z] in metres.
    tolerance : float
        Maximum allowed perpendicular deviation in metres. Larger values remove
        more points; smaller values are more conservative.

    Returns
    -------
    list of ndarray of shape (N_i, 2)
        Simplified segments with the same [R, Z] column layout.
    """
    if not _HAS_SHAPELY:
        raise ImportError(_SHAPELY_MSG)

    segments = _extract_segments(source)
    result = []
    for seg in segments:
        if len(seg) < 2:
            result.append(seg.copy())
            continue
        simplified = LineString(seg).simplify(tolerance, preserve_topology=False)
        result.append(np.array(simplified.coords))
    return result


def to_linestrings(source):
    """
    Convert manifold segments to shapely LineString objects.

    No interpolation is performed. Each LineString stores exactly the points
    from the corresponding segment, connected by straight edges.

    Parameters
    ----------
    source : Manifold or list of ndarray of shape (N, 2)
        Manifold object or list of (N, 2) [R, Z] arrays.

    Returns
    -------
    list of shapely.geometry.LineString
    """
    if not _HAS_SHAPELY:
        raise ImportError(_SHAPELY_MSG)

    return [LineString(seg) for seg in _extract_segments(source)]


def get_separatrix(hdf5_path, component=0, psi_n_level=0.99):
    """
    Extract the equilibrium separatrix from an M3DC1 HDF5 file.

    Reads the finite-element mesh and poloidal flux from the equilibrium group.
    No Fusion-IO evaluation is needed. Returns the longest connected contour
    segment at the requested ψ_N level.

    Parameters
    ----------
    hdf5_path : str or Path
        Path to the M3DC1 HDF5 file.
    component : int
        Toroidal Fourier component of psi to use. 0 = axisymmetric (n=0),
        which is the equilibrium flux. Default 0.
    psi_n_level : float
        Normalised poloidal flux level to contour. Default 0.99 rather than
        1.0 because the true separatrix (ψ_N=1) often does not close cleanly.

    Returns
    -------
    ndarray of shape (N, 2)
        Contour [R, Z] in metres at the requested ψ_N level, ordered along
        the curve. The longest connected segment is returned.
    """
    import h5py
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.tri import Triangulation

    with h5py.File(hdf5_path, 'r') as f:
        elements = f['equilibrium/mesh/elements'][:]
        psi      = f['equilibrium/fields/psi'][:]
        psi_lcfs = float(f['scalars/psi_lcfs'][0])
        psimin   = float(f['scalars/psimin'][0])

    R     = elements[:, 4].astype(float)
    Z     = elements[:, 5].astype(float)
    psi_N = (psi[:, component].astype(float) - psimin) / (psi_lcfs - psimin)

    tri = Triangulation(R, Z)
    fig = Figure()
    FigureCanvasAgg(fig)
    ax  = fig.add_subplot(111)
    cs  = ax.tricontour(tri, psi_N, levels=[psi_n_level])

    segs = cs.allsegs[0]
    if not segs:
        raise RuntimeError(
            f"No psi_N = {psi_n_level} contour found in the HDF5 file.")

    longest = max(segs, key=len)
    return np.asarray(longest)   # shape (N, 2): columns [R, Z]


def read_xpoint_from_hdf5(hdf5_path, null=None):
    """
    Read X-point initial guess from an M3DC1 C1.h5 file.

    M3DC1 stores X-point coordinates in /scalars/ as float32 datasets with
    shape (N,) where all elements repeat the same value (multiple time steps).
    Single-null configurations zero out xnull2/znull2 as a sentinel.

    Parameters
    ----------
    hdf5_path : str or Path
        Path to the M3DC1 HDF5 file (typically C1.h5).
    null : int or None
        1 — primary null (scalars/xnull, scalars/znull)
        2 — secondary null (scalars/xnull2, scalars/znull2)
        None — auto: selects the primary null for single-null configurations;
               raises ValueError for double-null configurations.

    Returns
    -------
    R, Z : float
        X-point coordinates in metres.

    Raises
    ------
    ValueError
        If null=2 is requested but no secondary X-point exists, or if
        null=None is used with a double-null configuration.
    """
    import h5py
    with h5py.File(hdf5_path, 'r') as f:
        R1 = float(f['scalars/xnull'][0])
        Z1 = float(f['scalars/znull'][0])
        has_secondary = (
            'scalars/xnull2' in f and
            'scalars/znull2' in f and
            float(f['scalars/xnull2'][0]) != 0.0
        )
        if has_secondary:
            R2 = float(f['scalars/xnull2'][0])
            Z2 = float(f['scalars/znull2'][0])

    if null == 1:
        return R1, Z1
    if null == 2:
        if not has_secondary:
            raise ValueError(
                "xpoint null=2 requested but no secondary X-point found "
                "(single-null configuration or xnull2/znull2 absent)."
            )
        return R2, Z2
    # auto
    if has_secondary:
        raise ValueError(
            f"Double-null configuration: two X-points found. "
            f"Specify null=1 (R={R1:.6f}, Z={Z1:.6f}) or "
            f"null=2 (R={R2:.6f}, Z={Z2:.6f})."
        )
    return R1, Z1


def trace_separatrix(tracer, phi, r_xpoint, z_xpoint, stability=0,
                     n_intervals=9, n_segments=30,
                     l_lim=0.005, theta_lim=20.0,
                     epsilon=1e-8, h=1e-8, tol=1e-14, max_iter=50,
                     precision_limit=1e-14, max_insertions=50,
                     verbose=False):
    """
    Trace the equilibrium separatrix by iterating the Poincaré map.

    In the unperturbed equilibrium field (no MPs), the stable and unstable
    manifolds of the X-point coincide with the separatrix. This function
    traces the manifold segment by segment and closes the resulting open
    curve by appending the X-point as the final row.

    The tracer must already be configured (step size, wall monitor) and
    loaded with the equilibrium field (timeslice=-1 for M3DC1Source).

    Parameters
    ----------
    tracer : Maglit
        Pre-configured field-line integrator.
    phi : float
        Toroidal angle of the Poincaré section (radians).
    r_xpoint, z_xpoint : float
        Initial guess for the X-point position (metres).
    stability : int
        0 = stable manifold (forward map), 1 = unstable manifold.
        Both coincide with the separatrix in the equilibrium field.
    n_intervals : int
        Intervals for the primary segment; produces n_intervals+1 points.
    n_segments : int
        Total number of segments including the primary. Adjust until the
        curve visually closes in a plot. Typical range: 20–80.
    l_lim : float
        Arc-length refinement threshold (metres).
    theta_lim : float
        Turning-angle refinement threshold (degrees).
    epsilon : float
        Distance from the X-point to the pivot point.
    h : float
        Step size for the numerical Jacobian.
    tol : float
        Newton convergence tolerance.
    max_iter : int
        Maximum Newton iterations.
    precision_limit : float
        Minimum arc length below which point insertion is skipped.
    max_insertions : int
        Maximum insertions per refinement pass.
    verbose : bool
        Print per-segment progress when True.

    Returns
    -------
    closed : ndarray of shape (N, 2)
        Closed separatrix curve [R, Z] in metres. The X-point is appended
        as the final row to close the loop.
    x_point : ndarray of shape (2,)
        [R, Z] of the converged X-point.

    Raises
    ------
    RuntimeError
        If the X-point Newton iteration does not converge.
    """
    from ._maglib import Manifold

    mf = Manifold(tracer, phi=phi, stability=stability)
    mf.configure(epsilon=epsilon, h=h, tol=tol, max_iter=max_iter,
                 precision_limit=precision_limit, max_insertions=max_insertions)

    if not mf.find_x_point(r_xpoint, z_xpoint):
        raise RuntimeError(
            f"X-point Newton iteration did not converge "
            f"(initial guess R={r_xpoint}, Z={z_xpoint})."
        )
    x_point = mf.x_point

    seg = mf.primary_segment(n_intervals)
    if verbose:
        print(f"  Segment 1 (primary): {seg.shape[0]} pts")

    for i in range(1, n_segments):
        _, seg = mf.new_segment(seg, l_lim, theta_lim)
        if verbose:
            print(f"  Segment {i+1:3d}: {seg.shape[0]:4d} pts")

    all_pts = np.vstack(mf.output_data)
    closed = np.vstack([all_pts, x_point])
    return closed, x_point


def _extract_segments(source):
    if hasattr(source, 'output_data'):
        return source.output_data
    if isinstance(source, list) and all(isinstance(s, np.ndarray) for s in source):
        return source
    raise TypeError("source must be a Manifold object or a list of (N, 2) numpy arrays.")
