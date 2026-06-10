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
    ndarray of shape (N_lobes, 6)
        Columns: Rmid, Zmid, angle [rad], perimeter [m], area [m²], h_param [m].
        Returns shape (0, 6) when fewer than two intersections are found.
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
        return np.empty((0, 6))

    gtype = inter_geom.geom_type
    if gtype == 'Point':
        pts = [inter_geom]
    elif gtype == 'MultiPoint':
        pts = list(inter_geom.geoms)
    else:
        pts = [g for g in inter_geom.geoms if g.geom_type == 'Point']

    if len(pts) < 2:
        return np.empty((0, 6))

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

        # Midpoint: centroid of equilibrium boundary (all points except p2)
        mid = eq_coords[:-1].mean(axis=0) if len(eq_coords) > 1 else eq_coords[0].copy()

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

        rows.append([mid[0], mid[1], angle, perimeter, area, h_param])

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


def _extract_segments(source):
    if hasattr(source, 'output_data'):
        return source.output_data
    if isinstance(source, list) and all(isinstance(s, np.ndarray) for s in source):
        return source
    raise TypeError("source must be a Manifold object or a list of (N, 2) numpy arrays.")
