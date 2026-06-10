import numpy as np

try:
    from shapely.geometry import LineString
    _HAS_SHAPELY = True
except ImportError:
    _HAS_SHAPELY = False

_SHAPELY_MSG = "shapely is required. Install with: pip install shapely"


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
