from ._maglib import Sode, SodeMethod, SodeStatus, Collider, FieldSource, M3DC1Source, Maglit, Footprint, Manifold
from .plot import plot_footprint, plot_manifold
from .manifold_tools import simplify, to_linestrings, lobe_map, get_separatrix, trace_separatrix

__all__ = [
    "Sode", "SodeMethod", "SodeStatus",
    "Collider",
    "FieldSource", "M3DC1Source",
    "Maglit",
    "Footprint", "plot_footprint",
    "Manifold", "plot_manifold", "simplify", "to_linestrings", "lobe_map",
    "get_separatrix", "trace_separatrix",
]
