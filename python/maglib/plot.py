import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter


def plot_footprint(source, which_plot="all", xaxis="rad", cmap="jet", cmap_key=10,
                   figsize=(10, 5), sizef=1, dpi=80, v_min=None, v_max=None,
                   turn_cap=None, psi_cap=False, savefig=None):
    """
    Plot footprint data from fpgen output or a maglib.Footprint object.

    Parameters
    ----------
    source : str, numpy.ndarray, or Footprint
        File path, (N, 6) array, or Footprint object.
        Columns: R0, Z0, phi0, connection_length, psiN_min, turns.
    which_plot : str
        Subset to plot: "cl", "psi", "turns", "au", or "all". Default "all".
    xaxis : str
        X-axis units: "rad" or "deg". Default "rad".
    cmap : str
        Matplotlib colormap name. Default "jet".
    cmap_key : int
        Number of discrete colors used for "turns" and "au" plots. Default 10.
    figsize : tuple of float
        Base figure size in inches. Default (10, 5).
    sizef : float
        Uniform scale factor applied to figsize. Default 1.
    dpi : int
        Figure resolution. Default 80.
    v_min, v_max : float, optional
        Manual colormap bounds for the "au" plot. Both must be provided together.
    turn_cap : tuple of float, optional
        (min, max) color scale bounds for the "turns" plot. Uses LogNorm when None.
    psi_cap : bool
        Mask points where psiN_min >= 1 or turns <= 1 (open field lines). Default False.
    savefig : str, optional
        Base path for saved figures. Files are written as <savefig>_cl.png, etc.
    """
    if isinstance(source, str):
        data = np.loadtxt(source)
    elif isinstance(source, np.ndarray):
        data = source
    elif hasattr(source, 'output_data'):
        data = source.output_data
    else:
        raise TypeError("source must be a file path, numpy array, or Footprint object.")

    if data.ndim != 2 or data.shape[1] < 6:
        raise ValueError("Footprint data must have shape (N, 6).")

    if which_plot not in ("cl", "psi", "turns", "au", "all"):
        raise ValueError("which_plot must be 'cl', 'psi', 'turns', 'au', or 'all'.")
    if xaxis not in ("rad", "deg"):
        raise ValueError("xaxis must be 'rad' or 'deg'.")

    # Detect plate orientation from uniqueness of R column
    if np.unique(data[:, 0]).size == 1:
        ylabel = "$Z$ ( m )"
        y_col = 1
    else:
        ylabel = "$R$ ( m )"
        y_col = 0

    x = data[:, 2]
    y = data[:, y_col]
    sort_idx = np.lexsort((y, x))
    data = data[sort_idx]
    x = data[:, 2]
    y = data[:, y_col]

    z_cl = data[:, 3]
    z_psi = data[:, 4]
    z_turns = data[:, 5]

    n_phi = len(np.unique(x))
    n_y = len(data) // n_phi
    z_cl = np.reshape(z_cl, (n_phi, n_y)).T
    z_psi = np.reshape(z_psi, (n_phi, n_y)).T
    z_turns = np.reshape(z_turns, (n_phi, n_y)).T

    if psi_cap:
        mask = (z_psi >= 1) | (z_turns <= 1)
        z_cl = np.ma.masked_where(mask, z_cl)
        z_psi = np.ma.masked_where(mask, z_psi)
        z_turns = np.ma.masked_where(mask, z_turns)

    cmap_cont = mpl.colormaps.get_cmap(cmap)
    cmap_disc = mpl.colormaps.get_cmap(cmap).resampled(cmap_key)
    for cm in (cmap_cont, cmap_disc):
        cm.set_under("white")
        if psi_cap:
            cm.set_bad("white")

    fsize = (figsize[0] * sizef, figsize[1] * sizef)
    extent = [0, 360, float(np.min(y)), float(np.max(y))]
    imshow_kw = dict(origin='lower', aspect='auto', extent=extent)

    def _xaxis(ax):
        if xaxis == "rad":
            def _fmt(val, pos):
                r = val * np.pi / 180
                return f"{r:.1f}" if r > 0 else "0"
            ax.xaxis.set_major_formatter(FuncFormatter(_fmt))
            ax.set_xlabel(r"$\phi$ ( rad )")
        else:
            ax.set_xlabel(r"$\phi$ ( deg )")
        ax.set_ylabel(ylabel)

    def _save(suffix):
        if savefig is not None:
            plt.savefig(f"{savefig}_{suffix}.png", dpi=dpi, bbox_inches='tight')

    if which_plot in ("all", "cl"):
        fig, ax = plt.subplots(figsize=fsize, dpi=dpi)
        im = ax.imshow(z_cl, cmap=cmap_cont, norm=LogNorm(), **imshow_kw)
        fig.colorbar(im, ax=ax).set_label("connection length ( m )")
        _xaxis(ax)
        _save("cl")
        plt.show(block=False)

    if which_plot in ("all", "psi"):
        fig, ax = plt.subplots(figsize=fsize, dpi=dpi)
        im = ax.imshow(z_psi, cmap=cmap_cont, **imshow_kw)
        fig.colorbar(im, ax=ax).set_label(r"$\psi_{N\,\,\mathrm{min}}$")
        _xaxis(ax)
        _save("psi")
        plt.show(block=False)

    if which_plot in ("all", "turns"):
        norm = None if turn_cap is not None else LogNorm()
        vmin = turn_cap[0] if turn_cap is not None else None
        vmax = turn_cap[1] if turn_cap is not None else None
        fig, ax = plt.subplots(figsize=fsize, dpi=dpi)
        im = ax.imshow(z_turns, cmap=cmap_disc, norm=norm, vmin=vmin, vmax=vmax, **imshow_kw)
        fig.colorbar(im, ax=ax).set_label("toroidal turns")
        _xaxis(ax)
        _save("turns")
        plt.show(block=False)

    if which_plot in ("all", "au"):
        if v_min is not None and v_max is not None:
            if v_min >= v_max:
                raise ValueError("v_min must be less than v_max.")
            z_norm = ((1.0 - 1.0 / (z_psi * z_cl)) - v_min) / (v_max - v_min)
        else:
            z_raw = 1.0 / (z_cl * z_psi)
            lo, hi = float(np.nanmin(z_raw)), float(np.nanmax(z_raw))
            if lo == hi:
                raise ValueError("All 1/(psiN_min * CL) values are identical; cannot normalize.")
            z_norm = (z_raw - lo) / (hi - lo)
        fig, ax = plt.subplots(figsize=fsize, dpi=dpi)
        im = ax.imshow(z_norm, cmap=cmap_disc, **imshow_kw)
        fig.colorbar(im, ax=ax).set_label(
            r"$(\psi_{N\,\,\mathrm{min}} \cdot \mathrm{CL})^{-1}_\mathrm{norm}$")
        _xaxis(ax)
        _save("au")
        plt.show(block=False)

    plt.show()


def plot_manifold(source, cuts=None, labels=None, wall=None, linestyle="line",
                  figsize=(5, 8), dpi=80, sizef=1, linewidth=2, scatter_size=5,
                  xlim=None, zlim=None, savefig=None):
    """
    Plot manifold curves from mfgen output or a maglib.Manifold object.

    Parameters
    ----------
    source : str, list of str, list of ndarray, or Manifold
        A single file path, list of file paths, list of (N, 2) arrays, or a
        Manifold object. For a Manifold object all segments share one color and
        label. For a list each entry gets its own color and label.
    cuts : list of int, optional
        Point limit per curve (file/array inputs only). Use 0 or None for all
        points; negative values slice from the end.
    labels : list of str, optional
        Legend labels, one per curve. Defaults to "Manifold N" for list inputs
        or "Manifold" for a Manifold object.
    wall : str, optional
        Path to a two-column (R, Z) wall contour file.
    linestyle : str
        "line", "scatter", or "both". Default "line".
    figsize : tuple of float
        Base figure size in inches. Default (5, 8).
    dpi : int
        Figure resolution. Default 80.
    sizef : float
        Uniform scale factor applied to figsize. Default 1.
    linewidth : float
        Line width. Default 2.
    scatter_size : float
        Marker size for scatter points. Default 5.
    xlim : tuple of float, optional
        R-axis limits.
    zlim : tuple of float, optional
        Z-axis limits.
    savefig : str, optional
        Path to save the figure.
    """
    if linestyle not in ("line", "scatter", "both"):
        raise ValueError("linestyle must be 'line', 'scatter', or 'both'.")

    from_manifold_object = hasattr(source, 'output_data')

    if isinstance(source, str):
        _cuts = cuts or [None]
        segments = [_load_file(source, _cuts[0])]
        labels = labels or ["Manifold"]
    elif isinstance(source, list) and len(source) > 0 and isinstance(source[0], str):
        _cuts = cuts or [None] * len(source)
        if len(_cuts) != len(source):
            raise ValueError("cuts length must match source length.")
        segments = [_load_file(f, c) for f, c in zip(source, _cuts)]
        labels = labels or [f"Manifold {i+1}" for i in range(len(source))]
    elif isinstance(source, list) and len(source) > 0 and isinstance(source[0], np.ndarray):
        _cuts = cuts or [None] * len(source)
        if len(_cuts) != len(source):
            raise ValueError("cuts length must match source length.")
        segments = [_apply_cut(s, c) for s, c in zip(source, _cuts)]
        labels = labels or [f"Manifold {i+1}" for i in range(len(source))]
    elif from_manifold_object:
        segments = source.output_data
        labels = labels or ["Manifold"]
    elif isinstance(source, list) and len(source) == 0:
        raise ValueError("source list is empty.")
    else:
        raise TypeError(
            "source must be a file path, list of file paths, list of arrays, or Manifold object.")

    fig, ax = plt.subplots(figsize=(figsize[0] * sizef, figsize[1] * sizef), dpi=dpi)
    ax.set_aspect('equal', adjustable='box')

    if from_manifold_object:
        color = None
        for i, seg in enumerate(segments):
            lbl = labels[0] if i == 0 else "_nolegend_"
            color = _draw_segment(ax, seg, lbl, linestyle, linewidth, scatter_size, color)
    else:
        for seg, lbl in zip(segments, labels):
            _draw_segment(ax, seg, lbl, linestyle, linewidth, scatter_size)

    if wall is not None:
        wall_data = np.loadtxt(wall, ndmin=2)
        ax.plot(wall_data[:, 0], wall_data[:, 1], color="black",
                linewidth=linewidth, label="First wall")

    if xlim is not None:
        ax.set_xlim(xlim)
    if zlim is not None:
        ax.set_ylim(zlim)

    ax.set_xlabel("R ( m )", fontsize=12)
    ax.set_ylabel("Z ( m )", fontsize=12)
    ax.legend()
    plt.tight_layout()

    if savefig is not None:
        plt.savefig(savefig, dpi=dpi, bbox_inches='tight')

    plt.show()


def _load_file(path, cut):
    return _apply_cut(np.loadtxt(path, ndmin=2), cut)


def _apply_cut(arr, cut):
    if cut is None or cut == 0:
        return arr
    return arr[:cut]


def _draw_segment(ax, seg, label, linestyle, linewidth, scatter_size, color=None):
    """Draw one segment; returns the color chosen by matplotlib."""
    R, Z = seg[:, 0], seg[:, 1]
    kw = {"color": color} if color is not None else {}

    if linestyle == "line":
        line, = ax.plot(R, Z, label=label, linewidth=linewidth, **kw)
        return line.get_color()
    elif linestyle == "scatter":
        sc = ax.scatter(R, Z, label=label, s=scatter_size, **kw)
        return sc.get_facecolor()[0]
    else:  # "both"
        line, = ax.plot(R, Z, label=label, linewidth=linewidth, **kw)
        used = line.get_color()
        ax.scatter(R, Z, s=scatter_size, color=used)
        return used
