# Import required modules
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter


def plot_fp(filename=None, plate="h", which_plot="all", xaxis="rad", cmap="jet", cmap_key=10, figsize=(10, 5), sizef=1, dpi=80, norm_f=6, v_min=None, v_max=None):
    """
    Function to plot footprints evaluated by the maglib code fpgen.
    
    Parameters
    ----------
    filename : str
        Path to the file containing the data.
    plate : str
        Divertor plate: "h" for horizontal (floor) or "v" for vertical (wall). Default is "h".
    which_plot : str
        Type of plot to generate: "cl" (connection length), "psi" (psi min), or "au" (arbitrary units), "all". Default is "all".
    xaxis : str
        Type of x-axis: "rad" for radians or "deg" for degrees. Default is "rad".
    cmap : str
        Colormap to use for the plots. Default is "jet".
    cmap_key : int
        Number of colors in the discretized version of the colormap. Default is 10.
    figsize : tuple
        Size of the figure in inches. Default is (10, 5).
    sizef : float
        Size factor for the figure (increase or decrease figure size without changig aspect ratio). Default is 1.
    dpi : int
        Dots per inch for the figure. Default is 80 (default is low for better performance).
    norm_f : float
        Normalization factor for the arbitrary units plot. Default is 6.
    v_min : float
        Minimum value for the normalized colormap in arbitrary units. Default is None.
    v_max : float
        Maximum value for the normalized colormap in arbitrary units. Default is None.

    Returns
    -------
    None
        Displays the plots.

    """


    # Load data
    if filename is None:
        raise ValueError("Filename cannot be None. Please provide a valid file path.")
    data = np.loadtxt(filename)

    # Extract relevant columns
    x = data[:, 2] #phi
    if plate == "h":
        y = data[:, 0] #r
    elif plate == "v":
        y = data[:, 1] #z
    else:
        raise ValueError("Invalid plate type. Use 'h' for horizontal or 'v' for vertical.")

    z_cl = data[:, 3] #connection length
    z_psi = data[:, 4] #psin

    # Determine the grid size
    num_rows = int(len(np.unique(x))) # Number of rows
    num_columns = int(len(data)/num_rows) # Number of columns
    z_cl = np.reshape(z_cl, (num_rows, num_columns)).transpose() # z column to matrix
    z_psi = np.reshape(z_psi, (num_rows, num_columns)).transpose() # z column to matrix

    # Custom x-axis formatter: convert degrees to radians
    def degrees_to_radians(x, pos):
        radians = x * np.pi / 180  # Convert degrees to radians
        return f"{radians:.1f}" if radians > 0 else "0"

    # Custom colormap
    cmap_cap = mpl.cm.get_cmap(cmap).copy()
    cmap_f = mpl.cm.get_cmap(cmap,cmap_key).copy()

    # Modify colormap to set under values to white
    cmap_cap.set_under("white")
    cmap_f.set_under("white")


    #### Connection Length plot ####
    if which_plot == "all" or which_plot == "cl":
        fsize=(figsize[0]*sizef, figsize[1]*sizef)
        plt.figure(figsize=fsize, dpi=dpi)
        plt.imshow(z_cl, cmap=cmap, extent=[0, 360, np.min(y), np.max(y)], origin='lower', aspect='auto', norm=LogNorm())
        plt.colorbar().set_label("connection length ( m )")
        if xaxis == "rad":
            plt.gca().xaxis.set_major_formatter(FuncFormatter(degrees_to_radians))
            plt.xlabel("$\phi$ ( rad )")
        elif xaxis == "deg":
            plt.xlabel("$\phi$ ( deg )")
            pass
        else:
            raise ValueError("Invalid x-axis type. Use 'rad' for radians or 'deg' for degrees.")
        plt.ylabel("$R$ ( m )")

        # Show plot
        plt.show(block=False)


    #### Psi min plot ####
    if which_plot == "all" or which_plot == "psi":
        fsize=(figsize[0]*sizef, figsize[1]*sizef)
        plt.figure(figsize=fsize, dpi=dpi)
        plt.imshow(z_psi, cmap=cmap, extent=[0, 360, np.min(y), np.max(y)], origin='lower', aspect='auto')
        plt.colorbar().set_label(r'$\psi_{n\,\,\mathrm{min}}$')
        if xaxis == "rad":
            plt.gca().xaxis.set_major_formatter(FuncFormatter(degrees_to_radians))
            plt.xlabel("$\phi$ ( rad )")
        elif xaxis == "deg":
            plt.xlabel("$\phi$ ( deg )")
            pass
        else:
            raise ValueError("Invalid x-axis type. Use 'rad' for radians or 'deg' for degrees.")
        plt.ylabel("$R$ ( m )")

        # Show plot
        plt.show(block=False)


    #### arbitrary units plot ####
    if which_plot == "all" or which_plot == "au":

        # Check if v_min and v_max are provided
        if v_min is not None and v_max is not None:
            if not isinstance(v_min, (int, float)) or not isinstance(v_max, (int, float)):
                raise ValueError("v_min and v_max must be numbers.")
            if v_min >= v_max:
                raise ValueError("v_min must be less than v_max.")
            # Normalize z_norm to the range [v_min, v_max]
            z_norm = ((1 - (1/(z_psi*z_cl))) - v_min) / (v_max - v_min)
        else:
            # default normalization
            z_norm = (1-(1/(z_psi*z_cl)-np.min(1/(z_psi*z_cl)))*norm_f/np.max(1/(z_psi*z_cl)))
        
        fsize=(figsize[0]*sizef, figsize[1]*sizef)
        plt.figure(figsize=fsize, dpi=dpi)
        plt.imshow(z_norm, cmap=cmap_f, extent=[0, 360, np.min(y), np.max(y)], origin='lower', aspect='auto', vmin=0, vmax=1)
        plt.colorbar().set_label('$( \psi_{n\,\,\mathrm{min}}$ . connection length ) $^{-1}_\mathrm{norm}$') 
        if xaxis == "rad":
            plt.gca().xaxis.set_major_formatter(FuncFormatter(degrees_to_radians))
            plt.xlabel("$\phi$ ( rad )")
        elif xaxis == "deg":
            plt.xlabel("$\phi$ ( deg )")
            pass
        else:
            raise ValueError("Invalid x-axis type. Use 'rad' for radians or 'deg' for degrees.")
        plt.ylabel("$R$ ( m )")

        # Show plot
        plt.show(block=False)