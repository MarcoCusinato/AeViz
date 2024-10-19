
"""
This module contains functions that are used to find and set the limits
of the plots. This is a utility module.

"""
import numpy as np

## GW limits
def GW_limit(GW_data):
    """
    Retursn the limits for the GW data.
    Usually the limits are set to the maximum absolute value of the data,
    however if crashes are to happen, the limits are set to 150.
    """   
    lim = np.max(np.abs(GW_data))
    lim *= 1.1
    if lim > 150 or np.isnan(lim):
            lim = 150
    return (-lim, lim)

def set2Dlims(ax, xlim, ylim, number, form_factor, sim_dim):
    """
    We set the limits given the axes and the number and positioning of
    the figure which is part of.
    Used for 2D plots ie contourf plots.
    """
    if number == 5 and form_factor == 1:
        ## This is the GWs decomposition
        if xlim == None:
            ax["E"].set_ylim(ylim)
        if ylim == None:
            ax["E"].set_xlim(xlim)
    else:
        if sim_dim == 2:
            ## Used for 2D simulations, we do not extend the grid
            set2dlims2Dsim(ax, xlim, ylim, number, form_factor)
        else:
            ## Used for 3D simulations
            set2Dlims3Dsim(ax, xlim, ylim, number, form_factor)
                
def set2dlims2Dsim(ax, xlim, ylim, number, form_factor):
    """
    Set the limits for contourf plots of 2D simulations.
    The aspect is set to equal to maintain the aspect ratio, and they share the
    same limits.
    Up to 2 plots the vertical limits are doouble the horizontal limits.
    From 3 plots and on, the limits are the same.
    """
    if xlim == None:
        ylim = list(ylim)
        ylim.sort()
        xlim = [0]
        if number == 4 and form_factor == 2:
            ylim[0] = 0
        else:
            ylim[0] = -ylim[1]
        xlim.append(ylim[1])
        
    if ylim == None:
        xlim = list(xlim)
        xlim.sort()
        xlim[0] = 0
        if number == 4 and form_factor == 2:
            ylim = xlim
        else:
            ylim = [-xlim[1], xlim[1]]

    if number == 1 and form_factor == 2:
        ax["A"].set(xlim=(xlim[0], xlim[1]), ylim=(ylim[0], ylim[1]), aspect=1)
    elif number == 2 and form_factor == 2:
        ax["A"].set(xlim=(xlim[1], xlim[0]), ylim=(ylim[0], ylim[1]), aspect=1)
        ax["B"].set(xlim=(xlim[0], xlim[1]), ylim=(ylim[0], ylim[1]), aspect=1)
    elif number == 3 and form_factor == 2:
        ax["B"].set(xlim=(xlim[1], xlim[0]), ylim=(ylim[0], ylim[1]), aspect=1)
        ax["C"].set(xlim=(xlim[0], xlim[1]), ylim=(ylim[0], ylim[1]), aspect=1)
    elif number == 4 and form_factor == 2:
        ax["A"].set(xlim=(xlim[1], xlim[0]), ylim=(ylim[0], ylim[1]), aspect=1)
        ax["B"].set(xlim=(xlim[1], xlim[0]), ylim=(-ylim[1], ylim[0]),
                    aspect=1)
        ax["C"].set(xlim=(xlim[0], xlim[1]), ylim=(ylim[0], ylim[1]), aspect=1)
        ax["D"].set(xlim=(xlim[0], xlim[1]), ylim=(-ylim[1], ylim[0]),
                    aspect=1)

def set2Dlims3Dsim(ax, xlim, ylim, number, form_factor):
    """
    Set the limits for contourf plots of 3D simulations.
    The aspect is set to equal to maintain the aspect ratio, and they share the
    same limits.
    In case of 3D simulations and one plot, an entire slice is plotted,
    so the limits are equal to preserve the square shape of the plot.
    From 2 plots and on, the limits are the same as in the 2D case.
    """
    if xlim == None:
        ylim = list(ylim)
        ylim.sort()
        ylim[0] = -ylim[1]
        xlim = ylim
    if ylim == None:
        xlim = list(xlim)
        xlim.sort()
        xlim[0] = -xlim[1]
        ylim = xlim
            
    if number == 1 and form_factor == 2:
        ax["A"].set(xlim=(xlim[0], xlim[1]), ylim=(ylim[0], ylim[1]), aspect=1)
    elif number == 2 and form_factor == 2:
        ax["A"].set(xlim=(xlim[0], 0), ylim=(ylim[0], ylim[1]), aspect=1)
        ax["B"].set(xlim=(0, xlim[1]), ylim=(ylim[0], ylim[1]), aspect=1)
    elif number == 3 and form_factor == 2:
        ax["B"].set(xlim=(xlim[0], 0), ylim=(ylim[0], ylim[1]), aspect=1)
        ax["C"].set(xlim=(0, xlim[1]), ylim=(ylim[0], ylim[1]), aspect=1)
    elif number == 4 and form_factor == 2:
        ax["A"].set(xlim=(xlim[0], 0), ylim=(0, ylim[1]), aspect=1)
        ax["B"].set(xlim=(xlim[0], 0), ylim=(ylim[0], 0), aspect=1)
        ax["C"].set(xlim=(0, xlim[1]), ylim=(0, ylim[1]), aspect=1)
        ax["D"].set(xlim=(0, xlim[1]), ylim=(ylim[0], 0), aspect=1)