from AeViz.simulation.methods import *
from AeViz.utils.physics.profiles import calculate_profile
from AeViz.utils.physics.tidal_love import solve_tidal_love_profile
from AeViz.utils.decorators.grid import profile_grid
import re
from typing import Literal

"""
Functions to calculate tidal deformation and radial profiles from a
simulation in spherical symmetry.
These functions are not meant to be used standalone, but rather to be
imported into the Simulation class.
"""

## -----------------------------------------------------------------
## TIDAL DATA
## -----------------------------------------------------------------

@smooth
@derive
@sum_tob
def tidal_deformability(self, tob_corrected=True, save_checkpoints=True,
                        comp:Literal['PNS_core', 'PNS']='PNS', **kwargs):
    """
    Returns the tidal deformability at every timestep.
    If tob_corrected is True, the time is corrected for the time of
    bounce. If save_checkpoints is True, the checkpoints are saved
    during the calculation.
    Returns: time, tidal deformability
    """
    lambda_pns, _, lambda_core, _ = solve_tidal_love_profile(self,
                                                             save_checkpoints)
    if comp == 'PNS_core':
        return lambda_core
    else:
        return lambda_pns

@smooth
@derive
@sum_tob
def love_number(self, tob_corrected=True, save_checkpoints=True,
                comp:Literal['PNS_core', 'PNS']='PNS', **kwargs):
    """
    Returns the Love number at every timestep.
    If tob_corrected is True, the time is corrected for the time of
    bounce. If save_checkpoints is True, the checkpoints are saved
    during the calculation.
    Returns: time, Love number
    """
    kappa_pns, _, kappa_core, _ = solve_tidal_love_profile(self,
                                                             save_checkpoints)
    if comp == 'PNS_core':
        return kappa_core
    else:
        return kappa_pns

## -----------------------------------------------------------------
## PROFILES
## -----------------------------------------------------------------

@profile_grid
@derive
@sum_tob
def radial_profile(self, quantity, tob_corrected=True,
                   save_checkpoints=True, **kwargs):
    """
    Calucaltes the radial profile of the selected quantity. Name of
    the quantity must be the same as the one of the simulation
    methods.
    Returns: time, radius, profile
    If the array contains polarization (i.e. the magnetic fields)
    besides the regular angular and radial dependence, this will be
    returned as the last axis.
    """
    return calculate_profile(self, quantity, save_checkpoints, **kwargs)

def mass_grid(self, file_name, **kwargs):
    """
    Calculates the grid in mass at a specific timestep
    """
    dm = np.sum(self.mass(file_name), axis=tuple(range(self.dim - 1)))
    mgrid = np.cumsum(dm.to(u.M_sun))
    mgrid.set(label=r"$M$", name="mass_grid", log=False, cmap=None,
              limits=[0, mgrid.max().value])
    return mgrid