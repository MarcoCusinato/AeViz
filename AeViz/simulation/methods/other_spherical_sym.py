from AeViz.simulation.methods import *
from AeViz.utils.profiles import calculate_profile
from AeViz.utils.tidal_love import solve_tidal_love_profile
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
def tidal_deformability(self, tob_corrected=True, save_checkpoints=True,
                        comp:Literal['PNS_core', 'PNS']='PNS', **kwargs):
    """
    Returns the tidal deformability at every timestep.
    If tob_corrected is True, the time is corrected for the time of
    bounce. If save_checkpoints is True, the checkpoints are saved
    during the calculation.
    Returns: time, tidal deformability
    """
    if comp == 'PNS_core':
        time, _, data = solve_tidal_love_profile(self, save_checkpoints)
    else:
        time, data, _ = solve_tidal_love_profile(self, save_checkpoints)
    if not tob_corrected:
        time += self.tob
    return [time, data['lambda']]

@smooth
@derive
def love_number(self, tob_corrected=True, save_checkpoints=True,
                comp:Literal['PNS_core', 'PNS']='PNS', **kwargs):
    """
    Returns the Love number at every timestep.
    If tob_corrected is True, the time is corrected for the time of
    bounce. If save_checkpoints is True, the checkpoints are saved
    during the calculation.
    Returns: time, Love number
    """
    if comp == 'PNS_core':
        time, _, data = solve_tidal_love_profile(self, save_checkpoints)
    else:
        time, data, _ = solve_tidal_love_profile(self, save_checkpoints)
    if not tob_corrected:
        time += self.tob
    return [time, data['kappa2']]


## -----------------------------------------------------------------
## PROFILES
## -----------------------------------------------------------------

def radial_profile(self, quantity, save_checkpoints=True, **kwargs):
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