from __future__ import annotations
from AeViz.simulation import Simulation
from AeViz.units import u, aerray
import numpy as np
from AeViz.utils.physics.radii_utils import shock_radius
from AeViz.utils.math_utils import function_average_radii
from AeViz.units.constants import constants as c

def declare_shock_dictionary(dim: int,
                           magdim: int) -> dict:
    """
    _summary_

    Parameters
    ----------
    dim : int
        _description_
    magdim : int
        _description_

    Returns
    -------
    dict
        _description_
    """