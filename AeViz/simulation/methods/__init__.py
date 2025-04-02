from AeViz.units import u
from AeViz.units.aeseries import aerray, aeseries
from AeViz.utils.decorators.simulation import (smooth, derive, hdf_isopen,
                                    subtract_tob, sum_tob)
from AeViz.utils.decorators.grid import get_grid
import numpy as np
from typing import Literal