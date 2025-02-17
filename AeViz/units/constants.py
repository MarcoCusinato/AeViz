from AeViz.units.aerray import aerray
from astropy import units as u
import numpy as np

class constants:
    c = aerray(2.99792458e10, u.cm/u.s, 'speed_of_light', label='$c$')
    G = aerray(6.6743e-8, u.cm**3/(u.g * u.s**2), 'gravitational_constant',
               label='$G$')
    h = aerray(6.62607015e-27, u.erg * u.s, 'planck_constant', label='$h$')
    hbar = aerray(1.0545718e-27, u.erg * u.s, 'reduced_planck_constant',
                  label='$\hbar$')
    