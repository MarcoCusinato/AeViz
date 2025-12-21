from AeViz.units.aerray import aerray
from AeViz.units import u
from AeViz.units.constants import constants as c
import numpy as np

def E_nu_left(self):
    """
    Method that returns an array with the left neutrino energy grid.
    Results:
        E_nu_left: (numpy array)
    """
    if self.path_grid is None:
        raise TypeError("You shouldn't be here.")
    return aerray(self._cell__nu_grid_file[:, 0], u.MeV, name='E_left',
                  label=r'$E_\mathrm{l}$')

def E_nu_right(self):
    """
    Method that returns an array with the right neutrino energy grid.
    Results:
        E_nu_right: (numpy array)
    """
    if self.path_grid is None:
        raise TypeError("You shouldn't be here.")
    return aerray(self._cell__nu_grid_file[:, 2], u.MeV, name='E_right',
                  label=r'$E_\mathrm{r}$')

def E_nu(self):
    """
    Method that returns an array with the neutrino energy grid.
    Results:
        E_nu: (numpy array)
    """
    if self.path_grid is None:
        raise TypeError("You shouldn't be here.")
    return aerray(self._cell__nu_grid_file[:, 1], u.MeV, name='E',
                  label=r'$E$')

def dE_nu(self):
    """
    Method that returns an array with the neutrino energy integration element.
    Results:
        dE_nu: (numpy array)
    """
    dE = self.E_nu_right() - self.E_nu_left()
    dE.set(name='dE', label=r'$\mathrm{d}E$')
    return dE

    
def e(self):
    """
    Computes the small e for the neutrino fluxes
    """
    return (c.h  * c.c / (self.E_nu() * 4 * np.pi)) ** 3