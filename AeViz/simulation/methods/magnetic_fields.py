from AeViz.simulation.methods import *
from AeViz.utils.math_utils import strfct2D
from AeViz.units.constants import constants as c
from typing import Literal

"""
Functions to handle magnetic fields data from a simulation.
These functions are not meant to be used standalone, but rather to be
imported into the Simulation class.
"""

## -----------------------------------------------------------------
## MAGNETIC FIELDS DATA
## -----------------------------------------------------------------
@hdf_isopen
def __CT_magnetic_fields(self, file_name, **kwargs):
    """
    Magnetic field at the cells border, use this ONLY to calculate
    streamlines. If you want to plot the actual magnetic fields use
    the 'magnetic_field' method.
    """
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['mag_CT/data'][...]), self.dim)
    return aerray(data, u.G, 'cell_magnetic_fields')

@get_grid
@smooth
@hdf_isopen
def magnetic_fields(self, file_name, comp:Literal['all', 'r', 'th', 'ph']='all',
                    **kwargs):
    """
    Magnetic field at the cells center.
    """
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['mag_vol/data'][...]), self.dim)
    if comp == 'all':
        return (aerray(data[..., 0], u.G, 'B_r', r'$B_r$', 'coolwarm', [-1e15, 1e15],
                    True), 
                aerray(data[..., 1], u.G, 'B_theta', r'$B_\theta$', 'coolwarm',
                    [-1e15, 1e15], True),
                aerray(data[..., 2], u.G, 'B_phi', r'$B_\phi$', 'coolwarm',
                    [-1e15, 1e15], True))
    elif comp == 'r':
        return aerray(data[..., 0], u.G, 'B_r', r'$B_r$', 'coolwarm', [-1e15, 1e15],
                    True)
    elif comp == 'th':
        return aerray(data[..., 1], u.G, 'B_theta', r'$B_\theta$', 'coolwarm',
                    [-1e15, 1e15], True)
    elif comp == 'ph':
        return aerray(data[..., 2], u.G, 'B_phi', r'$B_\phi$', 'coolwarm',
                    [-1e15, 1e15], True)

@get_grid
@smooth
@derive
def poloidal_magnetic_fields(self, file_name, **kwargs):
    Br, Btheta, _ = self.magnetic_fields(file_name)
    data = np.sqrt(Br ** 2 + Btheta ** 2)
    data.set(label=r'$B_\mathrm{pol}$', name='B_pol', limits=[-1e15, 1e15],
             log=True, cmap='inferno')
    return data

@get_grid
@smooth
@derive
def toroidal_magnetic_fields(self, file_name, **kwargs):
    _, _, Bphi = self.magnetic_fields(file_name)
    return Bphi

@get_grid
@smooth
def magnetic_energy(self, file_name, comp: Literal['all', 'tot', 'pol', 'tor']='all',
                    **kwargs):
    """
    Magnetic energy density. Total, poloidal and toroidal.
    """
    Br, Btheta, Bphi = self.magnetic_fields(file_name)
    tot_b = 0.5 * (Br ** 2 + Btheta ** 2 + Bphi ** 2) / c.mu0
    tot_b.set(label=r'$E_\mathrm{mag}$', name='E_mag_tot', limits=[1e20, 1e28],
             log=True, cmap='magma')
    pol_b = 0.5 * (Br ** 2 + Btheta ** 2) / c.mu0
    pol_b.set(label=r'$E_\mathrm{mag,pol}$', name='E_mag_pol',
              limits=[1e20, 1e28], log=True, cmap='magma')
    tor_b = 0.5 * Bphi / c.mu0
    tor_b.set(label=r'$E_\mathrm{mag,tor}$', name='E_mag_tor',
              limits=[1e20, 1e28], log=True, cmap='magma')
    if comp == 'all':
        return  tot_b, pol_b, tor_b
    elif comp == 'pol':
        return pol_b
    elif comp == 'tor':
        return tor_b

@get_grid
@smooth
def stream_function(self, file_name, plane):
    return strfct2D(self.__CT_magnetic_fields(file_name), self.cell, 
                    self.ghost, plane)

@get_grid
@smooth
@derive
def alfven_velocity(self, file_name, **kwargs):
    """
    Alfven velocity
    """
    if self.dim == 1:
        return None
    Br, Btheta, Bphi = self.magnetic_fields(file_name)
    data = np.sqrt((Br ** 2 + Btheta ** 2 + Bphi ** 2) / 
                   self.rho(file_name) / c.mu0)
    data.set(label=r'$v_\mathrm{A}$', name='alfven_velocity',
              limits=[1e4, 1e7], log=True, cmap='gnuplot')
    return data