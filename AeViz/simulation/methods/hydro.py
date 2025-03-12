from AeViz.simulation.methods import *


"""
This module contains all the local hydro quantities.
These functions are not supposed to be used standalone, but as methods
of the Simulation class
"""

## -----------------------------------------------------------------
## HYDRODYNAMICAL DATA
## -----------------------------------------------------------------

@get_grid
@smooth
@derive
@hdf_isopen
def rho(self, file_name, **kwargs):
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['hydro/data']
        [..., self.hydroTHD_index['hydro']['I_RH']]), self.dim)
    return aerray(data, u.g / u.cm**3, 'density', r'$\rho$', 'viridis',
                  [1e4, 1e15], log=True)

## ENERGY
@get_grid
@smooth
@hdf_isopen
def MHD_energy(self, file_name, **kwargs):
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['hydro/data']\
            [..., self.hydroTHD_index['hydro']['I_EN']]), self.dim)
    return aerray(data, u.erg / u.cm**3, 'internal_energy',
                  r'$E$', 'nipy_spectral', [1e24, 1e35], log=True)

## VELOCITY
@get_grid
@smooth
@derive
@hdf_isopen
def radial_velocity(self, file_name, **kwargs):
    if self.relativistic:
        ivx = 'I_VELX'
    else:
        ivx = 'I_VX'
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    [ivx]]), self.dim)
    return aerray(data, u.cm / u.s, 'velocity_radial', r'$v_r$', 'Spectral_r',
                  [-1e10, 5e10], log=True)

@get_grid
@smooth
@derive
@hdf_isopen
def theta_velocity(self, file_name, **kwargs):
    if self.relativistic:
        ivy = 'I_VELY'
    else:
        ivy = 'I_VY'
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    [ivy]]), self.dim)
    return aerray(data, u.cm / u.s, 'velocity_theta', r'$v_\theta$',
                  'Spectral_r', [-1e10, 5e10], log=True)

@get_grid
@smooth
@derive
@hdf_isopen
def phi_velocity(self, file_name, **kwargs):
    if self.relativistic:
        ivz = 'I_VELZ'
    else:
        ivz = 'I_VZ'
    data =  self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    [ivz]]), self.dim)
    return aerray(data, u.cm / u.s, 'velocity_phi', r'$v_\phi$',
                  'cividis', [1e7, 5e10], log=True)


@get_grid
@smooth
@derive
@hdf_isopen
def soundspeed(self, file_name, **kwargs):
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_CSND']]), self.dim)
    return aerray(data, u.cm / u.s, 'soundspeed', r'$c_\mathrm{s}$',
                  'nipy_spectral', [1e8, 1e10], log=True)

@get_grid
@smooth
@derive
def omega(self, file_name, **kwargs):
    assert self.dim in [2, 3], "No rotational velocity in 1D"
    if self.dim == 2:
        data = self.phi_velocity(file_name, **kwargs) / \
            (np.sin(self.cell.theta(self.ghost))[:, None] * \
                self.cell.radius(self.ghost)[None, :])
    else:
        data = self.phi_velocity(file_name, **kwargs) / \
            (np.cos(self.cell.phi(self.ghost))[:, None, None] *
                np.sin(self.cell.theta(self.ghost))[None, :, None] * \
                self.cell.radius(self.ghost)[None, None, :])
    return aerray(data, u.rad / u.s, 'omega', r'$\Omega$', 'magma',
                  [1e0, 1e4], log=True)