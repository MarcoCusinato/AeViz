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
@mask_points
@smooth
@derive
@finite_differences
@hdf_isopen
def rho(self, file_name, **kwargs):
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['hydro/data']
        [..., self.hydroTHD_index['hydro']['I_RH']]), self.dim)
    return aerray(data, u.g / u.cm**3, 'density', r'$\rho$', 'viridis',
                  [1e4, 1e15], log=True)

@get_grid
@mask_points
@smooth
@derive
def mass(self, file_name, **kwargs):
    data = self.rho(file_name) * self.cell.dVolume_integration(self.ghost)
    data = data.to(u.Msun)
    data.set(name='mass', label=r'$M$', cmap='cividis', limits=[1e-10, 1e-1],
             log=True)
    return data

## ENERGY
@get_grid
@smooth
@finite_differences
@hdf_isopen
def MHD_energy(self, file_name, **kwargs):
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['hydro/data']\
            [..., self.hydroTHD_index['hydro']['I_EN']]), self.dim)
    return aerray(data, u.erg / u.cm**3, 'MHD_energy',
                  r'$E$', 'nipy_spectral', [1e24, 1e35], log=True)

@get_grid
@smooth
@finite_differences
@hdf_isopen
def rotational_energy(self, file_name, **kwargs):
    data = self.rho(file_name) * self.phi_velocity(file_name) ** 2 * 0.5
    data.to(u.erg/u.cm**3)
    data.set(name='rot_ene', label=r'$E_\mathrm{rot}$', cmap='viridis', log=True,
             limits=[1e24, 1e35])
    return data

@get_grid
@smooth
@finite_differences
@hdf_isopen
def kinetic_energy(self, file_name, comp:Literal['tot', 'r', 'th', 'ph']='tot',
                   **kwargs):
    if comp == 'tot':
        v = self.radial_velocity(file_name) ** 2 + \
            self.theta_velocity(file_name) ** 2 + \
                self.phi_velocity(file_name) ** 2
        label=r'$E_\mathrm{kin,tot}$'
    elif comp == 'r':
        v = self.radial_velocity(file_name) ** 2
        label=r'$E_{\mathrm{kin},r}$'
    elif comp == 'th':
        v =  self.theta_velocity(file_name) ** 2
        label=r'$E_{\mathrm{kin},\theta}$'
    elif comp == 'ph':
        v = self.phi_velocity(file_name) ** 2
        label=r'$E_{\mathrm{kin},\phi}$'
    data = self.rho(file_name) * v ** 2 * 0.5
    data.to(u.erg / u.cm**3)
    data.set(name='kin_ene', label=label, cmap='turbo', log=True,
             limits=[1e24, 1e35])
    return data

## VELOCITY
@get_grid
@mask_points
@smooth
@derive
@finite_differences
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
@mask_points
@smooth
@derive
@finite_differences
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
@mask_points
@smooth
@derive
@finite_differences
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
@mask_points
@smooth
@derive
@finite_differences
@hdf_isopen
def soundspeed(self, file_name, **kwargs):
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_CSND']]), self.dim)
    return aerray(data, u.cm / u.s, 'soundspeed', r'$c_\mathrm{s}$',
                  'nipy_spectral', [1e8, 1e10], log=True)

@get_grid
@mask_points
@smooth
@finite_differences
@derive
def omega(self, file_name, **kwargs):
    assert self.dim in [2, 3], "No rotational velocity in 1D"
    if self.dim == 2:
        data = self.phi_velocity(file_name) / \
            (np.sin(self.cell.theta(self.ghost))[:, None] * \
                self.cell.radius(self.ghost)[None, :])
    else:
        data = self.phi_velocity(file_name) / \
            (np.sin(self.cell.theta(self.ghost))[None, :, None] * \
                self.cell.radius(self.ghost)[None, None, :])
    return aerray(data, u.rad / u.s, 'omega', r'$\Omega$', 'magma',
                  [1e0, 1e4], log=True)