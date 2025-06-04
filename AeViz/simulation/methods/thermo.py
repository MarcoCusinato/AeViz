from AeViz.simulation.methods import *

"""
Functions to handle thermodynamical data from a simulation.
These functions are not meant to be used standalone, but rather to be
imported into the Simulation class.
"""

## -----------------------------------------------------------------
## THERMODYNAMICAL DATA
## -----------------------------------------------------------------

## THERMODYNAMICAL
@get_grid
@smooth
@derive
@hdf_isopen
def gas_pressure(self, file_name, **kwargs):
    data =  self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_PGAS']]), self.dim)
    return aerray(data, u.Ba, name='gas_pressure', label=r'$P_\mathrm{gas}$',
                  cmap='gist_rainbow_r', limits=[1e25, 1e34], log=True)

@get_grid
@smooth
@derive
@hdf_isopen
def temperature(self, file_name, **kwargs):
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_TMPR']]), self.dim)
    return aerray(data, u.MeV, name='temperature', label=r'$T$',
                  cmap='inferno', limits=[0, 40], log=False)

@get_grid
@smooth
@derive
@hdf_isopen
def enthalpy(self, file_name, **kwargs):
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_ENTH']]), self.dim)
    return aerray(data, u.erg, 'enthalpy', r'$H$', 'gist_stern', [1e25, 1e36],
                  True)

@get_grid
@smooth
@derive
@hdf_isopen
def entropy(self, file_name, **kwargs):
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_ENTR']]), self.dim)
    return aerray(data, u.kBol / u.bry, 'entropy', r'$s$', 'gist_rainbow_r',
                  [1.5, 15], False)

@get_grid
@smooth
@derive
@hdf_isopen
def adiabatic_index(self, file_name, **kwargs):
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_GAMM']]), self.dim)
    return aerray(data, u.dimensionless_unscaled, 'adiabatic_index',
                  r'$\Gamma$', 'cividis', [0.5, 3.5], False)

## RELATIVITY AND GRAVITY
@get_grid
@smooth
@derive
@hdf_isopen
def lorentz(self, file_name, **kwargs):
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_LRTZ']]), self.dim)
    return aerray(data, u.dimensionless_unscaled, 'lorentz_factor',
                  r'$\gamma$', 'gist_rainbow', [1, 1.1], False)

@get_grid
@smooth
@derive
@hdf_isopen
def gravitational_potential(self, file_name, **kwargs):
    data = np.squeeze(
        self._Simulation__data_h5['gravpot/data'][...])
    ## The following IF is here to address problem in saving 3D
    ## arrays on MN4 
    if data.shape[-1] == 2:
        data = data[..., 0]
    data = self.ghost.remove_ghost_cells(data, self.dim)
    return aerray(data, u.erg / u.g, 'gravitational_potential', r'\Phi',
                  'magma', [-1e22, -1e15], True)

@get_grid
@smooth
@derive
def gravitational_energy(self, file_name, **kwargs):
    data = 0.5 * self.rho(file_name) * \
        self.gravitational_potential(file_name)
    data.set(label=r'$E_mathrm{grav}$', name='gravitational_energy', log=True,
             limits=[-1e22, -1e15], cmap='magma')
    return data

@get_grid
@smooth
@derive
def lapse_function(self, file_name, **kwargs):
    if self.lapse_form == 1:
        alpha = 1 + self.gravitational_potential(file_name) / c.c ** 2
    elif self.lapse_form == 2:
        alpha = np.exp(self.gravitational_potential(file_name) / c.c ** 2)
    elif self.lapse_form == 3:
        alpha = np.where(self.gravitational_potential(file_name) / c.c ** 2 >= -1,
                         np.exp(-(self.gravitational_potential(file_name) / c.c ** 2)**2))
    else:
        alpha = aerray(np.ones_like(self.gravitational_potential(file_name)),
                       u.dimensionless_unscaled)
    alpha.set(name='lapse', label=r'$\alpha$', cmap='rainbow', limits=[-1, 1],
              log=False)
    return alpha
    
## ENERGY
@get_grid
@smooth
@derive
@hdf_isopen
def internal_energy(self, file_name, **kwargs):
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_EINT']]), self.dim)
    return aerray(data, u.erg / u.cm**3, 'internal_energy',
                  r'$E_\mathrm{int}$', 'nipy_spectral', [1e24, 1e35], log=True)

@get_grid
@smooth
@derive
@hdf_isopen
def nu_heat(self, file_name, **kwargs):
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_HEAT']]), self.dim)
    return aerray(data, u.erg / u.cm ** 3, 'neutrino_heat', r'$Q_\nu$',
                  'Spectral_r', [-1e31, 1e32], True)