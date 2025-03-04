from AeViz.simulation.methods import *

"""
Functions to handle composition data from a simulation.
These functions are not meant to be used standalone, but rather to be
imported into the Simulation class.
"""

## -----------------------------------------------------------------
## COMPOSITION DATA
## -----------------------------------------------------------------

@smooth
@derive
@hdf_isopen
def Ye(self, file_name, **kwargs):
    data = aerray(self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['hydro/data'][..., self.hydroTHD_index['hydro']
        ['I_YE']]), self.dim), u.g / u.cm ** 3) / self.rho(file_name, **kwargs)
    data.set(name='Ye', label= r'$Y_\mathrm{e}$', limits=[0.0, 0.5],
             cmap='gist_rainbow', log=False)
    return data

@smooth
@derive
@hdf_isopen
def neutron_fraction(self, file_name, **kwargs):
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_COMP'][0]]), self.dim)
    return aerray(data, u.dimensionless_unscaled, 'neutron_fraction',
                  r'$X_\mathrm{n}$', 'cividis', [0, 1], False)

@smooth
@derive
@hdf_isopen
def proton_fraction(self, file_name, **kwargs):
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_COMP'][1]]), self.dim)
    return aerray(data, u.dimensionless_unscaled, 'proton_fraction',
                  r'$X_\mathrm{p}$', 'viridis', [0, 1], False)

@smooth
@derive
@hdf_isopen
def alpha_fraction(self, file_name, **kwargs):
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_COMP'][2]]), self.dim)
    return aerray(data, u.dimensionless_unscaled, 'alpha_fraction',
                  r'$X_\alpha$', 'plasma', [0, 1], False)

@smooth
@derive
@hdf_isopen
def heavy_fraction(self, file_name, **kwargs):
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_COMP'][3]]), self.dim)
    return aerray(data, u.dimensionless_unscaled, 'heavy_fraction',
                  r'$X_\mathrm{h}$', 'magma', [0, 1], False)

@smooth
@derive
@hdf_isopen
def Abar(self, file_name, **kwargs):
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_COMP'][4]]), self.dim)
    return aerray(data, u.dimensionless_unscaled, 'Abar',
                  r'$\bar{A}$', 'gist_rainbow_r', [4, 80], False)

@smooth
@derive
@hdf_isopen
def Zbar(self, file_name, **kwargs):
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_COMP'][5]]), self.dim)
    return aerray(data, u.dimensionless_unscaled, 'Zbar',
                  r'$\bar{Z}$', 'nipy_spectral', [1, 34], False)

## CHEMICAL POTENTIAL
@smooth
@derive
@hdf_isopen
def electron_chemical_potential(self, file_name, **kwargs):
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_CPOT'][0]]), self.dim)
    return aerray(data, u.erg / u.g, 'electron_chemical_potential',
                  r'$\mu_\mathrm{e}$', 'coolwarm', [0.1, 300], True)

@smooth
@derive
@hdf_isopen
def neutron_chemical_potential(self, file_name, **kwargs):
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_CPOT'][1]]), self.dim)
    return aerray(data, u.erg / u.g, 'neutron_chemical_potential',
                  r'$\mu_\mathrm{n}$', 'bwr', [-2e2, 2e3], True)

@smooth
@derive
@hdf_isopen
def proton_chemical_potential(self, file_name, **kwargs):
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_CPOT'][2]]), self.dim)
    return aerray(data, u.erg / u.g, 'proton_chemical_potential',
                  r'$\mu_\mathrm{p}$', 'seismic', [-2e2, 3e3], True)

@smooth
@derive
@hdf_isopen
def neutrino_chemical_potential(self, file_name, **kwargs):
    data = self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                        ['I_CPOT'][3]]), self.dim)
    return aerray(data, u.erg / u.g, 'neutrino_chemical_potential',
                  r'$\mu_\nu$', 'Spectral_r', [-2e2, 3e3], True)