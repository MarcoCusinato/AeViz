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
@smooth
@derive
@hdf_isopen
def gas_pressure(self, file_name, **kwargs):
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_PGAS']]), self.dim)

@smooth
@derive
@hdf_isopen
def temperature(self, file_name, **kwargs):
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_TMPR']]), self.dim)

@smooth
@derive
@hdf_isopen
def enthalpy(self, file_name, **kwargs):
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_ENTH']]), self.dim)

@smooth
@derive
@hdf_isopen
def entropy(self, file_name, **kwargs):
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_ENTR']]), self.dim)

@smooth
@derive
@hdf_isopen
def adiabatic_index(self, file_name, **kwargs):
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_GAMM']]), self.dim)

## RELATIVITY AND GRAVITY
@smooth
@derive
@hdf_isopen
def lorentz(self, file_name, **kwargs):
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_LRTZ']]), self.dim)

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
    return self.ghost.remove_ghost_cells(data, self.dim)

@smooth
@derive
def gravitational_energy(self, file_name, **kwargs):
    return 0.5 * self.rho(file_name, **kwargs) * \
        self.gravitational_potential(file_name, **kwargs)

## ENERGY
@smooth
@derive
@hdf_isopen
def internal_energy(self, file_name, **kwargs):
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_EINT']]), self.dim)

@smooth
@derive
@hdf_isopen
def nu_heat(self, file_name, **kwargs):
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_HEAT']]), self.dim)