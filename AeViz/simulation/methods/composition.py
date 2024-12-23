from AeViz.simulation.methods import *

## -----------------------------------------------------------------
## COMPOSITION DATA
## -----------------------------------------------------------------

@smooth
@derive
@hdf_isopen
def Ye(self, file_name, **kwargs):
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['hydro/data'][..., self.hydroTHD_index['hydro']
        ['I_YE']]), self.dim) / self.rho(file_name, **kwargs)

@smooth
@derive
@hdf_isopen
def neutron_fraction(self, file_name, **kwargs):
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_COMP'][0]]), self.dim)

@smooth
@derive
@hdf_isopen
def proton_fraction(self, file_name, **kwargs):
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_COMP'][1]]), self.dim)

@smooth
@derive
@hdf_isopen
def alpha_fraction(self, file_name, **kwargs):
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_COMP'][2]]), self.dim)

@smooth
@derive
@hdf_isopen
def heavy_fraction(self, file_name, **kwargs):
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_COMP'][3]]), self.dim)

@smooth
@derive
@hdf_isopen
def Abar(self, file_name, **kwargs):
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_COMP'][4]]), self.dim)

@smooth
@derive
@hdf_isopen
def Zbar(self, file_name, **kwargs):
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_COMP'][5]]), self.dim)

## CHEMICAL POTENTIAL
@smooth
@derive
@hdf_isopen
def electron_chemical_potential(self, file_name, **kwargs):
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_CPOT'][0]]), self.dim)

@smooth
@derive
@hdf_isopen
def neutron_chemical_potential(self, file_name, **kwargs):
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_CPOT'][1]]), self.dim)

@smooth
@derive
@hdf_isopen
def proton_chemical_potential(self, file_name, **kwargs):
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_CPOT'][2]]), self.dim)

@smooth
@derive
@hdf_isopen
def neutrino_chemical_potential(self, file_name, **kwargs):
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                        ['I_CPOT'][3]]), self.dim)