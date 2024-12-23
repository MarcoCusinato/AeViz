from AeViz.simulation.methods import *


"""
This module contains all the local hydro quantities.add()
These functions are not supposed to be used standalone, but as methods
of the Simulation class
"""

## -----------------------------------------------------------------
## HYDRODYNAMICAL DATA
## -----------------------------------------------------------------

@smooth
@derive
@hdf_isopen
def rho(self, file_name, **kwargs):
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['hydro/data']\
            [..., self.hydroTHD_index['hydro']['I_RH']]), self.dim)

## ENERGY
@smooth
@hdf_isopen
def MHD_energy(self, file_name, **kwargs):
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['hydro/data']\
            [..., self.hydroTHD_index['hydro']['I_EN']]), self.dim)

## VELOCITY
@smooth
@derive
@hdf_isopen
def radial_velocity(self, file_name, **kwargs):
    if self.relativistic:
        ivx = 'I_VELX'
    else:
        ivx = 'I_VX'
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    [ivx]]), self.dim)

@smooth
@derive
@hdf_isopen
def theta_velocity(self, file_name, **kwargs):
    if self.relativistic:
        ivy = 'I_VELY'
    else:
        ivy = 'I_VY'
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    [ivy]]), self.dim)

@smooth
@derive
@hdf_isopen
def phi_velocity(self, file_name, **kwargs):
    if self.relativistic:
        ivz = 'I_VELZ'
    else:
        ivz = 'I_VZ'
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    [ivz]]), self.dim)

@smooth
@derive
@hdf_isopen
def soundspeed(self, file_name, **kwargs):
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                    ['I_CSND']]), self.dim)

@smooth
@derive
def omega(self, file_name, **kwargs):
    assert self.dim in [2, 3], "No rotational velocity in 1D"
    if self.dim == 2:
        return self.phi_velocity(file_name, **kwargs) / \
            (np.sin(self.cell.theta(self.ghost))[:, None] * \
                self.cell.radius(self.ghost)[None, :])
    else:
        return self.phi_velocity(file_name, **kwargs) / \
            (np.cos(self.cell.phi(self.ghost))[:, None, None] *
                np.sin(self.cell.theta(self.ghost))[None, :, None] * \
                self.cell.radius(self.ghost)[None, None, :])
