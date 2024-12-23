from AeViz.simulation.methods import *
from AeViz.utils.math_utils import strfct2D

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
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['mag_CT/data'][...]), self.dim)

@smooth
@hdf_isopen
def magnetic_fields(self, file_name, **kwargs):
    """
    Magnetic field at the cells center.
    """
    return self.ghost.remove_ghost_cells(np.squeeze(
        self._Simulation__data_h5['mag_vol/data'][...]), self.dim)

@smooth
@derive
def poloidal_magnetic_fields(self, file_name, **kwargs):
    data = self.magnetic_fields(file_name, **kwargs)
    return np.sqrt(data[..., 0] ** 2 + data[..., 1] ** 2)

@smooth
@derive
def toroidal_magnetic_fields(self, file_name, **kwargs):
    data = self.magnetic_fields(file_name, **kwargs)
    return data[..., 2]

@smooth
def magnetic_energy(self, file_name, **kwargs):
    """
    Magnetic energy density. Total, poloidal and toroidal.
    """
    data = self.magnetic_fields(file_name, **kwargs)
    return  0.5 * (data[..., 0] ** 2 + data[..., 1] ** 2 \
                + data[..., 2] ** 2), \
            0.5 * (data[..., 0] ** 2 + data[..., 1] ** 2), \
            0.5 * data[..., 2] ** 2

@smooth
def stream_function(self, file_name, plane):
    return strfct2D(self.__CT_magnetic_fields(file_name), self.cell, 
                    self.ghost, plane)

@smooth
@derive
def alfven_velocity(self, file_name, **kwargs):
    """
    Alfven velocity
    """
    if self.dim == 1:
        return None
    B = self.magnetic_fields(file_name, **kwargs)
    return np.sqrt(B[..., 0] ** 2 + B[..., 1] ** 2 + B[..., 2] ** 2) / \
        np.sqrt(self.rho(file_name, **kwargs))