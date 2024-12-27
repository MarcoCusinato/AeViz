from AeViz.utils.file_utils import list_module_functions
import numpy as np
import os
import types

class cell:
    """
    Class that allows to get the grid related quantities
    Initialization parameters:
        path_folder: (string) path to the simulation folder
        dim: (int, optional) dimension of the supernova simulation (1, 2, or 3)
    """
    def __init__(self, path_folder=None, dim=None,
                 radius=None, theta=None, phi=None,
                 neu=None, geom=None):
        assert dim in (1, 2, 3, None), "Supernova simulation can either ' \
            'be 1D, 2D or 3D"
        assert (path_folder is not None) or \
            (radius is not None and theta is not None and phi is not None), \
            "Please provide either path to the simulation folder ' \
                'or the grid coordinates"
        assert geom is not None, "Please provide the geometry of the simulation"
        assert geom in [1, 2], "Geomtry must be either 1 or 2"
        if path_folder is not None:
            self.path_grid = os.path.join(path_folder, 'grid')
            if neu >= 1:
                self.__nu_grid_file = np.loadtxt(os.path.join(self.path_grid,
                                                          'grid.e.dat'))[:, 1:]
            self.__radius_file = np.loadtxt(os.path.join(self.path_grid,
                                                         'grid.x.dat'))[:, 1:]
            try:
                self.__theta_file = np.loadtxt(os.path.join(self.path_grid,
                                                            'grid.y.dat'))[:, 1:]
            except:
                self.__theta_file = np.loadtxt(os.path.join(self.path_grid,
                                                            'grid.y.dat'))[1:]
            try:
                self.__phi_file = np.loadtxt(os.path.join(self.path_grid,
                                                          'grid.z.dat'))[:, 1:]
            except:
                self.__phi_file = np.loadtxt(os.path.join(self.path_grid, 
                                                          'grid.z.dat'))[1:]
        else:
            assert radius.ndim == 2, "Radius must be a 3D array"
            assert theta.ndim == 2, "Theta must be a 3D array"
            assert phi.ndim == 2, "Phi must be a 3D array"
            self.__radius_file = radius
            self.__theta_file = theta
            self.__phi_file = phi

        if dim is None:
            dim = 1
            if (self.__theta_file).size > 4:
                dim += 1
            if (self.__phi_file).size > 4:
                dim += 1
        self.dim = dim
        if neu >= 1:
            self.__load_neu_methods()
        if geom == 1:
            self.__load_cartesian_methods()
        else:
            self.__load_spherical_methods()

    def simulation_dimension(self):
        """
        Method that return the simulation dimension
        result
            simulation dimension (int)
        """
        return self.dim
    
    def __load_neu_methods(self):
        """
        Method that loads the neutrino grid methods
        """
        import AeViz.cell.cell_methods.nu_methods
        funcs = list_module_functions(AeViz.cell.cell_methods.nu_methods)
        for name, obj in funcs:
            setattr(self, name, types.MethodType(obj, self))

    def __load_cartesian_methods(self):
        """
        Method that loads the cartesian grid methods
        """
        import AeViz.cell.cell_methods.cartesian_methods
        funcs = list_module_functions(AeViz.cell.cell_methods.cartesian_methods)
        for name, obj in funcs:
            setattr(self, name, types.MethodType(obj, self))

    def __load_spherical_methods(self):
        """
        Method that loads the spherical grid methods
        """
        import AeViz.cell.cell_methods.spherical_methods
        funcs = list_module_functions(AeViz.cell.cell_methods.spherical_methods)
        for name, obj in funcs:
            setattr(self, name, types.MethodType(obj, self))