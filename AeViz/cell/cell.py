import numpy as np
import os

class cell:
    """
    Class that allows to get the grid related quantities
    Initialization parameters:
        path_folder: (string) path to the simulation folder
        dim: (int, optional) dimension of the supernova simulation (1, 2, or 3)
    """
    def __init__(self, path_folder=None, dim=None,
                 radius=None, theta=None, phi=None):
        assert dim in (1, 2, 3, None), "Supernova simulation can either be 1D, 2D or 3D"
        assert (path_folder is not None) or (radius is not None and theta is not None and phi is not None), \
            "Please provide either path to the simulation folder or the grid coordinates"
        if path_folder is not None:
            self.path_grid = os.path.join(path_folder, 'grid')
            self.__nu_grid_file = np.loadtxt(os.path.join(self.path_grid, 'grid.e.dat'))[:, 1:]
            self.__radius_file = np.loadtxt(os.path.join(self.path_grid, 'grid.x.dat'))[:, 1:]
            try:
                self.__theta_file = np.loadtxt(os.path.join(self.path_grid, 'grid.y.dat'))[:, 1:]
            except:
                self.__theta_file = np.loadtxt(os.path.join(self.path_grid, 'grid.y.dat'))[1:]
            try:
                self.__phi_file = np.loadtxt(os.path.join(self.path_grid, 'grid.z.dat'))[:, 1:]
            except:
                self.__phi_file = np.loadtxt(os.path.join(self.path_grid, 'grid.z.dat'))[1:]
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

    def simulation_dimension(self):
        """
        Method that return the simulation dimension
        result
            simulatio dimension (int)
        """
        return self.dim
    #radius methods
    def radius_left(self, ghost):
        """
        Method that returns an array with the 'left' coordinates of the radius, with the 
        selected number of ghost cells
        Parameters:
            ghost: (object) ghost
        Results:
            left radius coordinates: (numpy array)
        """
        return ghost.remove_ghost_cells(self.__radius_file[:, 0], self.dim, 'radius')
   
    def radius_right(self, ghost):
        """
        Method that returns an array with the 'right' coordinates of the radius, with the 
        selected number of ghost cells
        Parameters:
            ghost: (object) ghost
        Results:
            right radius coordinates: (numpy array)
        """
        return ghost.remove_ghost_cells(self.__radius_file[:, 2], self.dim, 'radius')

    def radius(self, ghost):
        """
        Method that returns an array with the 'center' coordinates of the radius, with the 
        selected number of ghost cells
        Parameters:
            ghost: (object) ghost
        Results:
            center radius coordinates: (numpy array)
        """
        return ghost.remove_ghost_cells(self.__radius_file[:, 1], self.dim, 'radius')

    def dr(self, ghost):
        """
        Method that returns an array with the lenght of each radial cell
        Parameters:
            ghost: (object) ghost
        Results:
            radial lenght of a cell: (numpy array)
        """
        return self.radius_right(ghost) - self.radius_left(ghost)

    #theta angle methods
    def theta_left(self, ghost):
        """
        Method that returns an array with the 'left' coordinates of the theta angle, with the 
        selected number of ghost cells
        Parameters:
            ghost: (object) ghost
        Results:
            left theta coordinates: (numpy array)
        """
        try:
            return ghost.remove_ghost_cells(self.__theta_file[:, 0], self.dim, 'theta')
        except:
            return self.__theta_file[0]

    def theta_right(self, ghost):
        """
        Method that returns an array with the 'right' coordinates of the theta angle, with the 
        selected number of ghost cells
        Parameters:
            ghost: (object) ghost
        Results:
            right theta coordinates: (numpy array)
        """
        try:
            return ghost.remove_ghost_cells(self.__theta_file[:, 2], self.dim, 'theta')
        except:
            return self.__theta_file[2]

    def theta(self, ghost):
        """
        Method that returns an array with the 'central' coordinates of the theta angle, with the 
        selected number of ghost cells
        Parameters:
            ghost: (object) ghost
        Results:
            central theta coordinates: (numpy array)
        """
        try:
            return ghost.remove_ghost_cells(self.__theta_file[:, 1], self.dim, 'theta')
        except:
            return self.__theta_file[1]
    
    def dtheta(self, ghost):
        """
        Method that returns an array with the angular (theta) lenght of each cell
        Parameters:
            ghost: (object) ghost
        Results:
            angular (theta) lenght of a cell: (numpy array)
        """
        return self.theta_right(ghost) - self.theta_left(ghost)
    
     #phi angle methods
    def phi_left(self, ghost):
        """
        Method that returns an array with the 'left' coordinates of the phi angle, with the 
        selected number of ghost cells
        Parameters:
            ghost: (object) ghost
        Results:
            left phi coordinates: (numpy array)
        """
        try:
            return ghost.remove_ghost_cells(self.__phi_file[:, 0], self.dim, 'phi')
        except:
            return self.__phi_file[0]

    def phi_right(self, ghost):
        """
        Method that returns an array with the 'left' coordinates of the phi angle, with the 
        selected number of ghost cells
        Parameters:
            ghost: (object) ghost
        Results:
            right phi coordinates: (numpy array)
        """
        try:
            return ghost.remove_ghost_cells(self.__phi_file[:, 2], self.dim, 'phi')
        except:
            return self.__phi_file[2]
    
    def phi(self, ghost):
        """
        Method that returns an array with the 'left' coordinates of the phi angle, with the 
        selected number of ghost cells
        Parameters:
            ghost: (object) ghost
        Results:
            central phi coordinates: (numpy array)
        """
        try:
            return ghost.remove_ghost_cells(self.__phi_file[:, 1], self.dim, 'phi')
        except:
            return self.__phi_file[1]

    def dphi(self, ghost):
        """
        Method that returns an array with the integration angular (phi) element in
        every dimension.
        Parameters:
            ghost: (object) ghost
        Results:
            dtheta for integration: (numpy array)
        """
        if self.dim != 3:
            return 2 * np.pi
        return self.phi_right(ghost) - self.phi_left(ghost)
    #lenght methods
    def lx(self, ghost):
        """
        Method that gives an array of cell's radial lengths.
        Parameters:
            ghost: (object) ghost
        Results:
            radial cell's lengths: (numpy array)
        """
        dr = self.dr(ghost)
        if self.dim == 1:
            return dr
        theta = np.ones(self.theta(ghost).shape[0])
        if self.dim == 2:
            return theta[:, None] * dr[None, :]
        else:
            phi = np.ones(self.phi(ghost).shape[0])
            return phi[:, None, None] * theta[None, :, None] * dr[None, None, :]

    def ly(self, ghost):
        """
        Method that gives an array of cell's theta angle lengths.
        Parameters:
            ghost: (object) ghost
        Results:
            radial cell's angle lengths: (numpy array)
        """
        if self.dim < 2:
            return None
        dtheta = self.dtheta(ghost)
        r = self.radius_left(ghost)
        if self.dim == 2:
            return dtheta[:, None] * r[None, :]
        else:
            phi = np.ones(self.phi(ghost).shape[0])
            return phi[:, None, None] * dtheta[None, :, None] * r[None, None, :]

    def lz(self, ghost):
        """
        Method that gives an array of cell's phi angle lengths.
        Parameters:
            ghost: (object) ghost
        Results:
            radial cell's angle lengths: (numpy array)
        """
        if self.dim < 2:
            return None
        r = self.radius_left(ghost)
        theta = np.sin(self.theta_left(ghost))
        dphi = self.dphi(ghost)
        if self.dim == 2:
            return dphi * theta[:, None] * r[None, :]
        else:
            return dphi[:, None, None] * theta[None, :, None] * r[None, None, :]
    #surface methods
    def ax(self, ghost):
        """
        Method that gives an array of cell's surface normal to the radius.
        Parameters:
            ghost: (object) ghost
        Results:
            surfaces normal to the radius: (numpy array)
        """
        if self.dim < 2:
            return None
        r = self.radius_left(ghost)**2
        theta = np.sin(self.theta_left(ghost))
        dphi = self.dphi(ghost)
        if self.dim == 2:
            return dphi * theta[:, None] * r[None, :]
        else:
            return dphi[:, None, None] * theta[None, :, None] * r[None, None, :]

    def ay(self, ghost):
        """
        Method that gives an array of cell's surface normal to the theta angle.
        Parameters:
            ghost: (object) ghost
        Results:
            surfaces normal to the theta angle: (numpy array)
        """
        if self.dim < 2:
            return None
        r = self.radius_left(ghost) * self.dr(ghost)
        theta = np.sin(self.theta_left(ghost))
        dphi = self.dphi(ghost)
        if self.dim == 2:
            return dphi * theta[:, None] * r[None, :]
        else:
            return dphi[:, None, None] * theta[None, :, None] * r[None, None, :]

    def az(self, ghost):
        """
        Method that gives an array of cell's surface normal to the phi angle.
        Parameters:
            ghost: (object) ghost
        Results:
            surfaces normal to the phi angle: (numpy array)
        """
        if self.dim < 2:
            return None
        r = self.radius_left(ghost) * self.dr(ghost)
        dtheta = self.dtheta(ghost)
        if self.dim == 2:
            return dtheta[:, None] * r[None, :]
        else:
            phi = np.ones(self.phi(ghost))
            return phi[:, None, None] * dtheta[None, :, None] * r[None, None, :]
    
    #integration methods
    def dVolume_integration(self, ghost):
        dr = self.dr_integration(ghost)
        dtheta = self.dtheta_integration(ghost)
        dphi = self.dphi(ghost)
        if self.dim == 1:
            return dr*dtheta*dphi
        if self.dim == 2:
            return dphi * dtheta[:, None] * dr[None, :]
        return dphi[:, None, None] * dtheta[None, :, None] * dr[None, None, :]

    def dVolume_sum(self, ghost):
        r = self.radius(ghost) ** 2
        if self.dim == 1:
            return 4 * np.pi * r
        theta = np.sin(self.theta(ghost))
        if self.dim == 2:
            return 2 * theta[:, None] * r[None, :]
        phi = self.phi(ghost)
        return phi[:, None, None] * theta[None, :, None] * r[None, None, :]

    def dr_integration(self, ghost):
        """
        Method that returns an array with the integration radial element
        Parameters:
            ghost: (object) ghost
        Results:
            dr for integration: (numpy array)
        """
        return (self.radius_right(ghost)**3 - self.radius_left(ghost)**3) / 3
    
    def dtheta_integration(self, ghost):
        """
        Method that returns an array with the integration angular (theta) element in
        every dimension.
        Parameters:
            ghost: (object) ghost
        Results:
            dtheta for integration: (numpy array)
        """
        if self.dim < 2:
            return 2
        return np.cos(self.theta_left(ghost)) - np.cos(self.theta_right(ghost))
    
    def dOmega(self, ghost):
        """
        Method that returns an array with the solid angle integration element.
        Parameters:
            ghost: (object) ghost
        Results:
            dOmega: (numpy array)
        """
        return self.dtheta_integration(ghost) * self.dphi(ghost)
    
    def E_nu_left(self):
        """
        Method that returns an array with the left neutrino energy grid.
        Results:
            E_nu_left: (numpy array)
        """
        if self.path_grid is None:
            raise TypeError("You shouldn't be here.")
        return self.__nu_grid_file[:, 0]
    
    def E_nu_right(self):
        """
        Method that returns an array with the right neutrino energy grid.
        Results:
            E_nu_right: (numpy array)
        """
        if self.path_grid is None:
            raise TypeError("You shouldn't be here.")
        return self.__nu_grid_file[:, 2]
    
    def E_nu(self):
        """
        Method that returns an array with the neutrino energy grid.
        Results:
            E_nu: (numpy array)
        """
        if self.path_grid is None:
            raise TypeError("You shouldn't be here.")
        return self.__nu_grid_file[:, 1]

    def dE_nu(self):
        """
        Method that returns an array with the neutrino energy integration element.
        Results:
            dE_nu: (numpy array)
        """
        return self.E_nu_right() - self.E_nu_left()
