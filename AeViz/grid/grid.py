from scipy.interpolate import (griddata, interp1d)
import numpy as np

class grid:
    """
    Class to compute a cartesian grid from spherical coordinates.
    BEWARE: is also possible to obtain a finer (o larger) grid and then
    interpolate the desired quantity to it, HOWEVER, it works only for
    1D and 2D cases.
    parameters:
        dim: (int) simulation dimension, MUST be in (1, 2, 3)
        radius: (float, numpy array) radius of the spherical simulation
        theta: (float, numpy array) polar angle
        phi: (float, numpy array) azimutal angle
    """
    def __init__(self, dim, radius, theta = None, phi = None):
        assert dim in (1, 2, 3), "Supernova simulation MUST be 1, 2 or 3D."
        if dim in (2, 3) and theta is None:
            raise TypeError("Theta angle MUST not be None in 2D or 3D " \
                "simulations.")
        if dim == 3 and phi is None:
            raise TypeError("Phi angle MUST not be None in 3D simulations.")
        self.dim = dim
        self.radius = radius
        self.theta = theta
        self.phi = phi
        if self.dim == 1:
            self.__default_grid_parameters = {'interpolation_method': 'cubic',
                                    'r_0': self.radius[0],
                                    'r_0_log': self.radius[50],
                                    'r_max': self.radius[-1],
                                    'n_r_lin': int(0.5*self.radius.shape[0]),
                                    'n_r_log': int(0.5*self.radius.shape[0])}
        elif self.dim == 2:
            self.__default_grid_parameters = {'interpolation_method': 'cubic',
                                    'r_0': self.radius[0],
                                    'r_0_log': self.radius[50],
                                    'r_max': self.radius[-1],
                                    'n_r_lin': int(0.5*self.radius.shape[0]),
                                    'n_r_log': int(0.5*self.radius.shape[0]),
                                    'theta_0': self.theta[0],
                                    'theta_max': self.theta[-1],
                                    'n_theta': self.theta.shape[0]}
        else:
            self.__default_grid_parameters = {'interpolation_method': 'cubic',
                                    'r_0': self.radius[0],
                                    'r_0_log': self.radius[50],
                                    'r_max': self.radius[-1],
                                    'n_r_lin': int(0.5*self.radius.shape[0]),
                                    'n_r_log': int(0.5*self.radius.shape[0]),
                                    'theta_0': self.theta[0],
                                    'theta_max': self.theta[-1],
                                    'n_theta': self.theta.shape[0],
                                    'phi_0': self.phi[0],
                                    'phi_max': self.phi[-1],
                                    'n_phi': self.phi.shape[0]}

    def cartesian_grid(self):
        """
        Returns grid-like arrays (similar to the ones returned by
        mehgrid).
        Returns  1, 2 or 3 arrays depending on the simulation dimension.
        """
        if self.dim == 1:
            return self.__1D_cartesian_grid(self.radius)
        elif self.dim == 2:
            return self.__2D_cartesian_grid(self.radius, self.theta)
        else:
            return self.__3D_cartesian_grid(self.radius, self.theta, self.phi)
    
    def velocity_sph_to_cart(self, v_r=None, v_theta=None, v_phi=None):
        """
        Returns the cartesian components of the velocity.
        """
        if self.dim == 1:
            return v_r, 0, 0
        if self.dim == 2:
            return self.__2D_velocity_sph_to_cart(v_r, v_theta, v_phi)
        return self.__3D_velocity_sph_to_cart(v_r, v_theta, v_phi)

    def new_cartesian_grid(self, **kwargs):
        """
        Returns a cartesian grid with the options provided.
        Possible options:
        'r_0': minimum radius value
        'r_0_log': start value of logaritmic radius spacing
        'r_max': maximum radius value
        'n_r_lin': number of linear points
        'n_r_log': number of logaritmic points
        'theta_0': minimum value of the azimuthal angle
        'theta_max': maximum value of the azimuthal angle
        'n_theta': number of azimuthal points
        'phi_0': minimum value of the polar angle
        'phi_max': maximum value of the polar angle
        'n_phi': number of polar points
        """
        if not kwargs:
            return self.cartesian_grid()
        updated_parameters = self.__default_grid_parameters.copy()
        updated_parameters.update(kwargs) 
        r_new = self.__new_radius(updated_parameters)
        if self.dim == 1:
            return self.__1D_cartesian_grid(r_new)
        theta_new = self.__new_theta(updated_parameters)
        if self.dim == 2:
            return self.__2D_cartesian_grid(r_new, theta_new)
        phi_new = self.__new_phi(updated_parameters)
        return self.__3D_cartesian_grid(r_new, theta_new, phi_new)
    
    def interpolate_quantity(self, grid, quantity, new_grid):
        if self.dim == 1:
            f = interp1d(grid, quantity, 'cubic')
            return f(new_grid)
        if self.dim == 2:    
            points = np.array( (grid[0].flatten(), grid[1].flatten()) ).T
            values = (quantity.T).flatten()
            return griddata(points, values, (new_grid[0], new_grid[1]),
                            'cubic')
        if self.dim == 3:
            raise TypeError("Not implemented yet :\").")
    
    def cartesian_grid_2D(self, plane):
        """
        Returns a 2D cartesian grid.
        """
        if self.dim == 1:
            return self.__1D_cartesian_grid(self.radius)
        elif self.dim == 2:
            return self.__2D_cartesian_grid(self.radius, self.theta)
        else:
            return self.__2D_cartesian_grid_from_3D(self.radius, self.theta,
                                                     self.phi, plane)
    
    def __1D_cartesian_grid(self, radius):
        return radius

    def __2D_cartesian_grid(self, radius, theta):
        X = (radius[None, :] * np.sin(theta)[:, None])
        Y = (radius[None, :] * np.cos(theta)[:, None])
        return X, Y
    
    def __2D_cartesian_grid_from_3D(self, radius, theta, phi, plane):
        if plane == 'xy':
            X = radius[None, :] * np.sin(theta[len(theta) // 2]) * \
                np.cos(phi)[:, None]
            Y = radius[None, :] * np.sin(theta[len(theta) // 2]) * \
                np.sin(phi)[:, None]
        elif plane in ['xz', 'yz']:
            theta = np.concatenate([theta, theta + np.pi])
            X = radius[None, :] * np.sin(theta)[:, None]
            Y = radius[None, :] * np.cos(theta)[:, None]
        return X, Y

    def __3D_cartesian_grid(self, radius, theta, phi):
        X = radius[None, None, :] * np.sin(theta)[None, :, None] * \
            np.cos(phi)[:, None, None]
        Y = radius[None, None, :] * np.sin(theta)[None, :, None] * \
            np.sin(phi)[:, None, None]
        Z = radius[None, None, :] * np.cos(theta)[None, :, None] * \
            np.ones(phi.shape)[:, None, None]
        return X, Y, Z
    
    def __new_radius(par):
        r = np.linspace(par['r_0'], par['r_0_log'], par['n_r_lin'],
                        endpoint=False)
        r_log = 10**np.linspace(np.log10(par['r_0_log']),
                                np.log10(par['r_max_log']), 
                                par['n_r_log']+1)
        return np.concatenate((r, r_log))
    
    def __new_theta(par):
        return np.linspace(par['theta_0'], par['theta_max'], par['n_theta'])

    def __new_phi(par):
        return np.linspace(par['phi_0'], par['phi_max'], par['n_phi'])
    
    def __2D_velocity_sph_to_cart(self, v_r, v_theta, v_phi):
        Vx = v_r * np.sin(self.theta)[:, None] + \
            v_theta * np.cos(self.theta)[:, None]
        Vy = v_phi * np.sin(self.theta)[:, None]
        Vz = v_r * np.cos(self.theta)[:, None] - \
            v_theta * np.sin(self.theta)[:, None]
        return Vx, Vy, Vz
    
    def __3D_velocity_sph_to_cart(self, v_r, v_theta, v_phi):
        Vx = v_r * np.sin(self.theta)[None, :, None] * \
                            np.cos(self.phi)[:, None, None] + \
            v_theta * np.cos(self.theta)[None, :, None] * \
                            np.cos(self.phi)[:, None, None] - \
            v_phi * np.sin(self.phi)[:, None, None]
        Vy = v_r * np.sin(self.theta)[None, :, None] * \
                            np.sin(self.phi)[:, None, None] + \
            v_theta * np.cos(self.theta)[None, :, None] * \
                            np.sin(self.phi)[:, None, None] + \
            v_phi * np.cos(self.phi)[:, None, None]
        Vz = v_r * np.cos(self.theta)[None, :, None] - \
            v_theta * np.sin(self.theta)[None, :, None]
        return Vx, Vy, Vz