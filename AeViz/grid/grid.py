import numpy as np
from AeViz.units import aerray, u
from AeViz.utils.math_utils import IDL_derivative
from AeViz.utils.files.string_utils import merge_strings, apply_symbol
import copy

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
        self.__cart_grid = None
        self.__cell = None
        self.__ghost = None

    @classmethod
    def from_cell(cls, cl, gh):
        if cl.geom != 2:
            raise NotImplementedError("grid only supports spherical geometry.")
        gr = cls(cl.dim, cl.radius(gh),
                   cl.theta(gh) if cl.dim > 1 else None,
                   cl.phi(gh) if cl.dim > 2 else None)
        gr.__cell = cl
        gr.__ghost = copy(gh)
        return gr

    def _metric(self):
        if self.__cell is None:
            return TypeError("Grid built from raw coordinates. No cell support.")
        return self.__cell

    def cartesian_grid(self):
        """
        Returns grid-like arrays (similar to the ones returned by
        mehgrid).
        Returns  1, 2 or 3 arrays depending on the simulation dimension.
        """
        if self.__cart_grid is not None:
            return self.__cart_grid
        if self.dim == 1:
            self.__cart_grid = self.__1D_cartesian_grid(self.radius).value
        elif self.dim == 2:
            X, Y =  self.__2D_cartesian_grid(self.radius, self.theta, 'xz')
            self.__cart_grid = np.array([X.value, Y.value])
        else:
            X, Y, Z = self.__3D_cartesian_grid(self.radius, self.theta, self.phi)
            self.__cart_grid = np.array([X.value, Y.value, Z.value])
        return self.__cart_grid

    def spherical_to_cartesian(self, r=None, th=None, ph=None,
                               add_front=0, add_back=0):
        """
        Returns the cartesian components of the input array.
        input:
            r: (numpy array) radial component
            th: (numpy array) polar component
            ph: (numpy array) azimutal component
            add_front: (int) number of axis to add in front of the array
            add_back: (int) number of axis to add in back of the array
        """
        if self.dim == 1:
            return r
        if self.dim == 2:
            return self.__2D_spherical_to_cartesian(r, th, ph, add_front,
                                                     add_back)
        return self.__3D_spherical_to_cartesian(r, th, ph, add_front, add_back)
    
    def cartesian_grid_2D(self, plane, theta_points=None):
        """
        Returns a 2D cartesian grid.
        """
        if self.dim == 1:
            return self.__2D_cartesian_grid_1D(self.radius, theta_points,
                                               plane)
        elif self.dim == 2:
            return self.__2D_cartesian_grid(self.radius, self.theta, plane,
                                            theta_points)
        else:
            return self.__2D_cartesian_grid_from_3D(self.radius, self.theta,
                                                     self.phi, plane)
    
    def map_1D_to_2D(self, quantity, theta_points):
        """
        Maps a 1D quantity to a 2D grid.
        """
        assert theta_points is not None, "Theta angle MUST not be None."
        if self.dim == 1:
            return np.tile(quantity, (theta_points, 1))
        elif self.dim == 2:
            quantity = quantity[quantity.shape[0] // 2, :]
            return np.tile(quantity, (theta_points, 1))
        else:
            raise TypeError("Your simulation is 3D, you are not supposed to \
                            be here.")

    def gradient(self, quantity, concatenate=False):
        assert quantity.ndim == self.dim, "Dimension misatch."
        gradient = []
        gradient.append(IDL_derivative(self.radius, quantity, 'radius'))
        if self.dim > 1:
            gradient.append(1 / self.radius * IDL_derivative(self.theta, quantity, 'theta'))
        if self.dim > 2:
            gradient.append(1 / (self.radius * np.sin(self.theta)) * 
                                IDL_derivative(self.phi, quantity, 'phi'))
        if isinstance(quantity, aerray):
            symb = [r'$r$', r'$\theta$', r'$\phi$']
            for i in range(gradient):
                gradient[i].set(name=f'gradient_{quantity.name}_{symb[i]}',
                                label=merge_strings(r'$\nabla$', quantity.label, r'$\cdot$',
                                                    apply_symbol(symb[i], '\hat')),
                                log = quantity.log,
                                limits=quantity.limits,
                                cmap=quantity.cmap)
        if len(gradient) == 1:
            return gradient[0]
        if concatenate:
            np.concatenate(gradient, axis=0)
        else:
            return gradient

    def divergence(self, quantity, concatenate=False):
        if isinstance(quantity, list) or isinstance(quantity, tuple):
            if not len(quantity) == self.dim:
                raise ValueError("Dimension mismatch")
        elif quantity.ndim != self.dim + 1:
            raise ValueError("Dimension mismatch")
        else:
            raise NotImplementedError()
        if isinstance(aerray, quantity) or isinstance(quantity[0], aerray):
            isaerray = True
            
        else:
            False
        if (not (isinstance(quantity, list) or isinstance(quantity, tuple)) and 
            self.dim > 1):
            quantity = [quantity[i] for i in range(self.dim)]
        if self.dim == 1 and not isinstance(quantity, list):
            quantity = [quantity]
        if not isinstance(quantity, list):
            quantity = list(quantity)


    def __1D_cartesian_grid(self, radius):
        return radius

    def __2D_cartesian_grid_1D(self, radius, theta_points, plane='xz'):
        assert theta_points is not None, "Theta angle MUST not be None."
        if plane in ['xy', 'yx']:
            theta = np.linspace(0, 2*np.pi, theta_points, endpoint=True)
        else:
            theta = np.linspace(0, np.pi, theta_points, endpoint=True)
        X = radius[None, :] * np.sin(theta)[:, None]
        Y = radius[None, :] * np.cos(theta)[:, None]
        X.set(name='X', label=r'$X$', log=False, limits=[-1.5e7, 1.5e7])
        Y.set(name='Y', label=r'$Y$', log=False, limits=[-1.5e7, 1.5e7])
        return X, Y
    
    def __2D_cartesian_grid(self, radius, theta, plane, phi_points=None):
        if plane in ['xy', 'yx']:
            assert phi_points is not None, "Phi angle MUST NOT be None."
            theta = np.linspace(0, 2*np.pi, phi_points, endpoint=True)
        X = (radius[None, :] * np.sin(theta)[:, None])
        Y = (radius[None, :] * np.cos(theta)[:, None])
        if plane in ['xy', 'yx']:
            X.set(name='X', label=r'$X$', log=False, limits=[-1.5e7, 1.5e7])
            Y.set(name='Y', label=r'$Y$', log=False, limits=[-1.5e7, 1.5e7])
        else:
            X.set(name='X', label=r'$X$', log=False, limits=[0, 1.5e7])
            Y.set(name='Z', label=r'$Z$', log=False, limits=[-1.5e7, 1.5e7])
        return X, Y
    
    def __2D_cartesian_grid_from_3D(self, radius, theta, phi, plane):
        if plane in ['xy', 'yx']:
            X = radius[None, :] * np.sin(theta[len(theta) // 2]) * \
                np.cos(phi)[:, None]
            Y = radius[None, :] * np.sin(theta[len(theta) // 2]) * \
                np.sin(phi)[:, None]
            X.set(name='X', label=r'$X$', log=False, limits=[-1.5e7, 1.5e7])
            Y.set(name='Y', label=r'$Y$', log=False, limits=[-1.5e7, 1.5e7])
        elif plane in ['xz', 'zx']:
            theta = np.concatenate([theta, theta + (np.pi * u.radian)])
            X = radius[None, :] * np.sin(theta)[:, None]
            Y = radius[None, :] * np.cos(theta)[:, None]
            X.set(name='X', label=r'$X$', log=False, limits=[-1.5e7, 1.5e7])
            Y.set(name='Z', label=r'$Z$', log=False, limits=[-1.5e7, 1.5e7])
        elif plane in ['yz', 'zy']:
            theta = np.concatenate([theta, theta + (np.pi * u.radian)])
            X = radius[None, :] * np.sin(theta)[:, None]
            Y = radius[None, :] * np.cos(theta)[:, None]
            X.set(name='Y', label=r'$Y$', log=False, limits=[-1.5e7, 1.5e7])
            Y.set(name='Z', label=r'$Z$', log=False, limits=[-1.5e7, 1.5e7])
        elif plane in ['xz_phi_avg', 'zx_phi_avg', 'yz_phi_avg', 'zy_phi_avg']:
            X = radius[None, :] * np.sin(theta)[:, None]
            Y = radius[None, :] * np.cos(theta)[:, None]
            X.set(name='X', label=r'$X$', log=False, limits=[0, 1.5e7])
            Y.set(name='Z', label=r'$Z$', log=False, limits=[-1.5e7, 1.5e7])
        else:
            TypeError(f"Plane {plane} not available.")
        return X, Y

    def __3D_cartesian_grid(self, radius, theta, phi):
        X = radius[None, None, :] * np.sin(theta)[None, :, None] * \
            np.cos(phi)[:, None, None]
        Y = radius[None, None, :] * np.sin(theta)[None, :, None] * \
            np.sin(phi)[:, None, None]
        Z = radius[None, None, :] * np.cos(theta)[None, :, None] * \
            np.ones(phi.shape)[:, None, None]
        X.set(name='X', label=r'$X$', log=False, limits=[-1.5e7, 1.5e7])
        Y.set(name='Y', label=r'$Y$', log=False, limits=[-1.5e7, 1.5e7])
        Z.set(name='Z', label=r'$Z$', log=False, limits=[-1.5e7, 1.5e7])
        return X, Y, Z

    def __2D_spherical_to_cartesian(self, r, th, ph, add_front, add_back):
        sint = np.sin(self.theta)[:, None]
        cost = np.cos(self.theta)[:, None]
        for i in range(add_front):
            cost = cost[None, ...]
            sint = sint[None, ...]
        for i in range(add_back):
            cost = cost[..., None]
            sint = sint[..., None]
        X = r * sint + th * cost
        Z = r * cost - th * sint
        if ph is None:
            return X, Z
        Y = ph * sint
        return X, Y, Z
    
    def __3D_spherical_to_cartesian(self, r, th, ph, add_front, add_back):
        sint = np.sin(self.theta)[None, :, None]
        cost = np.cos(self.theta)[None, :, None]
        sinp = np.sin(self.phi)[:, None, None]
        cosp = np.cos(self.phi)[:, None, None]
        
        for i in range(add_front):
            cost = cost[None, ...]
            sint = sint[None, ...]
            cosp = cosp[None, ...]
            sinp = sinp[None, ...]
        for i in range(add_back):
            cost = cost[..., None]
            sint = sint[..., None]
            cosp = cosp[..., None]
            sinp = sinp[..., None]
        
        X = r * sint * cosp + \
            th * cost * cosp - \
            ph * sinp
        Y = r * sint * sinp + \
            th * cost * sinp + \
            ph * cosp
        Z = r * cost - \
            th * sint
        return X, Y, Z