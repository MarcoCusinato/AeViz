import numpy as np
from AeViz.units import aerray, u
from AeViz.utils.math_utils import IDL_derivative
from AeViz.utils.files.string_utils import merge_strings, apply_symbol
from copy import copy

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

    @property
    def shape(self):
        s = (len(self.radius),)
        if self.dim > 1:
            s = (len(self.theta),) + s
        if self.dim > 2:
            s = (len(self.phi),) + s
        return s

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
            raise TypeError("Grid built from raw coordinates. No cell support.")
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

    def gradient(self, quantity, stack=False):
        assert quantity.ndim == self.dim, "Dimension misatch."
        gradient = []
        gradient.append(IDL_derivative(self.radius, quantity, 'radius'))
        if self.dim > 1:
            gradient.append(1 / self.radius * IDL_derivative(self.theta, quantity, 'theta'))
        if self.dim > 2:
            gradient.append(1 / (self.radius * np.sin(self.theta)[None, :, None]) * 
                                IDL_derivative(self.phi, quantity, 'phi'))
        if isinstance(quantity, aerray):
            symb = [r'$r$', r'$\theta$', r'$\phi$']
            for i in range(len(gradient)):
                gradient[i].set(name=f'gradient_{quantity.name}_{symb[i]}',
                                label=merge_strings(r'$\nabla$', quantity.label, r'$\cdot$',
                                                    apply_symbol(symb[i], r'\hat')),
                                log = quantity.log,
                                limits=quantity.limits,
                                cmap=quantity.cmap)
        elif isinstance(gradient[0], aerray):
            gradient = [gg.value for gg in gradient]
        if len(gradient) == 1:
            return gradient[0]
        if stack:
            return np.stack(gradient, axis=-1)
        else:
            return gradient

    def divergence(self, quantity):
        if isinstance(quantity, list) or isinstance(quantity, tuple):
            if not len(quantity) == self.dim:
                raise ValueError("Dimension mismatch")
            if isinstance(quantity, tuple):
                quantity = list(quantity)
        elif (quantity.ndim != self.dim + 1) or self.dim == 1:
            raise ValueError("Dimension mismatch")
        
        if isinstance(quantity, aerray) or isinstance(quantity[0], aerray):
            isaerray = True
            if isinstance(quantity, aerray):
                qt = quantity
            else:
                qt = quantity[0]
            name = f'divergence_{qt.name}'
            label = merge_strings(r'$\nabla\cdot$', qt.label.replace('_r', ''))
            limits = qt.limits
            cmap = qt.cmap
            log = qt.log
        else:
            isaerray = False
        if not (isinstance(quantity, list) and self.dim > 1):
            quantity = [quantity[..., i] for i in range(3)]
        if self.dim == 1 and not isinstance(quantity, list):
            quantity = [quantity]
        div = 1 / self.radius ** 2 * IDL_derivative(self.radius, quantity[0] * self.radius ** 2,
                                                    'radius')
        if self.dim > 1:
            st = np.sin(self.theta)
            if self.dim == 2:
                st = st[:, None]
            else:
                st = st[None, :, None]
            div += 1 / (self.radius * st) * IDL_derivative(self.theta,
                                                           quantity[1] * st,
                                                           'theta')
        if self.dim > 2:
            div += 1 / (self.radius * st) * IDL_derivative(self.phi, quantity[2],
                                                           'phi')
        if isaerray:
            div.set(name=name, label=label, limits=limits, cmap=cmap, log=log)
        elif isinstance(div, aerray):
            div = div.value
        return div

    def curl(self, quantity, stack=False):
        if self.dim == 1:
            if isinstance(quantity, aerray):
                return np.zeros(len(quantity)) * quantity.unit / self.radius.unit
            else:
                return np.zeros(len(quantity))
        if isinstance(quantity, list) or isinstance(quantity, tuple):
            if not len(quantity) == self.dim:
                raise ValueError("Dimension mismatch")
            if isinstance(quantity, tuple):
                quantity = list(quantity)
        elif quantity.ndim != self.dim + 1:
            raise ValueError("Dimension mismatch")
        
        if isinstance(quantity, aerray) or isinstance(quantity[0], aerray):
            isaerray = True
            if isinstance(quantity, aerray):
                qt = quantity
            else:
                qt = quantity[0]
            name = f'curl_{qt.name}'
            label = merge_strings(r'$(\nabla\times$', qt.label.replace('_r', ''), r'$)$')
            limits = qt.limits
            cmap = qt.cmap
            log = qt.log
        else:
            isaerray = False
        if not (isinstance(quantity, list)):
            quantity = [quantity[..., i] for i in range(3)]
        st = np.sin(self.theta)
        if self.dim == 2:
            st = st[:, None]
        else:
            st = st[None, :, None]
        crl = [
            1 / (self.radius * st) * (IDL_derivative(self.theta, quantity[2] * st,
                                                     'theta')),
            - 1 / self.radius * IDL_derivative(self.radius, quantity[2] * self.radius,
                                               'radius'),
            1 / self.radius * (IDL_derivative(self.radius, quantity[1] * self.radius,
                                              'radius') -
                               IDL_derivative(self.theta, quantity[0], 'theta'))
        ]
        if self.dim > 2:
            crl[0] -= 1 / (self.radius * st) * IDL_derivative(self.phi, quantity[1], 'phi')
            crl[1] += 1 / (self.radius * st) * IDL_derivative(self.phi, quantity[0], 'phi')
        if isaerray:
            symb = [r'$r$', r'$\theta$', r'$\phi$']
            for i in range(len(crl)):
                crl[i].set(
                    name = f'{name}_{symb[i]}',
                    limits = limits,
                    cmap = cmap,
                    log = log,
                    label = merge_strings(label, apply_symbol(symb[i], r"\hat"))
                )
        elif isinstance(crl[0], aerray):
            crl = [cr.value for cr in crl]
        if stack:
            return np.stack(crl, axis=-1)
        return crl

    def volume_average(self, quantity):
        dV = self._metric().dVolume_integration(self.__ghost)
        mask = np.isfinite(quantity)
        av = np.sum(quantity[mask] * dV[mask]) / np.sum(dV[mask])
        if isinstance(quantity, aerray):
            av.set(name=merge_strings(quantity.name, '_vol_avg'), 
                   label=merge_strings(r'$\langle $', quantity.label, 
                                            r'$\rangle_V$'),
                   log=quantity.log, limits=quantity.limits)
        elif isinstance(av, aerray):
            av = av.value
        return av

    def omega_average(self, quantity):
        if self.dim < 2:
            return quantity
        dOmega = self._metric().dOmega(self.__ghost)
        dOmega = dOmega[..., None]
        mask = np.isfinite(quantity)
        av = np.nansum(quantity * dOmega * mask, axis=tuple(range(self.dim-1))) / \
            np.sum(dOmega * mask, axis=tuple(range(self.dim-1)))
        if isinstance(quantity, aerray):
            av.set(name=merge_strings(quantity.name, '_ang_avg'), 
                   label=merge_strings(r'$\langle $', quantity.label, 
                                        r'$\rangle_\Omega$'),
                   log=quantity.log, limits=quantity.limits)
        elif isinstance(av, aerray):
            av = av.value
        return av

    def project(self, quantity, axis):
        def axis_vector(axis):
            cart_axis = {'x': (1., 0., 0.),
                         'y': (0., 1., 0.),
                         'z': (0., 0., 1.)}
            if isinstance(axis, str):
                if axis.lower() not in cart_axis:
                    raise ValueError(f"Axis '{axis}' not valid, use 'x', 'y', 'z' "
                                    "or a tuple (nx, ny, nz).")
                return np.array(cart_axis[axis.lower()])
            n = np.asarray(axis, dtype=float)
            if n.shape != (3,):
                raise ValueError("A tuple axis must have 3 cartesian components.")
            norm = np.linalg.norm(n)
            if norm == 0:
                raise ValueError("The axis vector cannot be null.")
            return n / norm
        
        if self.dim == 1:
            return quantity
        n = axis_vector(axis)
        components = self.__vector_components(quantity)
        ccart = self.spherical_to_cartesian(components[0], components[1], components[2])
        proj = n[0] * ccart[0] + n[1] * ccart[1] + n[2] * ccart[3]
        if isinstance(components[0], aerray):
            tag = axis.lower() if isinstance(axis, str) else 'n'
            proj.set(name=f"{components[0].name}_proj_{tag}",
                    label=merge_strings(components[0].label.replace('_r', '') if components[0].label is not None else '',
                                        r'$\cdot$',
                                        apply_symbol(f'${tag}$', r'\hat')),
                    limits=components[0].limits, cmap=components[0].cmap,
                    log=False)       # projections can change sign
        elif isinstance(proj, aerray):
            proj = proj.value
        return proj

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
            raise TypeError(f"Plane {plane} not available.")
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

    def __vector_components(self, quantity):
        """List of components from a list/tuple or an array with components
        on the last axis. Checks count and shape against the grid."""
        if isinstance(quantity, (list, tuple)):
            comps = list(quantity)
        else:
            if quantity.ndim != self.dim + 1 and self.dim > 1:
                raise ValueError("Dimension mismatch")
            if self.dim == 1:
                if not isinstance(quantity, list):
                    comps = [quantity]
                else:
                    comps = quantity
            else:
                if quantity.shape[-1] != 3:
                    raise ValueError(f"Expected at 3 components, "
                                                f"got {quantity.shape[-1]}.")
                comps = [quantity[..., i] for i in range(3)]            
        for c in comps:
            if c.shape != self.shape:
                raise ValueError(f"Component shape {c.shape} does not match the "
                                f"grid {self.shape}. Was it trimmed with the "
                                "same ghost cells?")
        return comps