from AeViz.simulation.methods import *
from AeViz.utils.math_utils import IDL_derivative, function_average

"""
Function to process convection and instabilities data from a simulation
in spherical coordinates.
These functions are not meant to be used standalone, but rather to be
imported into the Simulation class.
"""

## -----------------------------------------------------------------
## CONVECTION AND INSTABILITIES DATA
## -----------------------------------------------------------------

@smooth
def BV_frequency(self, file_name, mode=1, **kwargs):
    """
    Returns the Brunt-Vaisala frequency at specific timestep.
    """
    rho = self.rho(file_name, **kwargs)
    radius = self.cell.radius(self.ghost)
    BV = (1 / self.soundspeed(file_name, **kwargs) ** 2 * \
            IDL_derivative(radius, self.gas_pressure(file_name, **kwargs)) - \
            IDL_derivative(radius, rho)) / rho
    if mode == 1:
        """
        check e.g. Gossan+20 `10.1093/mnras/stz3243`
        """
        geff = IDL_derivative(radius,
                                self.gravitational_potential(file_name,
                                                            **kwargs))
    elif mode == 2:
        """
        Check Fryer+21 ` 10.1134/S1063772921100103`
        """
        
        vr = self.radial_velocity(file_name, **kwargs)
        geff = IDL_derivative(radius,
                                self.gravitational_potential(file_name,
                                                            **kwargs)) - \
                vr * IDL_derivative(radius, vr) 
    else:
        raise ValueError("Mode not recognized.")
    return geff * BV

@smooth
def convective_velocity(self, file_name, **kwargs):
    """
    Returns the convective velocity at specific timestep. Defined as
    in `https://doi.org/10.3847/1538-4357/ac4507`:
    v_conv = <vr-vr_ave>_omega
    """
    dOmega = self.cell.dOmega(self.ghost)
    rho = self.rho(file_name, **kwargs)
    vr = self.radial_velocity(file_name, **kwargs)
    vrave = function_average(vr * rho, self.dim, 'Omega', dOmega) / \
        function_average(rho, self.dim, 'Omega', dOmega)
    return function_average((vr - vrave), self.dim, 'Omega', dOmega)

@smooth
def turbulent_velocity(self, file_name, **kwargs):
    """
    Returns the turbulent velocity at specific timestep. Defined as
    in `https://doi.org/10.3847/1538-4357/ac4507`:
    v_conv = <(v-v_ave)²>^0.5_omega
    """
    dOmega = self.cell.dOmega(self.ghost)
    rho = self.rho(file_name, **kwargs)
    rho_ave = function_average(rho, self.dim, 'Omega', dOmega)
    vr, vtheta, vphi = self.radial_velocity(file_name, **kwargs), \
        self.theta_velocity(file_name, **kwargs), self.phi_velocity(file_name, **kwargs)
    vrave, vthetaave, vphiave = \
        function_average(vr * rho, self.dim, 'Omega', dOmega) / rho_ave, \
            0, function_average(vphi, self.dim, 'Omega', dOmega) / rho_ave
    return function_average((vr - vrave) ** 2 + (vtheta - vthetaave) ** \
        2 + (vphi - vphiave) ** 2, self.dim, 'Omega', dOmega) ** 0.5

@smooth
def convective_flux(self, file_name, **kwargs):
    """
    Returns the convective flux at specific timestep. Defined as
    in `https://doi.org/10.3847/1538-4357/ac4507`:
    F_conv = <(0.5 rho v_turb² + e + P)v_conv>_omega
    """
    return function_average((0.5 * self.rho(file_name, **kwargs) * \
        self.turbulent_velocity(file_name, **kwargs) ** 2 + self.internal_energy(
            file_name, **kwargs) + self.gas_pressure(file_name, **kwargs)) * \
        self.convective_velocity(file_name, **kwargs), self.dim, 'Omega', 
        self.cell.dOmega(self.ghost))

@smooth
def Rossby_number(self, file_name, lenghtscale=True):
    """
    Returns the Rossby number at specific timestep. Defined as
    in `https://doi.org/10.3847/1538-4357/ac4507`:
    Ro = v_conv / (Omega R)
    If lenghtscale is True, the Rossby number is divided by a sort
    of typical lenghtscale of the convection.
    H = 1/|∂_r ρ/ρ|
    """
    if lenghtscale:
        H = 1 / np.abs(IDL_derivative(self.cell.radius(self.ghost), self.rho(
            file_name)) / self.rho(file_name))
    else:
        H = 1
    return self.convective_velocity(file_name) / (self.cell.radius(
        self.ghost) * self.omega(file_name) * H)

@smooth
def epicyclic_frequency(self, file_name, **kwargs):
    """
    Returns the epicyclic frequency at specific timestep. Defined as
    kappa² = ( (2Ω) /  R ) ∂(R²Ω)/∂R
    R is the cylindrical radius
    """
    if self.dim == 1:
        return None
    # get the spherical coordinates
    r = self.cell.radius(self.ghost)
    theta = self.cell.theta(self.ghost)
    # get the rotational frequency
    omega = self.omega(file_name)
    # We need to derive omega by r and theta
    domgdr = IDL_derivative(r, omega)
    domgdtheta = IDL_derivative(theta, omega, 'theta')
    # fix r and theta dimensions according to the simulation dimension
    if self.dim == 2:
        r = r[None, :]
        theta = theta[:, None]
    elif self.dim == 3:
        r = r[None, None, :]
        theta = theta[None, :, None]
    else:
        raise ValueError("Dimension not recognized.")
    # get the cylindrical radius        
    R = r * np.sin(theta)
    domgdr /= np.sin(theta)
    domgdtheta /= (r * np.cos(theta))
    kappa = (2 * omega ) / R * (R ** 2 * (domgdr + domgdtheta) + \
                                2 * R * omega)
    return kappa
