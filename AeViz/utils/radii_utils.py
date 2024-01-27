from AeViz.utils.file_utils import save_hdf
import numpy as np
from AeViz.utils.math_utils import IDL_derivative



def PNS_radius(simulation, file_name):
    """
    Calculates the radius of the PNS for each timestep.
    Employed method: radius at which the density drops below
                     10^11 g/cm^3
    """
    return simulation.cell.radius(simulation.ghost)[np.argmax(
        simulation.rho(file_name) < 1e11, axis=-1)]

def innercore_radius(simulation, file_name):
    """
    Calculates the radius of the inner core for each timestep. 
    The innercore is defined as the region in sonic contact with the
    centre of the star. This corresponds to the region where the
    sound speed is larger than the radial velocity.
    Employed method: comparing the module of the radial with the
                     sound speed.
    """
    return simulation.cell.radius(simulation.ghost)[np.argmax(
        simulation.radial_velocity(file_name) ** 2 
        >= simulation.soundspeed(file_name) ** 2, axis=-1)]

def gain_radius(simulation, file_name, PNS_radius):
    """
    Calculates the radius of the gain region for each timestep. Gain
    region is defined as the region where the neutrino heating becomes
    larger than the cooling.
    Employed method: finding the radius, outside the PNS where the 
                     neutrino energy depostion is larger than 0.
    """
    radius = simulation.cell.radius(simulation.ghost)
    neu_ene_dep = simulation.nu_heat(file_name)
    while radius.ndim < neu_ene_dep.ndim:
        radius = radius[None, ...]
    neu_ene_dep = np.where(radius >= PNS_radius[..., None], neu_ene_dep,  0)
    return simulation.cell.radius(simulation.ghost)[np.argmax(neu_ene_dep > 0,
                                                              axis=-1)]

def neutrino_sphere_radii(simulation, file_name):
    """
    Calculates the neutrino sphere for each neutrino flvaour radii at
    each timestep.
    Employed method: calculate k(r, θ, ϕ) as
    k(r, θ, ϕ) = (∑_ω p^α(r, θ, ϕ, ω) κ_α(r, θ, ϕ, ω)) / \
        ∑_ω p^α(r, θ, ϕ, ω)
    where p^α(r, θ, ϕ, ω) is the neutrino three-momentum of each 
    neutrino flavour,  ω is the neutrino energy and κ_α(r, θ, ϕ, ω) is
    the neutrino opacity.
    Then the neutrino sphere is defined as the radius where the integral
    τ = ∫^R_∞ dr k(r, θ, ϕ) is equal to 1.
    """
    tau = 1
    momenta = simulation.neutrino_momenta(file_name)
    kappa = simulation.neutrino_momenta_opacities(file_name)
    k = np.sum(momenta * kappa, axis=(-1, simulation.dim)) / \
               np.sum(momenta, axis=(-1, simulation.dim))
    dr = simulation.cell.dr(simulation.ghost)[..., None]
    while k.ndim < dr.ndim:
        k = k[..., None]
    k = np.nancumsum(np.flip(k * dr, axis=-2), axis=-2)
    return np.flip(simulation.cell.radius(simulation.ghost))\
                                [np.argmax(k >= tau, axis=-2)]

def PNS_nucleus(simulation, file_name):
    """
    Calculates the radius of the PNS nucleus for each timestep.
    Employed method: point at minimum Ye inside the PNS. If such point
                        does not exist it returns the point at which the
                        density is greater then 5x10^{12} g/cm³.
    """
    radius = simulation.cell.radius(simulation.ghost)
    R_30Km_index = np.argmax(radius >= 30e5)
    rho_thres_index = np.argmax(simulation.rho(file_name) <= 5e12, axis=-1)
    Ye = simulation.Ye(file_name)[..., :R_30Km_index]
    dYedr = IDL_derivative(radius[:R_30Km_index], Ye)
    minYe = np.argmax(Ye, axis=-1)
    mask =  (dYedr[np.arange(len(minYe)), minYe-1] < 0) & \
        (dYedr[np.arange(len(minYe)), minYe+1] > 0)
    minYe = np.where(mask, minYe, rho_thres_index)
    return radius[minYe]

def shock_radius(simulation, file_name):
    pass