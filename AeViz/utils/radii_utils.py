import numpy as np
from AeViz.utils.math_utils import IDL_derivative
from scipy.interpolate import griddata

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
    k = np.nansum(momenta * kappa, axis=(-1, simulation.dim)) / \
               np.nansum(momenta, axis=(-1, simulation.dim))
    dr = simulation.cell.dr(simulation.ghost)[..., None]
    while dr.ndim < k.ndim:
        dr = dr[None, ...]
    k = np.nancumsum(np.flip(k * dr, axis=-2), axis=-2)
    return np.flip(simulation.cell.radius(simulation.ghost))\
                                [np.argmax(k >= tau, axis=-2)]

def PNS_nucleus(simulation, file_name):
    """
    Calculates the radius of the PNS nucleus for each timestep.
    Employed method: entropy jump at s=4kb from the inside out.
    """
    radius = simulation.cell.radius(simulation.ghost)
    R_30Km_index = np.argmax(radius >= 30e5)
    s = simulation.entropy(file_name)[..., :R_30Km_index]
    minYe = np.argmax(s >= 4, axis=-1)
    return radius[minYe]



def shock_radius(simulation, file_name):
    """
    Calculates the shock radius for each timestep.
    Employed method: first jump in pressure and radial velocity after
                    the bounce, considered from infinite to the centre.
    """
    if simulation.time(file_name, True) <= 0:
        if simulation.dim == 1:
            return np.array([0])
        return np.zeros(simulation.cell.dVolume_integration(
            simulation.ghost).shape[:-1])
    if simulation.dim == 1:
        return shock_radius_1D(simulation, file_name)
    elif simulation.dim == 2:
        return interpol_1D(hampel_filter(
            shock_radius_2D(simulation, file_name)),
                           simulation.cell.theta(simulation.ghost))
    elif simulation.dim == 3:
        Theta, Phi = np.meshgrid(simulation.cell.theta(simulation.ghost), 
                                 simulation.cell.phi(simulation.ghost))
        return interpol_2D(shock_radius_3D(simulation,
                                                         file_name),
                           Theta, Phi)
    else:
        raise ValueError("Invalid dimension")
    
    
def shock_radius_1D(simulation, file_name):
    dP = IDL_derivative(simulation.cell.radius(simulation.ghost),
                        simulation.gas_pressure(file_name)) * \
                            simulation.cell.radius(simulation.ghost) / \
                            simulation.gas_pressure(file_name)
    dvr = IDL_derivative(simulation.cell.radius(simulation.ghost),
                         simulation.radial_velocity(file_name)) * \
                             simulation.cell.radius(simulation.ghost) / \
                             np.abs(simulation.radial_velocity(file_name))
    for ir in range(len(dP) - 1):
        if (dP[ir] < 10) and np.any(dvr[ir-5:ir+6] < -20):
            return np.array([simulation.cell.radius(simulation.ghost)[ir]])
    return np.array([0])
            

def shock_radius_2D(simulation, file_name):
    dP = IDL_derivative(simulation.cell.radius(simulation.ghost),
                        simulation.gas_pressure(file_name)) * \
                            simulation.cell.radius(simulation.ghost) / \
                            simulation.gas_pressure(file_name)
    dvr = IDL_derivative(simulation.cell.radius(simulation.ghost),
                         simulation.radial_velocity(file_name)) * \
                             simulation.cell.radius(simulation.ghost) / \
                             np.abs(simulation.radial_velocity(file_name))
    s = simulation.entropy(file_name)
    shock_r = np.empty(dP.shape[0])
    shock_r.fill(np.nan)
    for it in range(dP.shape[0]):
        for ir in reversed(range(dP.shape[1] - 1)):
            if (dP[it, ir] < -10) and \
                (np.any(dvr[it, max(0,ir-5):min(ir+6, dP.shape[1]-1)] < -20)):
                shock_r[it] = simulation.cell.radius(simulation.ghost)[ir]
                break
    ## COPY over the gcells
    if np.isnan(shock_r).all():
        return np.zeros(dP.shape[0])
    return shock_r

def shock_radius_3D(simulation, file_name):
    """
    Copied from Martin's IDL script.
    """
    dP = IDL_derivative(simulation.cell.radius(simulation.ghost),
                        simulation.gas_pressure(file_name)) * np.abs(
                            simulation.cell.radius(simulation.ghost) / \
                            simulation.gas_pressure(file_name))
    dvr = IDL_derivative(simulation.cell.radius(simulation.ghost),
                         simulation.radial_velocity(file_name)) * \
                             simulation.cell.radius(simulation.ghost) / \
                             np.abs(simulation.soundspeed(file_name))
    ds = IDL_derivative(simulation.cell.radius(simulation.ghost),
                        simulation.entropy(file_name)) * np.abs(
                            simulation.cell.radius(simulation.ghost) / \
                            simulation.entropy(file_name))
                        
    vr = simulation.radial_velocity(file_name)
    shock_r = np.empty((dP.shape[0], dP.shape[1]))
    shock_r.fill(np.nan)
    for ip in range(dP.shape[0]):
        for it in range(dP.shape[1]):
            for ir in range(dP.shape[2] - 1):
                if (ds[ip, it, ir] < -0.25) and \
                    (vr[ip, it, ir] >= 5e-7) and \
                    (dvr[ip, it, ir] <= -1) and \
                    (dP[ip, it, ir] <= -1) and \
                    (vr[ip, it, max(0, ir-10)] >= vr[ip, it, min(vr.shape[2],
                                                                 ir+10)]):
                    shock_r[ip, it] = simulation.cell.radius(simulation.ghost)[ir]
                    break
    return shock_r
    
def shock_radius_3D_OLD(simulation, file_name):
    dP = IDL_derivative(simulation.cell.radius(simulation.ghost),
                        simulation.gas_pressure(file_name)) * \
                            simulation.cell.radius(simulation.ghost) / \
                            simulation.gas_pressure(file_name)
    dvr = IDL_derivative(simulation.cell.radius(simulation.ghost),
                         simulation.radial_velocity(file_name)) * \
                             simulation.cell.radius(simulation.ghost) / \
                             np.abs(simulation.radial_velocity(file_name))
    s = simulation.entropy(file_name)
    shock_r = np.empty((dP.shape[0], dP.shape[1]))
    shock_r.fill(np.nan)
    for ip in range(dP.shape[0]):
        for it in range(dP.shape[1]):
            for ir in range(dP.shape[2] - 1):
                if (dP[ip, it, ir] < -500) and \
                (np.any(dvr[ip, it, max(0,ir-5):min(ir+6, dP.shape[2] - 1)] 
                    < -20)) and \
                (np.all(s[ip, it, max(0,ir-5):min(ir+6, dP.shape[2] - 1)]
                        < 100)):
                    shock_r[ip, it] = simulation.cell.radius(simulation.ghost)[ir]
                    break
    return shock_r

def hampel_filter(shock_radius, sigma=3):
    """
    Applies a hampel filter to the shock radius.
    """
    
    assert shock_radius.size > 1 and shock_radius.ndim >= 1, "Expected 1 or 2 \
            dimensional array, with at least 2 elements."
    for i in range(10):
        rmedian = np.nanmedian(shock_radius)
        diff = np.abs(shock_radius - rmedian)
        diffmedian = np.nanmedian(diff)
        if i == 0:
            threshold = np.minimum(sigma, np.nanmax(diff / (1.4826 * \
                diffmedian)) * 0.95)
        mask = diff > 1.4826 * diffmedian * threshold
        shock_radius[mask] = np.nan
        num_outliers = np.sum(mask)
        if num_outliers == 0:
            break
    return shock_radius
        
    
def interpol_1D(shock_radius, theta):
    mask = np.isnan(shock_radius)
    if mask.sum() == 0:
        return shock_radius
    shock_radius[mask] = np.interp(theta[mask], theta[~mask],
                                   shock_radius[~mask])
    return shock_radius

def interpol_2D(shock_radius, Theta, Phi):
    mask = np.isnan(shock_radius)
    if mask.sum() == 0:
        return shock_radius
    median = np.nanmedian(shock_radius)
    shock_radius = np.ma.masked_invalid(shock_radius)
    shock_radius[shock_radius.mask] = \
        griddata((Phi[~shock_radius.mask], Theta[~shock_radius.mask]),
                 shock_radius[~shock_radius.mask].ravel(),
                 (Phi[shock_radius.mask], Theta[shock_radius.mask]),
                    method='linear', fill_value=median)
    return shock_radius
    
