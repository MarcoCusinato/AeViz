import numpy as np
from AeViz.utils.math_utils import IDL_derivative
from scipy.interpolate import griddata
from AeViz.units import u

def PNS_radius(simulation, file_name):
    """
    Calculates the radius of the PNS for each timestep.
    Employed method: radius at which the density drops below
                     10^11 g/cm^3
    """
    return simulation.cell.radius(simulation.ghost)[np.argmax(
        simulation.rho(file_name) < (1e11 * u.g / u.cm ** 3), axis=-1)]

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
    if simulation.dim==1:
        axis = -1
    else:
        axis = (-1, simulation.dim)
    k = [np.nansum(mom * ka, axis=axis) /
          np.nansum(mom, axis=axis) for (mom, ka) in
          zip(momenta, kappa)] 
    dr = simulation.cell.dr(simulation.ghost)
    while dr.ndim < k[0].ndim:
        dr = dr[None, ...]
    k = [np.nancumsum(np.flip(ka * dr, axis=-1), axis=-1) for ka in k]
    return [np.flip(simulation.cell.radius(simulation.ghost))\
                                [np.argmax(ka >= tau, axis=-1)] for ka in k]

def PNS_nucleus(simulation, file_name):
    """
    Calculates the radius of the PNS nucleus for each timestep.
    Employed method: entropy jump at s=4kb from the inside out.
    """
    radius = simulation.cell.radius(simulation.ghost)
    R_30Km_index = np.argmax(radius >= (30 * u.km))
    s = simulation.entropy(file_name)[..., :R_30Km_index]
    S4Kb = np.argmax(s >= (4 * u.kBol / u.bry), axis=-1)
    return radius[S4Kb]

def shock_radius(simulation, file_name, rmax=None):
    """
    Calculates the shock radius for each timestep.
    Employed method: first jump in pressure and radial velocity after
                    the bounce, considered from infinite to the centre.
    """
    if rmax is None:
        rmax = simulation.cell.radius(simulation.ghost)[-1]
    if simulation.time(file_name, True) <= 0:
        if simulation.dim == 1:
            return 0.0 * simulation.cell.radius(simulation.ghost).unit
        return np.zeros(simulation.cell.dVolume_integration(
            simulation.ghost).shape[:-1]) * \
                simulation.cell.radius(simulation.ghost).unit
    if simulation.dim == 1:
        return shock_radius_1D(simulation, file_name, rmax)
    elif simulation.dim == 2:
        return interpol_1D(hampel_filter(
            shock_radius_2D(simulation, file_name, rmax)),
                           simulation.cell.theta(simulation.ghost))
    elif simulation.dim == 3:
        Theta, Phi = np.meshgrid(simulation.cell.theta(simulation.ghost), 
                                 simulation.cell.phi(simulation.ghost))
        return interpol_2D(shock_radius_3D(simulation, file_name, rmax),
                           Theta, Phi)
    else:
        raise ValueError("Invalid dimension")
   
def shock_radius_1D(simulation, file_name, rmax):
    r = simulation.cell.radius(simulation.ghost)
    vr = simulation.radial_velocity(file_name)
    s = simulation.entropy(file_name)
    dP = IDL_derivative(r, s) * r / s
    dvr = IDL_derivative(r, vr) * r / np.abs(vr)
    for ir in range(len(dP) - 1):
        if r[ir] > rmax:
            continue
        if (dP[ir] < -5) and np.any(dvr[ir-5:ir+6] < -20):
            return r[ir]
    return 0.0 * r.unit

def shock_radius_2D(simulation, file_name, rmax):
    vr = simulation.radial_velocity(file_name)
    r = simulation.cell.radius(simulation.ghost)
    p = simulation.gas_pressure(file_name)
    dP = IDL_derivative(r, p) * r / p
    dvr = IDL_derivative(r, vr) * r / np.abs(vr)
    s = simulation.entropy(file_name)
    shock_r = np.empty(dP.shape[0])
    shock_r.fill(np.nan)
    for it in range(dP.shape[0]):
        for ir in reversed(range(dP.shape[1] - 1)):
            if r[ir] > rmax:
                continue
            if (dP[it, ir] < -10) and \
                (np.any(dvr[it, max(0,ir-5):min(ir+6, dP.shape[1]-1)] < -20)) \
                and (np.abs(vr[it, ir]) > 1e8) and s[it, ir] < 400 :
                shock_r[it] = r[ir]
                break
    ## COPY over the gcells
    if np.isnan(shock_r).all():
        return np.zeros(dP.shape[0]) * r.unit
    return shock_r * r.unit

def shock_radius_3D(simulation, file_name, rmax):
    """
    Copied from Martin's IDL script.
    """
    r = simulation.cell.radius(simulation.ghost)
    p = simulation.gas_pressure(file_name)
    vr = simulation.radial_velocity(file_name)
    entr = simulation.entropy(file_name)
    dP = IDL_derivative(r, p) * np.abs(r / p)
    dvr = IDL_derivative(r, vr) * r / np.abs(simulation.soundspeed(file_name))
    ds = IDL_derivative(r, entr) * np.abs(r / entr)
    shock_r = np.empty((dP.shape[0], dP.shape[1]))
    shock_r.fill(np.nan)
    for ip in range(dP.shape[0]):
        for it in range(dP.shape[1]):
            for ir in range(dP.shape[2] - 1):
                if r[ir] > rmax:
                    continue
                if (ds[ip, it, ir] < -0.15) and \
                    (vr[ip, it, ir] >= 1) and \
                    (dvr[ip, it, ir] <= -0.7) and \
                    (dP[ip, it, ir] <= -0.7) and \
                    (vr[ip, it, max(0, ir-10)] >= vr[ip, it, min(vr.shape[2],
                                                                 ir+10)]):
                    shock_r[ip, it] = r[ir]
                    break
    return shock_r * r.unit
    
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
    shock_radius_copy = shock_radius.copy()
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
    if np.sum(np.isnan(shock_radius)) >= shock_radius.size / 2:
        return shock_radius_copy
    return shock_radius
 
def interpol_1D(shock_radius, theta):
    mask = np.isnan(shock_radius)
    if mask.sum() == 0:
        return shock_radius
    shock_radius[mask] = np.interp(theta[mask], theta[~mask],
                                   shock_radius[~mask])
    return shock_radius

def interpol_2D(shock_radius, Theta, Phi):
    unit = shock_radius.unit
    shock_radius = shock_radius.value
    mask = np.isnan(shock_radius)
    if mask.sum() == 0:
        return shock_radius * unit
    if np.all(mask):
        return np.zeros_like(shock_radius) * unit
    median = np.nanmedian(shock_radius)
    shock_radius = np.ma.masked_invalid(shock_radius)
    shock_radius[shock_radius.mask] = \
        griddata((Phi.value[~shock_radius.mask], Theta.value[~shock_radius.mask]),
                 shock_radius[~shock_radius.mask].ravel(),
                 (Phi.value[shock_radius.mask], Theta.value[shock_radius.mask]),
                    method='linear', fill_value=median)
    return shock_radius * unit
    
