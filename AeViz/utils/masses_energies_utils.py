from AeViz.utils.file_utils import save_hdf
from AeViz.utils.utils import progressBar, check_existence
from AeViz.units.units import units
import numpy as np
from AeViz.utils.radii_utils import gain_radius

u = units()

def innercore_mass_energy(simulation, file_name, innercore_radius, gcells, dV):
    """
    Calculates the mass, kinetic, magnetic, rotational, gravitational
    and total energy and T/|W| of the inner core for one timestep.
    """
    if simulation.dim == 1:
        mask = (simulation.cell.radius(simulation.ghost) <= \
            innercore_radius)
    else:
        mask = (simulation.cell.radius(simulation.ghost) <= \
                simulation.ghost.remove_ghost_cells_radii(innercore_radius, 
                                simulation.dim, **gcells)[..., None])
    rho = simulation.rho(file_name)[mask] * dV[mask]
    if simulation.dim == 1:
        ene_rot = 0
        ene_kin = np.sum(0.5 * rho * \
            simulation.radial_velocity(file_name)[mask] ** 2)
        ene_mag = 0
        
    else:
        ene_rot = np.sum(0.5 * rho * \
            simulation.phi_velocity(file_name)[mask] ** 2)
        ene_kin = np.sum(0.5 * rho * \
            (simulation.radial_velocity(file_name)[mask] + \
            simulation.theta_velocity(file_name)[mask]) ** 2)
        ene_mag = np.sum(simulation.magnetic_energy(file_name)[mask] \
            * dV[mask])
    ene_grav = np.sum(simulation.gravitational_energy(file_name)[mask] \
        * dV[mask])
    return u.convert_to_solar_masses(np.sum(rho)), ene_kin, ene_mag, ene_rot, \
        ene_grav, ene_kin + ene_rot, ene_kin / np.abs(ene_grav)
        
def  gain_region_mass_energy(simulation, file_name, shock_radius, sgcells, 
                             gain_radius, ggcells, dV):
    """
    Calculates the mass, neutrino heating energy of the gain region for
    one timestep.
    """
    if simulation.dim == 1:
        mask = (simulation.cell.radius(simulation.ghost) >= \
            np.minimum(shock_radius, gain_radius)) & \
                (simulation.cell.radius(simulation.ghost) <= \
            np.maximum(shock_radius, gain_radius))
    else:
        mask = (simulation.cell.radius(simulation.ghost) >= \
            np.minimum(simulation.ghost.remove_ghost_cells_radii(shock_radius, 
                                simulation.dim, **sgcells)[..., None], 
                       simulation.ghost.remove_ghost_cells_radii(gain_radius, 
                                simulation.dim, **ggcells)[..., None])) & \
                (simulation.cell.radius(simulation.ghost) <= \
            np.maximum(simulation.ghost.remove_ghost_cells_radii(shock_radius, 
                                simulation.dim, **sgcells)[..., None], 
                       simulation.ghost.remove_ghost_cells_radii(gain_radius, 
                                simulation.dim, **ggcells)[..., None]))
    mass = u.convert_to_solar_masses(np.sum(simulation.rho(file_name)[mask] \
        * dV[mask]))
    nu_heat = np.sum(simulation.nu_heat(file_name)[mask] * dV[mask])
    return mass, nu_heat

def PNS_mass_energy(simulation, file_name, PNS_radius, gcells, dV):
    """
    Calculates the mass, kinetic, magnetic, rotational, gravitational
    and convective energy of the PNS for one timestep.
    """
    pass
        
    