from AeViz.units.units import units
import numpy as np

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
        ene_mag = np.sum(simulation.magnetic_energy(file_name)[0][mask] * \
             dV[mask])
    ene_grav = np.sum(simulation.gravitational_energy(file_name)[mask] * \
        dV[mask])
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

def PNS_mass_energy(simulation, file_name, PNS_radius, gcells, dV, grid, gr):
    """
    Calculates the mass, kinetic, magnetic, rotational, gravitational
    and convective energy of the PNS for one timestep.
    """
    if simulation.dim == 1:
        mask = (simulation.cell.radius(simulation.ghost) <= \
            PNS_radius)
    else:
        mask = (simulation.cell.radius(simulation.ghost) <= \
            simulation.ghost.remove_ghost_cells_radii(PNS_radius,
                                                      simulation.dim,
                                                  **gcells)[..., None])
    rho = simulation.rho(file_name)[mask] * dV[mask]
    ene_grav = np.sum(simulation.gravitational_energy(file_name)[mask] \
        * dV[mask])
    if simulation.dim == 1:
        ene_rot = 0
        ene_kin = np.sum(0.5 * rho * \
            simulation.radial_velocity(file_name)[mask] ** 2)
        ene_mag = 0
        conv_ene = 0
        Lx, Ly, Lz = 0, 0, 0
        L_tot = 0
    else:
        vr = simulation.radial_velocity(file_name)
        vt = simulation.theta_velocity(file_name)
        vp = simulation.phi_velocity(file_name)
        ene_rot = np.sum(0.5 * rho * \
            vp[mask] ** 2)
        ene_kin = np.sum(0.5 * rho * \
            (vr[mask] ** 2 + \
            vt[mask] ** 2))
        ene_mag = np.sum(simulation.magnetic_energy(file_name)[0][mask] \
            * dV[mask])
        conv_ene = np.sum(vt[mask] ** 2 \
            * rho)
        vx, vy, vz = gr.velocity_sph_to_cart(vr, vt, vp)
        masked_grid = [i[mask] if not type(i) == int else i for i in grid] 
        Lx, Ly, Lz = (rho * (vz[mask] * masked_grid[1] - vy[mask] * \
            masked_grid[2])).sum(), \
            (rho * (vx[mask] * masked_grid[2] - vz[mask] * \
                masked_grid[0])).sum(), \
            (rho * (vy[mask] * masked_grid[0] - vx[mask] * \
                masked_grid[1])).sum()
        L_tot = np.sqrt(Lx ** 2 + Ly ** 2 + Lz ** 2)
    return u.convert_to_solar_masses(np.sum(rho)), ene_kin, ene_mag, ene_rot, \
        ene_grav, ene_kin + ene_rot, conv_ene, Lx, Ly, Lz, L_tot


def unbound_mass_energy(simulation, file_name, dV):
    """
    Calculates the explosion energy and the unbound mass
    for one timestep. It also returns the ejecta kinetic and magnetic
    energy.
    """
    rho = simulation.rho(file_name)
    mhd_ene = simulation.MHD_energy(file_name) + \
        simulation.gravitational_energy(file_name)
    mask = (mhd_ene > 0) & (simulation.cell.radius(simulation.ghost) < 1e10)
    ej_mass = u.convert_to_solar_masses(np.sum(rho[mask] * dV[mask]))
    expl_ene = np.sum(mhd_ene[mask] * dV[mask])
    if simulation.dim == 1:
        ej_kin = np.sum(0.5 * rho[mask] * \
            simulation.radial_velocity(file_name)[mask] ** 2)
        ej_mag = 0
    else:
        ej_kin = np.sum(0.5 * dV[mask] * rho[mask] * \
            (simulation.radial_velocity(file_name)[mask] ** 2 + \
            simulation.theta_velocity(file_name)[mask] ** 2 + \
            simulation.phi_velocity(file_name)[mask] ** 2))
        ej_mag = np.sum(simulation.magnetic_energy(file_name)[0][mask] * \
            dV[mask])
    return ej_mass, expl_ene, ej_kin, ej_mag

def mass_flux(simulation, file_name, dOmega, radius_index):
    """
    Callulates the flow of matter at the index
    """
    if simulation.dim == 1:
        return -u.convert_to_solar_masses(4 * np.pi * \
                simulation.cell.radius(simulation.ghost)[radius_index] ** 2 * \
                np.sum(simulation.radial_velocity(file_name)[radius_index] * \
                simulation.rho(file_name)[radius_index]))
    
    return  -u.convert_to_solar_masses(
        simulation.cell.radius(simulation.ghost)[radius_index] ** 2 * \
            np.sum(dOmega *\
                simulation.radial_velocity(file_name)[..., radius_index] * \
            simulation.rho(file_name)[..., radius_index]))
    
