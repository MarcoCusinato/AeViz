import numpy as np
import h5py, os
from AeViz.utils.utils import (check_existence, progressBar, checkpoints)
from AeViz.utils.file_utils import save_hdf
from AeViz.units.units import units
from AeViz.grid.grid import grid

u = units()

def velocity_kick(simulation, file_name, PNS_radius, gcells, dV, dOmega, r400):
    """
    Calculates the velocity kick of the PNS for one timestep.
    """
    if simulation.time(file_name) <= 0 or simulation.dim == 1:
        return 0., 0., 0., [0., 0., 0.]
    
    mask = (simulation.cell.radius(simulation.ghost) >= \
        simulation.ghost.remove_ghost_cells_radii(PNS_radius,
                                                    simulation.dim,
                                                **gcells)[..., None]) & \
            (simulation.cell.radius(simulation.ghost) <= 1e9)
    rho = simulation.rho(file_name)[mask] * dV[mask]
    gh = grid(simulation.dim, simulation.cell.radius(simulation.ghost),
            simulation.cell.theta(simulation.ghost),
            simulation.cell.phi(simulation.ghost))
    vx, vy, vz = gh.velocity_sph_to_cart(simulation.radial_velocity(file_name),
                                          simulation.theta_velocity(file_name),
                                          simulation.phi_velocity(file_name))
    nu_flux = u.convert_to_solar_masses(
        np.sum(simulation.neutrino_momenta_grey(file_name)[..., r400, :, 0] / \
        u.speed_light * 4e7 ** 2 * dOmega[..., None], axis=-1))
    vz = u.convert_to_solar_masses(np.sum(vz[mask] * rho))
    if simulation.dim == 2:
        vy = 0.0
        vx = 0.0
    else:
        vy = np.sum(vy[mask] * rho)
        vx = np.sum(vx[mask] * rho)
    
    return -vx, -vy, -vz, -nu_flux


def calculate_kick(simulation, save_checkpoints=True):
    if check_existence(simulation, 'kick_velocity.h5'):
        time, vx, vy, vz, nu_flux  = \
            read_kick(simulation)
        if len(simulation.hdf_file_list) == len(time):
            return time, vx, vy, vz, nu_flux
        else:
            start_point = len(time)
            print('Checkpoint found for the kick file, starting' \
                  ' from checkpoint.\nPlease wait...')
    else:
        start_point = 0
        print('No checkpoint found for the kick file, starting' \
              ' from the beginning.\nPlease wait...')
    if (checkpoints[simulation.dim] == False) or (not save_checkpoints):
        checkpoint = len(simulation.hdf_file_list)
    else:
        checkpoint = checkpoints[simulation.dim]
    ## Get the radii
    _, PNS_radius, _, _, _, pgcells = simulation.PNS_radius()
    ## Get the grid
    dV = simulation.cell.dVolume_integration(simulation.ghost)
    dOmega = simulation.cell.dOmega(simulation.ghost)
    r400 = np.argmax(simulation.cell.radius(simulation.ghost) >= \
                     u.convert_to_cm(400))
    ## Indices
    findex = start_point
    check_index = 0
    progress_index = 0
    total_points = len(simulation.hdf_file_list) - start_point
    nu_flux = {}
    for file in simulation.hdf_file_list[start_point:]:
        progressBar(progress_index, total_points, suffix='Computing...')
        vx_f, vy_f, vz_f, nu_flux_f = velocity_kick(simulation, file, 
                                                    PNS_radius[..., findex],
                                                    pgcells, dV, dOmega, r400)
        try:
            time = np.append(time, simulation.time(file))
            vx = np.append(vx, vx_f)
            vy = np.append(vy, vy_f)
            vz = np.append(vz, vz_f)
            nu_flux['nue'] = np.append(nu_flux['nue'], nu_flux_f[0])
            nu_flux['nua'] = np.append(nu_flux['nua'], nu_flux_f[1])
            nu_flux['nux'] = np.append(nu_flux['nux'], nu_flux_f[2])
        except:
            time = np.array([simulation.time(file)])
            vx = np.array([vx_f])
            vy = np.array([vy_f])
            vz = np.array([vz_f])
            nu_flux['nue'] = np.array([nu_flux_f[0]])
            nu_flux['nua'] = np.array([nu_flux_f[1]])
            nu_flux['nux'] = np.array([nu_flux_f[2]])

        if (check_index >= checkpoint) and save_checkpoints:
            print('Checkpoint reached, saving...\n')
            save_hdf(os.path.join(simulation.storage_path, 
                                  'kick_velocity.h5'),
                     ['time', 'vx', 'vy', 'vz', 'nu_flux'],
                     [time, vx, vy, vz, nu_flux])
            check_index = 0
        check_index += 1
        progress_index += 1
        findex += 1
    print('Computation completed, saving...')
    save_hdf(os.path.join(simulation.storage_path, 
                                  'kick_velocity.h5'),
                     ['time', 'vx', 'vy', 'vz', 'nu_flux'],
                     [time, vx, vy, vz, nu_flux])
    return time, vx, vy, vz, nu_flux   
    

def read_kick(simulation):
    """
    Reads the kick file.
    """
    with h5py.File(os.path.join(simulation.storage_path, 
                                  'kick_velocity.h5'), 'r') as f:
        time = f['time'][...]
        vx = f['vx'][...]
        vy = f['vy'][...]
        vz = f['vz'][...]
        nu_flux = {'nue': f['nu_flux']['nue'][...],
                   'nua': f['nu_flux']['nua'][...],
                   'nux': f['nu_flux']['nux'][...]}
    return time, vx, vy, vz, nu_flux
    