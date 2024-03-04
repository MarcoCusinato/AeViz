import numpy as np
import h5py, os
from AeViz.utils.utils import (check_existence, progressBar, checkpoints)
from AeViz.utils.file_utils import save_hdf
from AeViz.units.units import units

u = units()

def velocity_kick(simulation, file_name, PNS_radius, gcells, dV, dOmega, r400):
    """
    Calculates the velocity kick of the PNS for one timestep.
    """
    if simulation.time(file_name) <= 0:
        return 0, 0, 0, [0, 0, 0]
    if simulation.dim == 1:
        mask = (simulation.cell.radius(simulation.ghost) >= \
            PNS_radius) & (simulation.cell.radius(simulation.ghost) <= 1e10)
    else:
        mask = (simulation.cell.radius(simulation.ghost) >= \
            simulation.ghost.remove_ghost_cells_radii(PNS_radius,
                                                      simulation.dim,
                                                  **gcells)[..., None]) & \
                (simulation.cell.radius(simulation.ghost) <= 1e10)
    rho = simulation.rho(file_name)[mask] * dV[mask]
    vr = np.sum(simulation.radial_velocity(file_name)[mask] * rho)
    if simulation.dim == 1:
        vtheta = 0
        vphi = 0
    else:
        vtheta = np.sum(simulation.theta_velocity(file_name)[mask] * rho)
        vphi = np.sum(simulation.phi_velocity(file_name)[mask] * rho)
    nu_flux = np.sum(simulation.neutrino_momenta_grey(file_name)[..., r400, :, 0] / \
        u.speed_light * 4e7 ** 2 * dOmega[..., None], axis=-1)
    return -vr, -vtheta, -vphi, -nu_flux


def calculate_kick(simulation, save_checkpoints=True):
    if check_existence(simulation, 'kick_velocity.h5'):
        time, vr, vtheta, vphi, nu_flux  = \
            read_kick(simulation)
        if len(simulation.hdf_file_list) == len(time):
            return time, vr, vtheta, vphi, nu_flux
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
        vr_f, vtheta_f, vphi_f, nu_flux_f = velocity_kick(simulation, file, 
                                                    PNS_radius[..., findex],
                                                    pgcells, dV, dOmega, r400)
        try:
            time = np.append(time, simulation.time(file))
            vr = np.append(vr, vr_f)
            vtheta = np.append(vtheta, vtheta_f)
            vphi = np.append(vphi, vphi_f)
            nu_flux['nue'] = np.append(nu_flux['nue'], nu_flux_f[0])
            nu_flux['nua'] = np.append(nu_flux['nua'], nu_flux_f[1])
            nu_flux['nux'] = np.append(nu_flux['nux'], nu_flux_f[2])
        except:
            time = np.array([simulation.time(file)])
            vr = np.array([vr_f])
            vtheta = np.array([vtheta_f])
            vphi = np.array([vphi_f])
            nu_flux['nue'] = np.array([nu_flux_f[0]])
            nu_flux['nua'] = np.array([nu_flux_f[1]])
            nu_flux['nux'] = np.array([nu_flux_f[2]])

        if (check_index >= checkpoint) and save_checkpoints:
            print('Checkpoint reached, saving...\n')
            save_hdf(os.path.join(simulation.storage_path, 
                                  'kick_velocity.h5'),
                     ['time', 'vr', 'vtheta', 'vphi', 'nu_flux'],
                     [time, vr, vtheta, vphi, nu_flux])
            check_index = 0
        check_index += 1
        progress_index += 1
        findex += 1
    print('Computation completed, saving...')
    save_hdf(os.path.join(simulation.storage_path, 
                                  'kick_velocity.h5'),
                     ['time', 'vr', 'vtheta', 'vphi', 'nu_flux'],
                     [time, vr, vtheta, vphi, nu_flux])
    return time, vr, vtheta, vphi, nu_flux   
    

def read_kick(simulation):
    """
    Reads the kick file.
    """
    with h5py.File(os.path.join(simulation.storage_path, 
                                  'kick_velocity.h5'), 'r') as f:
        time = f['time'][...]
        vr = f['vr'][...]
        vtheta = f['vtheta'][...]
        vphi = f['vphi'][...]
        nu_flux = {'nue': f['nu_flux']['nue'][...],
                   'nua': f['nu_flux']['nua'][...],
                   'nux': f['nu_flux']['nux'][...]}
    return time, vr, vtheta, vphi, nu_flux
    