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
        return [0., 0., 0.], [0., 0., 0.], [0., 0., 0.], [0., 0., 0.]
    
    outer_radius = 1e10
    if outer_radius >= simulation.cell.radius(simulation.ghost)[-1]:
        outer_radius = simulation.cell.radius(simulation.ghost)[-2]
    
    mask = (simulation.cell.radius(simulation.ghost) >= \
        simulation.ghost.remove_ghost_cells_radii(PNS_radius,
                                                    simulation.dim,
                                                    **gcells)[..., None]) & \
            (simulation.cell.radius(simulation.ghost) <= outer_radius)
    rho = simulation.rho(file_name)[mask] * dV[mask]
    gh = grid(simulation.dim, simulation.cell.radius(simulation.ghost),
            simulation.cell.theta(simulation.ghost),
            simulation.cell.phi(simulation.ghost))
    vx, vy, vz = gh.velocity_sph_to_cart(simulation.radial_velocity(file_name),
                                          simulation.theta_velocity(file_name),
                                          simulation.phi_velocity(file_name))
    nu_flux = simulation.neutrino_momenta_grey(file_name)[..., r400, :, :] * \
                                          dOmega[..., None, None]
    nue_flux = list(gh.velocity_sph_to_cart(
                                        nu_flux[..., 0, 0][..., None],
                                        nu_flux[..., 0, 1][..., None],
                                        nu_flux[..., 0, 2][..., None]
                                       )
                    )
    nua_flux = list(gh.velocity_sph_to_cart(
                                        nu_flux[..., 1, 0][..., None],
                                        nu_flux[..., 1, 1][..., None],
                                        nu_flux[..., 1, 2][..., None]
                                        )
                    )
    nux_flux = list(gh.velocity_sph_to_cart(
                                        nu_flux[..., 2, 0][..., None],
                                        nu_flux[..., 2, 1][..., None],
                                        nu_flux[..., 2, 2][..., None]
                                        )
                    )
    nue_flux = [-u.convert_to_solar_masses(np.sum(comp)) \
                * (4e7 ** 2 / u.speed_light) for comp in nue_flux]
    nua_flux = [-u.convert_to_solar_masses(np.sum(comp)) \
                * (4e7 ** 2 / u.speed_light) for comp in nua_flux]
    nux_flux = [-u.convert_to_solar_masses(np.sum(comp)) \
                * (4 * 4e7 ** 2 / u.speed_light) for comp in nux_flux]
    
    vz = u.convert_to_solar_masses(np.sum(vz[mask] * rho))
    if simulation.dim == 2:
        vy = 0.0
        vx = 0.0
        for nu in [nue_flux, nua_flux, nux_flux]:
            nu[0] = 0.0
            nu[1] = 0.0
    else:
        vy = u.convert_to_solar_masses(np.sum(vy[mask] * rho))
        vx = u.convert_to_solar_masses(np.sum(vx[mask] * rho))
    
    return [-vx, -vy, -vz], nue_flux, nua_flux, nux_flux


def calculate_kick(simulation, save_checkpoints=True):
    if check_existence(simulation, 'kick_velocity.h5'):
        time, hydro_v, nu_flux, processed_hdf = \
            read_kick(simulation)
        ## Retrocompatibility with old files
        if processed_hdf is None:
            if len(simulation.hdf_file_list) == len(time):
                save_hdf(os.path.join(simulation.storage_path, 
                                  'kick_velocity.h5'),
                     ['time', 'hydro', 'nu_flux', 'processed'],
                     [time, hydro_v, nu_flux,
                      simulation.hdf_file_list])
                return time, hydro_v, nu_flux
            else:
                start_point = 0
                hydro_v = 0
                nu_flux = 0
                processed_hdf = []
        elif processed_hdf[-1].decode("utf-8") == simulation.hdf_file_list[-1]:
            return time, hydro_v, nu_flux
        else:
            start_point = len(processed_hdf)
            processed_hdf = [ff.decode("utf-8") for ff in processed_hdf]
            print('Checkpoint found for the kick file, starting' \
                  ' from checkpoint.\nPlease wait...')
    else:
        start_point = 0
        processed_hdf = []
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
    for file in simulation.hdf_file_list[start_point:]:
        progressBar(progress_index, total_points, suffix='Computing...')
        hydro_v_f, nue_flux_f, nua_flux_f, nux_flux_f = \
                                            velocity_kick(simulation, file, 
                                                    PNS_radius[..., findex],
                                                    pgcells, dV, dOmega, r400)
        try:
            time = np.append(time, simulation.time(file))
            hydro_v = {'x': np.append(hydro_v['x'], hydro_v_f[0]),
                       'y': np.append(hydro_v['y'], hydro_v_f[1]),
                       'z': np.append(hydro_v['z'], hydro_v_f[2])}
            nu_flux['nue'] = {'x': np.append(nu_flux['nue']['x'], nue_flux_f[0]),
                              'y': np.append(nu_flux['nue']['y'], nue_flux_f[1]),
                              'z': np.append(nu_flux['nue']['z'], nue_flux_f[2])}
            nu_flux['nua'] = {'x': np.append(nu_flux['nua']['x'], nua_flux_f[0]),
                              'y': np.append(nu_flux['nua']['y'], nua_flux_f[1]),
                              'z': np.append(nu_flux['nua']['z'], nua_flux_f[2])}
            nu_flux['nux'] = {'x': np.append(nu_flux['nux']['x'], nux_flux_f[0]),
                              'y': np.append(nu_flux['nux']['y'], nux_flux_f[1]),
                              'z': np.append(nu_flux['nux']['z'], nux_flux_f[2])}
        except:
            time = np.array([simulation.time(file)])
            hydro_v = {'x': np.array([hydro_v_f[0]]),
                       'y': np.array([hydro_v_f[1]]),
                       'z': np.array([hydro_v_f[2]])}
            nu_flux = {'nue': {'x': np.array([nue_flux_f[0]]),
                              'y': np.array([nue_flux_f[1]]),
                              'z': np.array([nue_flux_f[2]])},
                       'nua': {'x': np.array([nua_flux_f[0]]),
                              'y': np.array([nua_flux_f[1]]),
                              'z': np.array([nua_flux_f[2]])},
                       'nux': {'x': np.array([nux_flux_f[0]]),
                              'y': np.array([nux_flux_f[1]]),
                              'z': np.array([nux_flux_f[2]])}
                        }
        processed_hdf.append(file)
        if (check_index >= checkpoint) and save_checkpoints:
            print('Checkpoint reached, saving...\n')
            save_hdf(os.path.join(simulation.storage_path, 
                                  'kick_velocity.h5'),
                     ['time', 'hydro', 'nu_flux', 'processed'],
                     [time, hydro_v, nu_flux, processed_hdf])
            check_index = 0
        check_index += 1
        progress_index += 1
        findex += 1
    print('Computation completed, saving...')
    save_hdf(os.path.join(simulation.storage_path, 
                                  'kick_velocity.h5'),
                     ['time', 'hydro', 'nu_flux', 'processed'],
                     [time, hydro_v, nu_flux, processed_hdf])
    return time, hydro_v, nu_flux
    

def read_kick(simulation):
    """
    Reads the kick file.
    """
    with h5py.File(os.path.join(simulation.storage_path, 
                                  'kick_velocity.h5'), 'r') as f:
        time = f['time'][...]
        hydro_v = {'x': f['hydro/x'][...],
                   'y': f['hydro/y'][...],
                   'z': f['hydro/z'][...]}
        nu_flux = {'nue': {'x': f['nu_flux/nue/x'][...],
                           'y': f['nu_flux/nue/y'][...],
                           'z': f['nu_flux/nue/z'][...]},
                   'nua': {'x': f['nu_flux/nua/x'][...],
                           'y': f['nu_flux/nua/y'][...],
                           'z': f['nu_flux/nua/z'][...]},
                   'nux': {'x': f['nu_flux/nux/x'][...],
                           'y': f['nu_flux/nux/y'][...],
                           'z': f['nu_flux/nux/z'][...]}}
        if 'processed' in f:
            processed = f['processed'][...]
        else:
            processed = None
    return time, hydro_v, nu_flux, processed
    