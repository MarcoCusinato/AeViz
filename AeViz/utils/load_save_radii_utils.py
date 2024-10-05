from AeViz.utils.radii_utils import (PNS_radius, innercore_radius, gain_radius,
                                     neutrino_sphere_radii, PNS_nucleus,
                                     shock_radius)
from AeViz.utils.file_utils import save_hdf
import numpy as np
from typing import Literal
import os, h5py
from AeViz.utils.math_utils import function_average_radii
from AeViz.utils.utils import progressBar, check_existence, checkpoints


functions = {
    'PNS': PNS_radius,
    'innercore': innercore_radius,
    'gain': gain_radius,
    'neutrino': neutrino_sphere_radii,
    'shock': shock_radius,
    'nucleus': PNS_nucleus
}

save_names = {
    'PNS': 'PNS_radius.h5',
    'innercore': 'innercore_radius.h5',
    'gain': 'gain_radius.h5',
    'neutrino': 'neutrino_sphere_radii.h5',
    'shock': 'shock_radius.h5',
    'nucleus': 'PNS_nucleus.h5'
}


def calculate_radius(simulation, radius:Literal['PNS', 'innercore', 'gain', 
                                             'neutrino', 'shock', 'nucleus'],
                     save_checkpoints=True):
    """
    Calculates the selected radius for each timestep of the simulation.
    In case of neutrinos, since in some cases are not saved for each
    timestep, the timestep for which they are not saved is skipped.
    If the save checkpoint flag is on, every 200 timesteps (for 3D
    simulations) or 600 timesteps (for 2D simulations) the data is saved
    in a hdf file.
    Input:

    """
    if check_existence(simulation, save_names[radius]):
        time, full_radius, max_radius, min_radius, avg_radius, ghost_cells, \
            processed_hdf = \
            read_radius(simulation, radius)
        ## Retrocompatibility option
        if processed_hdf is None:
            if len(simulation.hdf_file_list) == len(time):
                save_hdf(os.path.join(simulation.storage_path, save_names[radius]),
                    ['time', 'radii', 'max', 'min', 'avg', 'gcells', 'processed'],
                    [time, full_radius, max_radius, min_radius, avg_radius,
                    simulation.ghost.return_ghost_dictionary(),
                    simulation.hdf_file_list])
                return time, full_radius, max_radius, min_radius, avg_radius, \
                    ghost_cells
            else:
                start_point = 0
        elif processed_hdf[-1] == simulation.hdf_file_list[-1]:
            return time, full_radius, max_radius, min_radius, avg_radius, \
                ghost_cells
        else:
            start_point = len(processed_hdf)
            print('Checkpoint found for ' + radius + ' radius, starting' \
                ' from checkpoint.\nPlease wait...')
    else:
        start_point = 0
        print('No checkpoint found for ' + radius + ' radius, starting from' \
            ' the beginning.\nPlease wait...')
    if (checkpoints[simulation.dim] == False) or (not save_checkpoints):
        checkpoint = len(simulation.hdf_file_list)
    else:
        checkpoint = checkpoints[simulation.dim]
    g = 0
    dOmega = simulation.cell.dOmega(simulation.ghost)
    simulation.ghost.update_ghost_cells(t_l = g, t_r = g, p_l = g, p_r = g)
    check_index = 0
    progress_index = 0
    total_points = len(simulation.hdf_file_list) - start_point
    if radius == 'gain':
        _, PNS_rad, _, _, _, _ = calculate_radius(simulation, 'PNS')
    processed_hdf = np.array([])
    for file in simulation.hdf_file_list[start_point:]:
        progressBar(progress_index, total_points, suffix='Computing...')
        if radius == 'gain':
            rad_step = functions[radius](simulation, file,
                                         PNS_rad[..., check_index])   
        else:
            try:
                rad_step = functions[radius](simulation, file)
            except Exception as ex:
                print('Error in file ' + file + ', skipping...')
                print(ex)
                print(ex == 'Unable to open object (component not found)')
                check_index += 1
                progress_index += 1
                continue
        if radius == 'neutrino':
            try:
                time = np.concatenate((time, simulation.time(file)))
                full_radius = {key: np.concatenate((full_radius[key], 
                                                    rad_step[..., i, None]),
                                                    axis=-1)
                               for (key, i) in zip(['nue', 'nua', 'nux'],
                                                   range(3))}
                nog_rad_step = {key: simulation.ghost.remove_ghost_cells_radii(
                    rad_step[..., i], simulation.dim)
                                for (key, i) in zip(['nue', 'nua', 'nux'],
                                                    range(3))}
                max_radius = {key: np.concatenate((max_radius[key],
                                                    np.array([np.nanmax(
                                                        nog_rad_step[key])])))
                              for key in ['nue', 'nua', 'nux']}
                min_radius = {key: np.concatenate((min_radius[key],
                                                    np.array([np.nanmin(
                                                        nog_rad_step[key])])))
                              for key in ['nue', 'nua', 'nux']}
                avg_radius = {key: np.concatenate((avg_radius[key],
                                np.array([function_average_radii(
                                    nog_rad_step[key], simulation.dim,
                                    dOmega)]))) for key in
                                ['nue', 'nua', 'nux']}
                
            except:
                time = simulation.time(file)
                full_radius = {key: rad_step[..., i, None] for (key, i) in
                               zip(['nue', 'nua', 'nux'], range(3))}
                nog_rad_step = {key: simulation.ghost.remove_ghost_cells_radii(
                    rad_step[..., i], simulation.dim)
                                for (key, i) in zip(['nue', 'nua', 'nux'],
                                                    range(3))}
                max_radius = {key: np.array([np.nanmax(nog_rad_step[key])])
                              for key in ['nue', 'nua', 'nux']}
                min_radius = {key: np.array([np.nanmin(nog_rad_step[key])])
                              for key in ['nue', 'nua', 'nux']}
                avg_radius = {key: np.array([function_average_radii(
                                    nog_rad_step[key], simulation.dim,
                                    dOmega)]) for key in
                                ['nue', 'nua', 'nux']}
        else:
            nog_rad_step = simulation.ghost.remove_ghost_cells_radii(rad_step,
                                                               simulation.dim)
            max_rad_step = np.array([np.nanmax(nog_rad_step)])
            min_rad_step = np.array([np.nanmin(nog_rad_step)])
            avg_rad_step = np.array([function_average_radii(nog_rad_step, 
                                                simulation.dim, dOmega)])
            try:
                time = np.concatenate((time, simulation.time(file)))
                full_radius = np.concatenate((full_radius,
                                            rad_step[..., None]),
                                            axis=-1)
                max_radius = np.concatenate((max_radius, max_rad_step))
                min_radius = np.concatenate((min_radius, min_rad_step))
                avg_radius = np.concatenate((avg_radius, avg_rad_step))
            except:
                time = simulation.time(file)
                full_radius = rad_step[..., None]
                max_radius = max_rad_step
                min_radius = min_rad_step
                avg_radius = avg_rad_step
        processed_hdf = np.append(processed_hdf, file)
        if (check_index >= checkpoint and save_checkpoints):
            print('Checkpoint reached, saving...\n')
            save_hdf(os.path.join(simulation.storage_path, save_names[radius]),
                     ['time', 'radii', 'max', 'min', 'avg', 'gcells', 'processed'],
                     [time, full_radius, max_radius, min_radius, avg_radius,
                      simulation.ghost.return_ghost_dictionary(), processed_hdf])
            check_index = 0
        check_index += 1
        progress_index += 1
    print('Computation completed, saving...')
    save_hdf(os.path.join(simulation.storage_path, save_names[radius]),
                ['time', 'radii', 'max', 'min', 'avg', 'gcells', 'processed'],
                [time, full_radius, max_radius, min_radius, avg_radius,
                simulation.ghost.return_ghost_dictionary(), processed_hdf])
    print('Done!')
    simulation.ghost.restore_default()
    return time, full_radius, max_radius, min_radius, avg_radius, \
        simulation.ghost.return_ghost_dictionary()

    
def read_radius(simulation, radius:Literal['PNS', 'innercore', 'gain', 
                                             'neutrino', 'shock', 'nucleus']):
    """
    Reads the data from the hdf file. Returns a tuple containing:
    (time, radii, max, min, avg, ghost_cells)
    In case of neutrinos, radii, max, min and avg are dictionaries.
    """
    radius_data = h5py.File(os.path.join(simulation.storage_path, 
                                            save_names[radius]), 'r')
    if radius == 'neutrino':
        data = [
                radius_data['time'][...],
                {'nue': radius_data['radii/nue'][...],
                'nua': radius_data['radii/nua'][...],
                'nux': radius_data['radii/nux'][...]},
                {'nue': radius_data['max/nue'][...],
                    'nua': radius_data['max/nua'][...],
                    'nux': radius_data['max/nux'][...]},
                {'nue': radius_data['min/nue'][...],
                'nua': radius_data['min/nua'][...],
                'nux': radius_data['min/nux'][...]},
                {'nue': radius_data['avg/nue'][...],
                    'nua': radius_data['avg/nua'][...],
                    'nux': radius_data['avg/nux'][...]},
                {'p_l': list(radius_data['gcells/phi'])[0],
                    'p_r':  list(radius_data['gcells/phi'])[1],
                    't_l':  list(radius_data['gcells/theta'])[0],
                    't_r':  list(radius_data['gcells/theta'])[1],
                    'r_l':  list(radius_data['gcells/radius'])[0],
                    'r_r':  list(radius_data['gcells/radius'])[1]}]
    else:
        data = [
                radius_data['time'][...],
                radius_data['radii'][...],
                radius_data['max'][...],
                radius_data['min'][...],
                radius_data['avg'][...],
                {'p_l': list(radius_data['gcells/phi'])[0],
                    'p_r':  list(radius_data['gcells/phi'])[1],
                    't_l':  list(radius_data['gcells/theta'])[0],
                    't_r':  list(radius_data['gcells/theta'])[1],
                    'r_l':  list(radius_data['gcells/radius'])[0],
                    'r_r':  list(radius_data['gcells/radius'])[1]}]
    if 'processed' in radius_data:
        data.append(radius_data['processed'][...])
    else:
        data.append(None)
    radius_data.close()
    return data
        

