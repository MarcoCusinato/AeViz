from AeViz.utils.physics.radii_utils import (PNS_radius, innercore_radius, gain_radius,
                                     neutrino_sphere_radii, PNS_nucleus,
                                     shock_radius, isodensity_radii)
from AeViz.utils.files.file_utils import save_hdf, create_series
import numpy as np
from typing import Literal
import os, h5py
from AeViz.utils.math_utils import function_average_radii
from AeViz.utils.utils import progressBar, check_existence, checkpoints
from AeViz.units.aeseries import aerray, aeseries
from AeViz.units import u
from AeViz.utils.files.string_utils import merge_strings

functions = {
    'PNS': PNS_radius,
    'innercore': innercore_radius,
    'gain': gain_radius,
    'neutrino': neutrino_sphere_radii,
    'shock': shock_radius,
    'nucleus': PNS_nucleus,
    'isodensity': isodensity_radii
}

save_names = {
    'PNS': 'PNS_radius.h5',
    'innercore': 'innercore_radius.h5',
    'gain': 'gain_radius.h5',
    'neutrino': 'neutrino_sphere_radii.h5',
    'shock': 'shock_radius.h5',
    'nucleus': 'PNS_nucleus.h5',
    'isodensity': 'isodensities_radii.h5'
}

def calculate_radius(simulation, radius:Literal['PNS', 'innercore', 'gain', 
                                             'neutrino', 'shock', 'nucleus',
                                             'isodensity'],
                     save_checkpoints=True, rmax=None):
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
                save_hdf(os.path.join(simulation.storage_path,
                                      save_names[radius]),
                        ['time', 'radii', 'max', 'min', 'avg', 'gcells',
                         'processed'],
                        [time, full_radius, max_radius, min_radius, avg_radius,
                        simulation.ghost.return_ghost_dictionary(),
                        simulation.hdf_file_list])
                return create_series(time, full_radius, max_radius,
                                            min_radius, avg_radius, ghost_cells)
            else:
                start_point = 0
                time = 0
                full_radius = 0
                max_radius = 0
                min_radius = 0
                avg_radius = 0
                ghost_cells = 0
                processed_hdf = []
        elif processed_hdf[-1].decode("utf-8") == simulation.hdf_file_list[-1]:
            return create_series(time, full_radius, max_radius,
                                        min_radius, avg_radius, ghost_cells)
        else:
            start_point = len(processed_hdf)
            processed_hdf = [ff.decode("utf-8") for ff in processed_hdf]
            print('Checkpoint found for ' + radius + ' radius, starting' \
                ' from checkpoint.\nPlease wait...')
    else:
        start_point = 0
        processed_hdf = []
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
        PNS_rad, _ = simulation.PNS_radius(rad='full')
    for file in simulation.hdf_file_list[start_point:]:
        progressBar(progress_index, total_points, suffix='Computing...')
        if radius == 'gain':
            rad_step = functions[radius](simulation, file,
                                         PNS_rad.data[..., check_index])   
        else:
            try:
                if rmax is None:
                    rad_step = functions[radius](simulation, file)
                else:
                    rad_step = functions[radius](simulation, file, rmax)
            except KeyError:
                print('Missing dataset in file ' + file + \
                    ', skipping but adding as processed...')
                check_index += 1
                progress_index += 1
                processed_hdf.append(file)
                continue
            except Exception as e:
                print('Error in file ' + file)
                raise e
        if radius == 'neutrino':
            try:
                time = np.concatenate((time, simulation.time(file)))
                full_radius = {key: np.concatenate((full_radius[key], 
                                                    rad_step[i][..., None]),
                                                    axis=-1)
                               for (key, i) in zip(['nue', 'nua', 'nux'],
                                                   range(3))}
                nog_rad_step = {key: simulation.ghost.remove_ghost_cells_radii(
                    rad_step[i], simulation.dim) for (key, i) in
                    zip(['nue', 'nua', 'nux'], range(3))}
                max_radius = {key: np.concatenate((max_radius[key],
                                                   np.nanmax(nog_rad_step[key])))
                              for key in ['nue', 'nua', 'nux']}
                min_radius = {key: np.concatenate((min_radius[key], np.nanmin(nog_rad_step[key])))
                              for key in ['nue', 'nua', 'nux']}
                avg_radius = {key: np.concatenate((avg_radius[key],
                                function_average_radii(nog_rad_step[key],
                                                       simulation.dim, dOmega)))
                                for key in ['nue', 'nua', 'nux']}
            except Exception as e:
                print(e)
                time = simulation.time(file)
                full_radius = {key: rad_step[i][..., None] for (key, i) in
                               zip(['nue', 'nua', 'nux'], range(3))}
                nog_rad_step = {key: simulation.ghost.remove_ghost_cells_radii(
                    rad_step[i], simulation.dim) for (key, i) in
                    zip(['nue', 'nua', 'nux'], range(3))}
                max_radius = {key: np.nanmax(nog_rad_step[key])
                              for key in ['nue', 'nua', 'nux']}
                min_radius = {key: np.nanmin(nog_rad_step[key])
                              for key in ['nue', 'nua', 'nux']}
                avg_radius = {key: function_average_radii(
                                    nog_rad_step[key], simulation.dim,
                                    dOmega) for key in
                                ['nue', 'nua', 'nux']}
        elif radius == 'isodensity':
            try:
                time = np.concatenate((time, simulation.time(file)))
                full_radius = {key: np.concatenate((full_radius[key], 
                                                    rad_step[key][..., None]),
                                                    axis=-1)
                               for key in rad_step.keys()}
                nog_rad_step = {key: simulation.ghost.remove_ghost_cells_radii(
                    rad_step[key], simulation.dim) for key in rad_step.keys()}
                max_radius = {key: np.concatenate((max_radius[key],
                                                   np.nanmax(nog_rad_step[key])))
                              for key in nog_rad_step.keys()}
                min_radius = {key: np.concatenate((min_radius[key], np.nanmin(nog_rad_step[key])))
                              for key in nog_rad_step.keys()}
                avg_radius = {key: np.concatenate((avg_radius[key],
                                function_average_radii(nog_rad_step[key],
                                                       simulation.dim, dOmega)))
                                for key in nog_rad_step.keys()}
            except Exception as e:
                print(e)
                time = simulation.time(file)
                full_radius = {key: rad_step[key][..., None] for key in 
                               rad_step.keys()}
                nog_rad_step = {key: simulation.ghost.remove_ghost_cells_radii(
                    rad_step[key], simulation.dim) for key in rad_step.keys()}
                max_radius = {key: np.nanmax(nog_rad_step[key])
                              for key in nog_rad_step.keys()}
                min_radius = {key: np.nanmin(nog_rad_step[key])
                              for key in nog_rad_step.keys()}
                avg_radius = {key: function_average_radii(nog_rad_step[key],
                                                          simulation.dim,
                                                          dOmega)
                              for key in nog_rad_step.keys()}
        else:
            nog_rad_step = simulation.ghost.remove_ghost_cells_radii(rad_step,
                                                               simulation.dim)
            max_rad_step = np.nanmax(nog_rad_step)
            min_rad_step = np.nanmin(nog_rad_step)
            avg_rad_step = function_average_radii(nog_rad_step, 
                                                simulation.dim, dOmega)
            try:
                time = np.concatenate((time, simulation.time(file)))
                full_radius = np.concatenate((full_radius,
                                            rad_step[..., None]),
                                            axis=-1)
                max_radius = np.concatenate((max_radius, max_rad_step))
                min_radius = np.concatenate((min_radius, min_rad_step))
                avg_radius = np.concatenate((avg_radius, avg_rad_step))
            except Exception as e:
                print(e)
                time = simulation.time(file)
                full_radius = rad_step[..., None]
                max_radius = max_rad_step
                min_radius = min_rad_step
                avg_radius = avg_rad_step
        processed_hdf.append(file)
        if (check_index >= checkpoint and save_checkpoints):
            print('Checkpoint reached, saving...\n')
            save_hdf(os.path.join(simulation.storage_path, save_names[radius]),
                     ['time', 'radii', 'max', 'min', 'avg', 'gcells',
                      'processed'],
                     [time, full_radius, max_radius, min_radius, avg_radius,
                      simulation.ghost.return_ghost_dictionary(),
                      processed_hdf])
            check_index = 0
        check_index += 1
        progress_index += 1
    print('Computation completed, saving...')
    save_hdf(os.path.join(simulation.storage_path, save_names[radius]),
                ['time', 'radii', 'max', 'min', 'avg', 'gcells', 'processed'],
                [time, full_radius, max_radius, min_radius, avg_radius,
                simulation.ghost.return_ghost_dictionary(), processed_hdf])
    print('Done!')
    ## Unpack the ghost cells
    out_gcells = {}
    for key in simulation.ghost.return_ghost_dictionary().keys():
        out_gcells[key[0]+'_l'] = simulation.ghost.return_ghost_dictionary()[key][0]
        out_gcells[key[0]+'_r'] = simulation.ghost.return_ghost_dictionary()[key][1]
    simulation.ghost.restore_default()
    time, full_radius, max_radius, min_radius, avg_radius, ghost_cells, \
            _ = read_radius(simulation, radius)
    return create_series(time, full_radius, max_radius, min_radius,
                                avg_radius, out_gcells)
   
def read_radius(simulation,
                radius:Literal['PNS', 'innercore', 'gain',
                               'neutrino', 'shock', 'nucleus', 'isodensity'],
                ):
    """
    Reads the data from the hdf file. Returns a tuple containing:
    (time, radii, max, min, avg, ghost_cells)
    In case of neutrinos, radii, max, min and avg are dictionaries.
    """
    radius_data = h5py.File(os.path.join(simulation.storage_path, 
                                            save_names[radius]), 'r')
    if radius == 'neutrino':
        data = [
                aerray(radius_data['time'][...], u.s, 'time',
                       r'$t-t_\mathrm{b}$', None,
                       [-0.005, radius_data['time'][-1]]),
                {'nue': aerray(radius_data['radii/nue'][...], u.cm, 'Rnue',
                               r'$R_{\nu_e}$', None, [0, 1.5e7]),
                 'nua': aerray(radius_data['radii/nua'][...], u.cm, 'Rnux',
                               r'$R_{\overline{\nu}_e}$', None, [0, 1.5e7]),
                 'nux': aerray(radius_data['radii/nux'][...], u.cm, 'Rnux',
                               r'$R_{\overline{\nu}_x}$', None, [0, 1.5e7])},
                {'nue': aerray(radius_data['max/nue'][...], u.cm, 'Rnue_max',
                               r'$R_{\nu_e,max}$', None, [0, 1.5e7]),
                 'nua': aerray(radius_data['max/nua'][...], u.cm, 'Rnua_max',
                               r'$R_{\overline{\nu}_e,max}$', None, [0, 1.5e7]),
                 'nux': aerray(radius_data['max/nux'][...], u.cm, 'Rnuxmax',
                               r'$R_{\nu_x,max}$', None, [0, 1.5e7])},
                {'nue': aerray(radius_data['min/nue'][...], u.cm, 'Rnue_min',
                               r'$R_{\nu_e,min}$', None, [0, 1.5e7]),
                 'nua': aerray(radius_data['min/nua'][...], u.cm, 'Rnua_min',
                               r'$R_{\overline{\nu}_e,min}$', None, [0, 1.5e7]),
                 'nux': aerray(radius_data['min/nux'][...], u.cm, 'Rnux_min',
                               r'$R_{\nu_x,min}$', None, [0, 1.5e7])},
                {'nue': aerray(radius_data['avg/nue'][...], u.cm, 'Rnue_avg',
                               r'$R_{\nu_e,avg}$', None, [0, 1.5e7]),
                 'nua': aerray(radius_data['avg/nua'][...], u.cm, 'Rnua_avg',
                               r'$R_{\overline{\nu}_e,avg}$', None, [0, 1.5e7]),
                 'nux': aerray(radius_data['avg/nux'][...], u.cm, 'Rnux_avg',
                               r'$R_{\nu_x,avg}$', None, [0, 1.5e7])},
                {'p_l': list(radius_data['gcells/phi'])[0],
                 'p_r':  list(radius_data['gcells/phi'])[1],
                 't_l':  list(radius_data['gcells/theta'])[0],
                 't_r':  list(radius_data['gcells/theta'])[1],
                 'r_l':  list(radius_data['gcells/radius'])[0],
                 'r_r':  list(radius_data['gcells/radius'])[1]}]
    elif radius == 'isodensity':
        keys = ['1e+14', '1e+13', '1e+12', '1e+10', '1e+09', '1e+08']
        data = [
                aerray(radius_data['time'][...], u.s, 'time',
                       r'$t-t_\mathrm{b}$', None,
                       [-0.005, radius_data['time'][-1]]),
                {key: aerray(radius_data[f'radii/{key}'], u.cm, f'R{key}',
                             merge_strings(r'$R_\mathrm{',
                                           key.replace('+', '')),
                             None, [0, 1.5e7], False) for key in keys},
                {key: aerray(radius_data[f'max/{key}'], u.cm, f'R{key}',
                             merge_strings(r'$R_\mathrm{',
                                           key.replace('+', ''), r',max}$'),
                             None, [0, 1.5e7], False) for key in keys},
                {key: aerray(radius_data[f'min/{key}'], u.cm, f'R{key}',
                             merge_strings(r'$R_\mathrm{',
                                           key.replace('+', ''), r',min}$'),
                             None, [0, 1.5e7], False) for key in keys},
                {key: aerray(radius_data[f'avg/{key}'], u.cm, f'R{key}',
                             merge_strings(r'$R_\mathrm{',
                                           key.replace('+', ''), r',avg}$'),
                             None, [0, 1.5e7], False) for key in keys},
                {'p_l': list(radius_data['gcells/phi'])[0],
                 'p_r':  list(radius_data['gcells/phi'])[1],
                 't_l':  list(radius_data['gcells/theta'])[0],
                 't_r':  list(radius_data['gcells/theta'])[1],
                 'r_l':  list(radius_data['gcells/radius'])[0],
                 'r_r':  list(radius_data['gcells/radius'])[1]}
        ]
    else:
        lab = radius if len(radius) <=5 else radius[:3]
        lm = [0, 1.5e7] if radius in ['PNS', 'innercore'] else [1e6, 1e9] if radius in ['gain', 'shock'] else [0, 4e6]
        lg = True if radius in ['gain', 'shock'] else False
        data = [
                aerray(radius_data['time'][...], u.s, 'time',
                       r'$t-t_\mathrm{b}$', None,
                       [-0.005, radius_data['time'][-1]]),
                aerray(radius_data['radii'][...], u.cm, 'R_'+lab,
                      merge_strings(r'$R_\mathrm{', lab, r'}$'), None, lm, lg),
                aerray(radius_data['max'][...], u.cm, 'R_'+lab+'_max',
                      merge_strings(r'$R_\mathrm{', lab, r',max}$'), None, lm, lg),
                aerray(radius_data['min'][...], u.cm, 'R_'+lab+'_min',
                      merge_strings(r'$R_\mathrm{', lab, r',min}$'), None, lm, lg),
                aerray(radius_data['avg'][...], u.cm, 'R_'+lab+'_avg',
                      merge_strings(r'$R_\mathrm{', lab, r',avg}$'), None, lm, lg),
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

