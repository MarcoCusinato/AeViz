import numpy as np
from AeViz.utils.math_utils import function_average
import os, h5py
from AeViz.utils.utils import check_existence, progressBar, checkpoints
from AeViz.utils.file_utils import save_hdf


def calculate_profile(simulation, profile, save_checkpoints, **kwargs):
    """
    Returns the radial profile of the desired quantity. Quantity names
    must be the same as the ones in the simulation class.
    """
    if profile in ['Rossby_number', 'Ye', 'temperature', 'rho',
                   'entropy', 'convective_flux', 'gas_pressure']:
        return read_profile(simulation, profile, save_checkpoints)
    elif profile == 'BV_frequency':
        if 'mode' not in kwargs:
            kwargs['mode'] = 1
        if kwargs['mode'] == 1:
            return read_profile(simulation, 'BV_frequency', save_checkpoints)
        elif kwargs['mode'] == 2:
            return derive_profile(simulation, 'BV_frequency',
                                  **kwargs)
    else:
        return derive_profile(simulation, profile, **kwargs)

def read_profile(simulation, profile, save_checkpoints):
    """
    Reads the radial profile saved in the profiles.h5 file.
    """
    if check_existence(simulation, 'profiles.h5'):
        data = h5py.File(os.path.join(simulation.storage_path, 'profiles.h5'), 'r')
        if len(data['time'][...])  == len(simulation.hdf_file_list):
            t, pr = data['time'][...], data['profiles/' + profile][...]
            data.close()
            return t, simulation.cell.radius(simulation.ghost), pr
        else:
            derive_profiles(simulation, data, save_checkpoints)
            return read_profile(simulation, profile, save_checkpoints)
    else:
        derive_profiles(simulation, None, save_checkpoints)
        return read_profile(simulation, profile, save_checkpoints)

def derive_profile(simulation, profile, **kwargs):
    """
    Calculates a single profile and returns it.
    """
    qt = getattr(simulation, profile)
    radius = len(simulation.cell.radius(simulation.ghost))
    dOmega = simulation.cell.dOmega(simulation.ghost)
    
    for (file, progress) in zip(simulation.hdf_file_list,
                         range(len(simulation.hdf_file_list))):
        try:
            time = np.concatenate((time, simulation.time(file)))
            qt_local = qt(file, **kwargs)
            if qt_local.shape[-1] != radius:
                qt_av = np.zeros((radius, qt_local.shape[-1]))
                for i in range(qt_local.shape[-1]):
                    qt_av[:, i] = function_average(qt_local[..., i], 
                                                   simulation.dim, 'Omega',
                                                   dOmega)
            else:                                     
                qt_av = function_average(qt(file), simulation.dim, 'Omega',
                                            dOmega)
            profiles = np.concatenate((profiles, qt_av[..., None]), axis=-1)
        except:
            time = simulation.time(file)
            qt_local = qt(file, **kwargs)
            if qt_local.shape[-1] != radius:
                qt_av = np.zeros((radius, qt_local.shape[-1]))
                for i in range(qt_local.shape[-1]):
                    qt_av[:, i] = function_average(qt_local[..., i], 
                                                   simulation.dim, 'Omega',
                                                   dOmega)
            else:                                     
                qt_av = function_average(qt(file), simulation.dim, 'Omega',
                                            dOmega)
            profiles = qt_av[..., None]
        progressBar(progress, len(simulation.hdf_file_list),
                    'Calculating profile')
    if profiles.ndim == 3:
        profiles = profiles.swapaxes(-2, -1)
    return time, simulation.cell.radius(simulation.ghost), profiles

def derive_profiles(simulation, data, save_checkpoints):
    """
    Calculates and saves the pressure, temperature, Ye, entropy, 
    density, Rossby number, BV frequency, and convective flux radial
    profiles in an hdf file.
    """
    if data is None:
        start_point = 0
    else:
        time = data['time'][...]
        start_point = len(time)
        profiles = {'BV_frequency': data['profiles/BV_frequency'][...],
                    'Rossby_number': data['profiles/Rossby_number'][...],
                    'Ye': data['profiles/Ye'][...],
                    'temperature': data['profiles/temperature'][...],
                    'rho': data['profiles/rho'][...],
                    'entropy': data['profiles/entropy'][...],
                    'convective_flux': data['profiles/convective_flux'][...],
                    'gas_pressure': data['profiles/gas_pressure'][...]}
        data.close()
    
    if (checkpoints[simulation.dim] == False) or (not save_checkpoints):
        checkpoint = len(simulation.hdf_file_list)
    else:
        checkpoint = checkpoints[simulation.dim]
    dOmega = simulation.cell.dOmega(simulation.ghost)
    checkpoint_index = 0
    progress_index = 0
    total_points = len(simulation.hdf_file_list) - start_point
    for file in simulation.hdf_file_list[start_point:]:
        t_file = simulation.time(file, True)
        T_av = function_average(simulation.temperature(file), simulation.dim,
                                'Omega', dOmega)[..., None]
        rho_av = function_average(simulation.rho(file), simulation.dim,
                                  'Omega', dOmega)[..., None]
        S_av = function_average(simulation.entropy(file), simulation.dim,
                                'Omega', dOmega)[..., None]
        P_av = function_average(simulation.gas_pressure(file), simulation.dim,
                                'Omega', dOmega)[..., None]
        Ye_av = function_average(simulation.Ye(file), simulation.dim, 'Omega',
                                 dOmega)[..., None]
        BV_av = function_average(simulation.BV_frequency(file), simulation.dim,
                                 'Omega', dOmega)[..., None]
        if simulation.dim == 1:
            Ro_av = np.zeros(BV_av.shape)
            Fc_av = np.zeros(BV_av.shape)
        else:
            Ro_av = function_average(simulation.Rossby_number(file), 
                                    simulation.dim, 'Omega', dOmega)[..., None]
            
            Fc_av = simulation.convective_flux(file)[..., None]
        

        try:
            time = np.concatenate((time, t_file))
            profiles = {
                'BV_frequency': np.concatenate((profiles['BV_frequency'],
                                                BV_av), axis=-1),
                'Rossby_number': np.concatenate((profiles['Rossby_number'],
                                                 Ro_av), axis=-1),
                'Ye': np.concatenate((profiles['Ye'], Ye_av), axis=-1),
                'temperature': np.concatenate((profiles['temperature'],
                                               T_av), axis=-1),
                'rho': np.concatenate((profiles['rho'], rho_av), axis=-1),
                'entropy': np.concatenate((profiles['entropy'], S_av),
                                          axis=-1),
                'convective_flux': np.concatenate((profiles['convective_flux'],
                                                    Fc_av), axis=-1),
                'gas_pressure': np.concatenate((profiles['gas_pressure'],
                                                    P_av), axis=-1)
            }
        except:
            time = t_file
            profiles = {
                'BV_frequency': BV_av,
                'Rossby_number': Ro_av,
                'Ye': Ye_av,
                'temperature': T_av,
                'rho': rho_av,
                'entropy': S_av,
                'convective_flux': Fc_av,
                'gas_pressure': P_av
            }
        if checkpoint_index >= checkpoint:
            checkpoint_index = 0
            print('Saving checkpoint...')
            save_hdf(os.path.join(simulation.storage_path, 'profiles.h5'),
                     ['time', 'profiles'],  [time, profiles])
        
        progressBar(progress_index, total_points, 'Calculating profiles')
        progress_index += 1
        checkpoint_index += 1
    save_hdf(os.path.join(simulation.storage_path, 'profiles.h5'),
             ['time', 'profiles'],  [time, profiles])
    print('Profiles saved.')
    
    
        
        
        
        