import numpy as np
from AeViz.utils.math_utils import function_average
import os, h5py
from AeViz.utils.utils import check_existence, progressBar, checkpoints
from AeViz.utils.files.file_utils import save_hdf
from AeViz.units.aeseries import aerray, aeseries
from AeViz.units import u
from AeViz.utils.files.string_utils import merge_strings


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
            return derive_profile(simulation, 'BV_frequency', mode = 2)
    else:
        return derive_profile(simulation, profile)

def read_profile(simulation, profile, save_checkpoints):
    """
    Reads the radial profile saved in the profiles.h5 file.
    """
    if check_existence(simulation, 'profiles.h5'):
        data = h5py.File(os.path.join(simulation.storage_path, 'profiles.h5'), 'r')
        if not 'processed' in data.keys():
            data.close()
            data = h5py.File(os.path.join(simulation.storage_path,
                                          'profiles.h5'), 'r+')
            data.create_dataset('processed',
                                data=simulation.hdf_file_list[:len(data['time'][...])])
            data.close()
            data = h5py.File(os.path.join(simulation.storage_path,
                                          'profiles.h5'), 'r')
        if data['processed'][-1].decode("utf-8") == simulation.hdf_file_list[-1]:
            t, pr = data['time'][...], data['profiles/' + profile][...]
            data.close()
            return make_series(t, simulation.cell.radius(simulation.ghost), pr,
                               profile)
        else:
            derive_profiles(simulation, data, save_checkpoints)
            return read_profile(simulation, profile, save_checkpoints)
    else:
        derive_profiles(simulation, None, save_checkpoints)
        return read_profile(simulation, profile, save_checkpoints)

def derive_profile(simulation, profile):
    """
    Calculates a single profile and returns it.
    """
    qt = getattr(simulation, profile)
    radius = len(simulation.cell.radius(simulation.ghost))
    dOmega = simulation.cell.dOmega(simulation.ghost)
    zeroth_step = qt(0)
    nm, lb, lm, = zeroth_step.name, zeroth_step.label, zeroth_step.limits
    lg, cm = zeroth_step.log, zeroth_step.cmap
    nm = nm + '_profile'
    lb = merge_strings(r'$\langle $', lb, r'$\rangle_\Omega$')
    for (file, progress) in zip(simulation.hdf_file_list,
                         range(len(simulation.hdf_file_list))):
        try:
            time = np.concatenate((time, simulation.time(file)))
            qt_local = qt(file)
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
            qt_local = qt(file)
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
    profiles.set(name=nm, label=lb, limits=lm, cmap=cm, log=lg)
    time = aerray(time, u.s, 'time', r'$t-t_\mathrm{b}$', None,
                  [-0.005, time[-1]])
    return aeseries(profiles, time=time,
                    radius=simulation.cell.radius(simulation.ghost)) 

def derive_profiles(simulation, data, save_checkpoints):
    """
    Calculates and saves the pressure, temperature, Ye, entropy, 
    density, Rossby number, BV frequency, and convective flux radial
    profiles in an hdf file.
    """
    if data is None:
        start_point = 0
        processed_hdf = []
    else:
        time = data['time'][...] * u.s
        start_point = len(data['processed'][...])
        processed_hdf = [ff.decode("utf-8") for ff in data['processed'][...]]
        print('Checkpoint found. Starting from timestep', start_point)
        profiles = {'BV_frequency': data['profiles/BV_frequency'][...] * u.s ** (-2),
                    'Rossby_number': data['profiles/Rossby_number'][...] * \
                    u.dimensionless_unscaled,
                    'Ye': data['profiles/Ye'][...] * u.dimensionless_unscaled,
                    'temperature': data['profiles/temperature'][...] * u.MeV,
                    'rho': data['profiles/rho'][...] * u.g / u.cm ** 3,
                    'entropy': data['profiles/entropy'][...] * u.kBol / u.bry,
                    'convective_flux': data['profiles/convective_flux'][...] * \
                    u.erg / u.s / u.cm ** 2,
                    'gas_pressure': data['profiles/gas_pressure'][...] * u.Ba}
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
        except Exception as e:
            print(e)
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
        processed_hdf.append(file)
        if checkpoint_index >= checkpoint:
            checkpoint_index = 0
            print('Saving checkpoint...')
            save_hdf(os.path.join(simulation.storage_path, 'profiles.h5'),
                     ['time', 'profiles', 'processed'],
                     [time, profiles, processed_hdf])
        
        progressBar(progress_index, total_points, 'Calculating profiles...')
        progress_index += 1
        checkpoint_index += 1
    save_hdf(os.path.join(simulation.storage_path, 'profiles.h5'),
             ['time', 'profiles', 'processed'],  [time, profiles,
                                                  processed_hdf])
    print('Profiles saved.')
    
def make_series(time, radius, prof, name):
    """
    Returns a series of profiles for a given time.
    """
    t = aerray(time, u.s, 'time', r'$t-t_\mathrm{b}$', None, [-0.005, time[-1]])
    if name == 'Rossby_number':
        pr = aerray(prof, u.dimensionless_unscaled, 'Rossby_number_profile',
                    r'$\langle Ro\rangle_\Omega$', 'RdYlBu_r', [-1e-4, 1e-4],
                    True)
    elif name == 'Ye':
        pr = aerray(prof, u.dimensionless_unscaled, 'Ye_profile',
                    r'$\langle Y_\mathrm{e}\rangle_\Omega$', 'gist_rainbow',
                    [0.0, 0.5], False)
    elif name == 'temperature':
        pr = aerray(prof, u.MeV, 'temperature_profile',
                    r'$\langle T\rangle_\Omega$', 'inferno', [0.0, 40], False)
    elif name == 'rho':
        pr = aerray(prof, u.g / u.cm ** 3, 'rho_profile',
                    r'$\langle \rho\rangle_\Omega$', 'viridis',  [1e4, 1e15],
                    True)
    elif name == 'entropy':
        pr = aerray(prof, u.kBol / u.bry, 'rho_profile',
                    r'$\langle s\rangle_\Omega$', 'gist_rainbow_r',  [1.5, 15],
                    False)
    elif name == 'convective_flux':
        pr = aerray(prof, u.erg / u.s / u.cm ** 2, 'Fconv_profile',
                    r'$\langle F_\mathrm{conv}\rangle_\Omega$', 'RdYlGn_r',
                    [-1e40, 1e40], True)
    elif name == 'gas_pressure':
        pr = aerray(prof, u.Ba, 'gas_pressure_profile',
                    r'$\langle P_\mathrm{gas}\rangle_\Omega$',
                    'gist_rainbow_r',  [1e25, 1e34], True)
    elif name == 'BV_frequency':
        pr = aerray(prof, u.s ** (-2), 'BV_profile',
                    r'$\langle \omega_\mathrm{BV}\rangle_\Omega$', 'coolwarm',  [-1e6, 1e6],
                    True)
    return aeseries(pr, time=t, radius=radius)