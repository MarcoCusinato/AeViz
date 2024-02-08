from AeViz.utils.masses_energies_utils import (innercore_mass_energy, 
                                               gain_region_mass_energy,
                                               PNS_mass_energy,
                                               unbound_mass_energy,
                                               mass_flux)
from AeViz.utils.utils import (check_existence, progressBar, checkpoints)
from AeViz.utils.file_utils import save_hdf
import os, h5py
import numpy as np


def calculate_masses_energies(simulation, save_checkpoints=True):
    if check_existence(simulation, 'masses_energies.h5'):
        time, mdot, inner_me, gain_me, PNS_me, unb_me = \
            read_masses_energies(simulation)
        if len(simulation.hdf_file_list) == len(time):
            return time, mdot, inner_me, gain_me, PNS_me, unb_me
        else:
            start_point = len(time)
            print('Checkpoint found for the mass and energy file, starting' \
                  ' from checkpoint.\nPlease wait...')
    else:
        start_point = 0
        print('No checkpoint found for the mass and energy file, starting' \
              ' from the beginning.\nPlease wait...')
    if (checkpoints[simulation.dim] == False) or (not save_checkpoints):
        checkpoint = len(simulation.hdf_file_list)
    else:
        checkpoint = checkpoints[simulation.dim]
    ## Get the radii
    _, innercore_radius, _, _, _, igcells  = simulation.innercore_radius()
    _, shock_radius, _, _, _, sgcells = simulation.shock_radius()
    _, gain_radius, _, _, _, ggcells = simulation.gain_radius()
    _, PNS_radius, _, _, _, pgcells = simulation.PNS_radius()

    ## Get the volume elements
    dV = simulation.cell.dVolume_integration(simulation.ghost)
    dOmega = simulation.cell.dOmega(simulation.ghost)
    radius_index = np.argmax(simulation.cell.radius(simulation.ghost) >= 5e7)
    ## Indices
    findex = start_point
    check_index = 0
    progress_index = 0
    total_points = len(simulation.hdf_file_list) - start_point
    for file in simulation.hdf_file_list[start_point:]:
        progressBar(progress_index, total_points, suffix='Computing...')
        in_data = innercore_mass_energy(simulation, file,
                                        innercore_radius[..., findex], 
                                        igcells, dV)
        gr_data = gain_region_mass_energy(simulation, file,
                                          shock_radius[..., findex], sgcells,
                                          gain_radius[..., findex], ggcells,
                                          dV)
        PNS_data = PNS_mass_energy(simulation, file, PNS_radius[..., findex],
                                   pgcells, dV)
        unb_data = unbound_mass_energy(simulation, file, dV)
        try:
            time = np.concatenate((time, simulation.time(file)))
            mdot = np.concatenate((mdot, np.array([mass_flux(simulation, file, 
                                                    dOmega, radius_index)])))
            inner_me['mass'] = np.concatenate((inner_me['mass'], 
                                                np.array([in_data[0]])))
            
            inner_me['kinetic_ene'] = np.concatenate((inner_me['kinetic_ene'],
                                                np.array([in_data[1]])))
            inner_me['magnetic_ene'] = np.concatenate((inner_me['magnetic_ene'],
                                                np.array([in_data[2]])))
            inner_me['rotational_ene'] = np.concatenate((
                inner_me['rotational_ene'], np.array([in_data[3]])))
            inner_me['grav_ene'] = np.concatenate((
                inner_me['grav_ene'], np.array([in_data[4]])))
            inner_me['total_ene'] = np.concatenate((
                inner_me['total_ene'], np.array([in_data[5]])))
            inner_me['T_W'] = np.concatenate((
                inner_me['T_W'], np.array([in_data[6]])))
            
            gain_me['mass'] = np.concatenate((gain_me['mass'], 
                                                np.array([gr_data[0]])))
            gain_me['heating_ene'] = np.concatenate((gain_me['heating_ene'],
                                                np.array([gr_data[1]])))
            PNS_me['mass'] = np.concatenate((PNS_me['mass'],
                                                np.array([PNS_data[0]])))
            PNS_me['kinetic_ene'] = np.concatenate((PNS_me['kinetic_ene'],
                                                np.array([PNS_data[1]])))
            PNS_me['magnetic_ene'] = np.concatenate((PNS_me['magnetic_ene'],
                                                np.array([PNS_data[2]])))
            PNS_me['rotational_ene'] = np.concatenate((
                PNS_me['rotational_ene'], np.array([PNS_data[3]])))
            PNS_me['grav_ene'] = np.concatenate((
                PNS_me['grav_ene'], np.array([PNS_data[4]])))
            PNS_me['total_ene'] = np.concatenate((
                PNS_me['total_ene'], np.array([PNS_data[5]])))
            PNS_me['convective_ene'] = np.concatenate((
                PNS_me['convective_ene'], np.array([PNS_data[6]])))
            unb_me['mass'] = np.concatenate((unb_me['mass'],
                                                np.array([unb_data[0]])))
            unb_me['energy'] = np.concatenate((unb_me['energy'],
                                                np.array([unb_data[1]])))
        except:
            time = simulation.time(file)
            mdot = np.array([mass_flux(simulation, file, dOmega,
                                       radius_index)])
            inner_me = {
                'mass': np.array([in_data[0]]),
                'kinetic_ene': np.array([in_data[1]]),
                'magnetic_ene': np.array([in_data[2]]),
                'rotational_ene': np.array([in_data[3]]),
                'grav_ene': np.array([in_data[4]]),
                'total_ene': np.array([in_data[5]]),
                'T_W': np.array([in_data[6]])
            }
            gain_me = {
                'mass': np.array([gr_data[0]]),
                'heating_ene': np.array([gr_data[1]])
            }
            PNS_me = {
                'mass': np.array([PNS_data[0]]),
                'kinetic_ene': np.array([PNS_data[1]]),
                'magnetic_ene': np.array([PNS_data[2]]),
                'rotational_ene': np.array([PNS_data[3]]),
                'grav_ene': np.array([PNS_data[4]]),
                'total_ene': np.array([PNS_data[5]]),
                'convective_ene': np.array([PNS_data[6]])
            }
            unb_me = {
                'mass': np.array([unb_data[0]]),
                'energy': np.array([unb_data[1]])
            }
        if (check_index >= checkpoint) and save_checkpoints:
            print('Checkpoint reached, saving...\n')
            save_hdf(os.path.join(simulation.storage_path,
                    'masses_energies.h5'),
                     ['time', 'mass_flux', 'innercore', 'gain_region', 'PNS',
                      'unbound'],
                     [time, mdot, inner_me, gain_me, PNS_me, unb_me])
            
            check_index = 0
        check_index += 1
        progress_index += 1
        findex += 1
    print('Computation completed, saving...')
    save_hdf(os.path.join(simulation.storage_path, 'masses_energies.h5'),
                    ['time', 'mass_flux', 'innercore', 'gain_region', 'PNS',
                    'unbound'],
                    [time, mdot, inner_me, gain_me, PNS_me, unb_me])
    return time, mdot, inner_me, gain_me, PNS_me, unb_me
        

def read_masses_energies(simulation):
    masses_energies_data = h5py.File(os.path.join(simulation.storage_path, 
                                            'masses_energies.h5'), 'r')
    data = [
        masses_energies_data['time'][...],
        masses_energies_data['mass_flux'][...],
        {
            'mass': masses_energies_data['innercore/mass'][...],
            'kinetic_ene': masses_energies_data['innercore/kinetic_ene']\
                [...],
            'magnetic_ene': masses_energies_data['innercore/magnetic_ene']\
                [...],
            'rotational_ene': masses_energies_data['innercore/rotational_ene']\
                [...],
            'grav_ene': masses_energies_data['innercore/grav_ene']\
                [...],
            'total_ene': masses_energies_data['innercore/total_ene'][...],
            'T_W': masses_energies_data['innercore/T_W'][...]
        },
        {
            'mass': masses_energies_data['gain_region/mass'][...],
            'heating_ene': masses_energies_data['gain_region/heating_ene'][...]
        },
        {
            'mass': masses_energies_data['PNS/mass'][...],
            'kinetic_ene': masses_energies_data['PNS/kinetic_ene'][...],
            'magnetic_ene': masses_energies_data['PNS/magnetic_ene'][...],
            'rotational_ene': masses_energies_data['PNS/rotational_ene'][...],
            'grav_ene': masses_energies_data['PNS/grav_ene'][...],
            'total_ene': masses_energies_data['PNS/total_ene'][...],
            'convective_ene': masses_energies_data['PNS/convective_ene']\
                [...]
        },
        {
            'mass': masses_energies_data['unbound/mass'][...],
            'energy': masses_energies_data['unbound/energy'][...]
        }
    ]
    return data