from AeViz.utils.physics.masses_energies_utils import (standard_mass_energy, 
                                               gain_region_mass_energy,
                                               PNS_mass_energy,
                                               unbound_mass_energy,
                                               mass_flux)
from AeViz.utils.utils import (check_existence, progressBar, checkpoints)
from AeViz.utils.files.file_utils import save_hdf, create_series
from AeViz.grid.grid import grid
import os, h5py
import numpy as np
from AeViz.units.aerray import aerray
from AeViz.units import u


def calculate_masses_energies(simulation, save_checkpoints=True):
    if check_existence(simulation, 'masses_energies.h5'):
        time, mdot, inner_me, gain_me, PNS_me, unb_me, nuc_me, \
            processed_hdf = \
            read_masses_energies(simulation)
        ## Retrocompatibility option
        if processed_hdf is None:
            if len(simulation.hdf_file_list) == len(time):
                save_hdf(os.path.join(simulation.storage_path,
                    'masses_energies.h5'),
                    ['time', 'mass_flux', 'innercore', 'gain_region', 'PNS',
                      'unbound', 'PNS_core', 'processed'],
                    [time, mdot, inner_me, gain_me, PNS_me, unb_me, nuc_me,
                    simulation.hdf_file_list])
                return create_series(time, mdot, inner_me, gain_me, PNS_me,
                                     unb_me, nuc_me)
            else:
                start_point = 0
                time = 0
                mdot = 0
                nuc_me = 0
                inner_me = 0
                gain_me = 0
                PNS_me = 0
                unb_data = 0
                processed_hdf = []
        elif processed_hdf[-1].decode("utf-8") == simulation.hdf_file_list[-1]:
            return create_series(time, mdot, inner_me, gain_me, PNS_me, unb_me,
                                 nuc_me)
        else:
            start_point = len(processed_hdf)
            processed_hdf = [ff.decode("utf-8") for ff in processed_hdf]
            print('Checkpoint found for the mass and energy file, starting' \
                  ' from checkpoint.\nPlease wait...')
    else:
        start_point = 0
        processed_hdf = []
        print('No checkpoint found for the mass and energy file, starting' \
              ' from the beginning.\nPlease wait...')
    if (checkpoints[simulation.dim] == False) or (not save_checkpoints):
        checkpoint = len(simulation.hdf_file_list)
    else:
        checkpoint = checkpoints[simulation.dim]
    ## Get the radii
    innercore_radius, igcells  = simulation.innercore_radius(rad='full')
    PNS_core_radius, pcgcells  = simulation.PNS_nucleus_radius(rad='full')
    shock_radius, sgcells      = simulation.shock_radius(rad='full')
    gain_radius, ggcells       = simulation.gain_radius(rad='full')
    PNS_radius, pgcells        = simulation.PNS_radius(rad='full')
    innercore_radius = innercore_radius.data
    PNS_core_radius = PNS_core_radius.data
    shock_radius = shock_radius.data
    gain_radius = gain_radius.data
    PNS_radius = PNS_radius.data
    ## Get the grid
    if simulation.dim == 1:
        gr = grid(1, simulation.cell.radius(simulation.ghost))
        X, Y, Z = (gr.cartesian_grid(), 0, 0)
    elif simulation.dim == 2:
        gr = grid(2, simulation.cell.radius(simulation.ghost), 
                  simulation.cell.theta(simulation.ghost))
        X, Z = gr.cartesian_grid()
        Y = 0 * X.unit
    else:
        gr = grid(3, simulation.cell.radius(simulation.ghost), 
                  simulation.cell.theta(simulation.ghost),
                  simulation.cell.phi(simulation.ghost))
        X, Y, Z = gr.cartesian_grid()

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
        in_data = standard_mass_energy(simulation, file,
                                        innercore_radius[..., findex], 
                                        igcells, dV)
        nuc_data = standard_mass_energy(simulation, file,
                                        PNS_core_radius[..., findex], 
                                        pcgcells, dV)
        gr_data = gain_region_mass_energy(simulation, file,
                                          shock_radius[..., findex], sgcells,
                                          gain_radius[..., findex], ggcells,
                                          dV)
        PNS_data = PNS_mass_energy(simulation, file, PNS_radius[..., findex],
                                   pgcells, dV, (X, Y, Z), gr)
        unb_data = unbound_mass_energy(simulation, file, dV)
        try:
            time = np.concatenate((time, simulation.time(file)))
            mdot = np.concatenate((mdot, mass_flux(simulation, file, 
                                                    dOmega, radius_index)))
            ## PNS CORE
            nuc_me['mass'] = np.concatenate((nuc_me['mass'], nuc_data[0]))
            
            nuc_me['kinetic_ene'] = np.concatenate((nuc_me['kinetic_ene'],
                                                nuc_data[1]))
            nuc_me['magnetic_ene'] = np.concatenate((nuc_me['magnetic_ene'],
                                                nuc_data[2]))
            nuc_me['rotational_ene'] = np.concatenate((
                nuc_me['rotational_ene'], nuc_data[3]))
            nuc_me['grav_ene'] = np.concatenate((
                nuc_me['grav_ene'], nuc_data[4]))
            nuc_me['total_ene'] = np.concatenate((
                nuc_me['total_ene'], nuc_data[5]))
            nuc_me['T_W'] = np.concatenate((
                nuc_me['T_W'], nuc_data[6]))
            ## INNERCORE
            inner_me['mass'] = np.concatenate((inner_me['mass'], in_data[0]))
            
            inner_me['kinetic_ene'] = np.concatenate((inner_me['kinetic_ene'],
                                                in_data[1]))
            inner_me['magnetic_ene'] = np.concatenate((inner_me['magnetic_ene'],
                                                in_data[2]))
            inner_me['rotational_ene'] = np.concatenate((
                inner_me['rotational_ene'], in_data[3]))
            inner_me['grav_ene'] = np.concatenate((
                inner_me['grav_ene'], in_data[4]))
            inner_me['total_ene'] = np.concatenate((
                inner_me['total_ene'], in_data[5]))
            inner_me['T_W'] = np.concatenate((
                inner_me['T_W'], in_data[6]))
            ## GAIN
            gain_me['mass'] = np.concatenate((gain_me['mass'], gr_data[0]))
            gain_me['heating_ene'] = np.concatenate((gain_me['heating_ene'],
                                                gr_data[1]))
            ## PNS
            PNS_me['mass'] = np.concatenate((PNS_me['mass'], PNS_data[0]))
            PNS_me['kinetic_ene'] = np.concatenate((PNS_me['kinetic_ene'],
                                                PNS_data[1]))
            PNS_me['magnetic_ene'] = np.concatenate((PNS_me['magnetic_ene'],
                                                PNS_data[2]))
            PNS_me['rotational_ene'] = np.concatenate((
                PNS_me['rotational_ene'], PNS_data[3]))
            PNS_me['grav_ene'] = np.concatenate((
                PNS_me['grav_ene'], PNS_data[4]))
            PNS_me['total_ene'] = np.concatenate((
                PNS_me['total_ene'], PNS_data[5]))
            PNS_me['convective_ene'] = np.concatenate((
                PNS_me['convective_ene'], PNS_data[6]))
            PNS_me['L']['Lx'] = np.concatenate((PNS_me['L']['Lx'],
                                                PNS_data[7]))
            PNS_me['L']['Ly'] = np.concatenate((PNS_me['L']['Ly'],
                                                PNS_data[8]))
            PNS_me['L']['Lz'] = np.concatenate((PNS_me['L']['Lz'],
                                                PNS_data[9]))
            PNS_me['L']['L_tot'] = np.concatenate((PNS_me['L']['L_tot'],
                                                PNS_data[10]))
            ## EXPLOSION
            unb_me['mass'] = np.concatenate((unb_me['mass'], unb_data[0]))
            unb_me['energy'] = np.concatenate((unb_me['energy'], unb_data[1]))
            unb_me['kinetic_ene'] = np.concatenate((unb_me['kinetic_ene'],
                                                unb_data[2]))
            unb_me['magnetic_ene'] = np.concatenate((unb_me['magnetic_ene'],
                                                unb_data[3]))
        except Exception as e:
            print(e)
            time = simulation.time(file)
            mdot = mass_flux(simulation, file, dOmega, radius_index)
            inner_me = {
                'mass': in_data[0],
                'kinetic_ene': in_data[1],
                'magnetic_ene': in_data[2],
                'rotational_ene': in_data[3],
                'grav_ene': in_data[4],
                'total_ene': in_data[5],
                'T_W': in_data[6]
            }
            nuc_me = {
                'mass': nuc_data[0],
                'kinetic_ene': nuc_data[1],
                'magnetic_ene': nuc_data[2],
                'rotational_ene': nuc_data[3],
                'grav_ene': nuc_data[4],
                'total_ene': nuc_data[5],
                'T_W': nuc_data[6]
            }
            gain_me = {
                'mass': gr_data[0],
                'heating_ene': gr_data[1]
            }
            PNS_me = {
                'mass': PNS_data[0],
                'kinetic_ene': PNS_data[1],
                'magnetic_ene': PNS_data[2],
                'rotational_ene': PNS_data[3],
                'grav_ene': PNS_data[4],
                'total_ene': PNS_data[5],
                'convective_ene': PNS_data[6], 
                'L': {
                    'Lx': PNS_data[7],
                    'Ly': PNS_data[8],
                    'Lz': PNS_data[9],
                    'L_tot': PNS_data[10]
                }
            }
            unb_me = {
                'mass': unb_data[0],
                'energy': unb_data[1],
                'kinetic_ene': unb_data[2],
                'magnetic_ene': unb_data[3]
            }
        processed_hdf.append(file)
        if (check_index >= checkpoint) and save_checkpoints:
            print('Checkpoint reached, saving...\n')
            save_hdf(os.path.join(simulation.storage_path,
                    'masses_energies.h5'),
                     ['time', 'mass_flux', 'innercore', 'gain_region', 'PNS',
                      'unbound', 'PNS_core', 'processed'],
                     [time, mdot, inner_me, gain_me, PNS_me, unb_me, nuc_me,
                      processed_hdf])
            
            check_index = 0
        check_index += 1
        progress_index += 1
        findex += 1
    print('Computation completed, saving...')
    save_hdf(os.path.join(simulation.storage_path, 'masses_energies.h5'),
                    ['time', 'mass_flux', 'innercore', 'gain_region', 'PNS',
                    'unbound', 'PNS_core', 'processed'],
                    [time, mdot, inner_me, gain_me, PNS_me, unb_me, nuc_me,
                     processed_hdf])
    time, mdot, inner_me, gain_me, PNS_me, unb_me, nuc_me, _ = \
        read_masses_energies(simulation)
    return create_series(time, mdot, inner_me, gain_me, PNS_me, unb_me, nuc_me)

def read_masses_energies(simulation):
    masses_energies_data = h5py.File(os.path.join(simulation.storage_path, 
                                            'masses_energies.h5'), 'r')
    data = [
        aerray(masses_energies_data['time'][...], u.s, 'time',
               r'$t-t_\mathrm{b}$', None, [-0.005,
                                           masses_energies_data['time'][-1]]),
        aerray(masses_energies_data['mass_flux'][...], u.M_sun/u.s, 'M_accr',
               r'$\dot{M}_{500km}$', None, [0, 3]),
        {
            'mass': aerray(masses_energies_data['innercore/mass'][...], u.M_sun,
                           'M_innercore', r'$M_\mathrm{innercore}$', None,
                           [0,2], False),
            'kinetic_ene': aerray(
                masses_energies_data['innercore/kinetic_ene'][...], u.erg,
                'kin_innercore', r'$E_\mathrm{kin,inn}$', None, [1e47, 1e53],
                True),
            'magnetic_ene': aerray(
                masses_energies_data['innercore/magnetic_ene'][...], u.erg,
                'mag_innercore', r'$E_\mathrm{mag,inn}$', None, [1e47, 1e52], True),
            'rotational_ene': aerray(
                masses_energies_data['innercore/rotational_ene'][...], u.erg,
                'rot_innercore',  r'$E_\mathrm{rot,inn}$', None, [1e47, 1e52], True),
            'grav_ene': aerray(masses_energies_data['innercore/grav_ene'][...],
                               u.erg, 'grav_innercore', r'$E_\mathrm{grav,inn}$',
                               None, [-1e54, 1e48], True),
            'total_ene': aerray(masses_energies_data['innercore/total_ene'][...],
                                u.erg, 'tot_innercore', r'$E_\mathrm{tot,inn}$',
                                None, [1e47, 1e53], True),
            'T_W': aerray(masses_energies_data['innercore/T_W'][...],
                          u.dimensionless_unscaled, 'T/W_innercore',
                          r'$T/|W|_\mathrm{inn}$', None, [0, 0.05], False)
        },
        {
            'mass': aerray(masses_energies_data['gain_region/mass'][...],
                           u.M_sun, 'gain_mass', r'$M_\mathrm{gain}$', None,
                           [0, 1], False),
            'heating_ene': aerray(
                masses_energies_data['gain_region/heating_ene'][...], u.erg,
                'gain_energy', r'$Q_\nu$', None, [-1e51, 3e51])
        },
        {
            'mass': aerray(masses_energies_data['PNS/mass'][...], u.M_sun,
                           'PNS_mass', r'$M_\mathrm{PNS}$', None, [0, 2]),
            'kinetic_ene': aerray(masses_energies_data['PNS/kinetic_ene'][...],
                                  u.erg, 'PNS_kin', r'$E_\mathrm{kin,PNS}$',
                                  None, [1e47, 1e52], True),
            'magnetic_ene': aerray(masses_energies_data['PNS/magnetic_ene'][...],
                                   u.erg, 'PNS_mag', r'$E_\mathrm{mag,PNS}$',
                                   None, [1e47, 1e51], True),
            'rotational_ene': aerray(
                masses_energies_data['PNS/rotational_ene'][...], u.erg, 'PNS_rot',
                r'$E_\mathrm{rot,PNS}$', None, [1e47, 1e51], True),
            'grav_ene': aerray(masses_energies_data['PNS/grav_ene'][...], u.erg,
                               'PNS_grav', r'$E_\mathrm{grav,PNS}$', None,
                               [-1e54, -1e49], True),
            'total_ene': aerray(masses_energies_data['PNS/total_ene'][...],
                                u.erg, 'PNS_tot', r'$E_\mathrm{tot,PNS}$', None,
                                [1e47, 1e52], True),                                
            'convective_ene': aerray(
                masses_energies_data['PNS/convective_ene'][...], u.erg,
                'PNS_conv', r'$E_\mathrm{conv,PNS}$', None, [1e47, 1e51], True),
            'L': {
                'Lx': aerray(masses_energies_data['PNS/L/Lx'][...],
                             u.erg*u.s, 'PNS_Jx', r'$J_\mathrm{PNS,x}$', None,
                             [-1e47, 1e47], False),
                'Ly': aerray(masses_energies_data['PNS/L/Ly'][...],
                             u.erg*u.s, 'PNS_Jy', r'$J_\mathrm{PNS,y}$', None,
                             [-1e47, 1e47], False),
                'Lz': aerray(masses_energies_data['PNS/L/Lz'][...],
                             u.erg*u.s, 'PNS_Jz', r'$J_\mathrm{PNS,z}$', None,
                             [0, 1e49], False),
                'L_tot': aerray(masses_energies_data['PNS/L/L_tot'][...],
                             u.erg*u.s, 'PNS_Jtot', r'$J_\mathrm{PNS,tot}$', None,
                             [0, 1e49], False)
            }
        },
        {
            'mass': aerray(masses_energies_data['unbound/mass'][...], u.M_sun,
                           'unbound_mass', r'$M_\mathrm{unb}$', None, [1e-4, 1],
                           True),
            'energy': aerray(masses_energies_data['unbound/energy'][...], u.erg,
                             'expl_ene', r'$E_\mathrm{expl}$', None, [1e47, 1e51],
                             True),
            'kinetic_ene': aerray(
                masses_energies_data['unbound/kinetic_ene'][...], u.erg,
                'expl_kin', r'$E_\mathrm{expl,kin}$', None, [1e47, 1e52],True),
            'magnetic_ene': aerray(
                masses_energies_data['unbound/magnetic_ene'][...], u.erg,
                'expl_mag', r'$E_\mathrm{expl,mag}$', None, [1e47, 1e52], True)
        },
        {
            'mass': aerray(masses_energies_data['PNS_core/mass'][...], u.M_sun,
                           'PNS_core_mass', r'$M_\mathrm{PNS core}$', None,
                           [0, 2], False),
            'kinetic_ene': aerray(
                masses_energies_data['PNS_core/kinetic_ene'][...], u.erg,
                'PNS_core_kin', r'$E_\mathrm{kin,PNS core}$', None, [1e47, 1e53],
                False),
            'magnetic_ene': aerray(
                masses_energies_data['PNS_core/magnetic_ene'][...], u.erg,
                'PNS_core_mag', r'$E_\mathrm{mag,PNS core}$', None, [1e47, 1e53],
                False),
            'rotational_ene': aerray(
                masses_energies_data['PNS_core/rotational_ene'][...], u.erg,
                'PNS_core_rot', r'$E_\mathrm{rot,PNS core}$', None, [1e47, 1e53],
                False),
            'grav_ene': aerray(masses_energies_data['PNS_core/grav_ene'][...],
                               u.erg, 'PNS_core_grav',
                               r'$E_\mathrm{grav,PNS core}$', None,
                               [-1e54, -1e48], False),
            'total_ene': aerray(masses_energies_data['PNS_core/total_ene'][...],
                                u.erg, 'PNS_core_tot', r'$E_\mathrm{tot,PNS core}$',
                                None, [1e47, 1e53], False),
            'T_W': aerray(masses_energies_data['PNS_core/T_W'][...],
                          u.dimensionless_unscaled, 'T/W_core',
                          r'$T/|W|_\mathrm{PNS core}$', None, [0, 0.05], False)
        }
    ]
    if 'processed' in masses_energies_data:
        data.append(masses_energies_data['processed'][...])
    else:
        data.append(None)
    return data