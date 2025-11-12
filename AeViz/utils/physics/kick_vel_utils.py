import numpy as np
import h5py, os
from AeViz.utils.utils import (check_existence, progressBar, checkpoints)
from AeViz.utils.files.file_utils import save_hdf, create_series
from AeViz.units import u
from AeViz.units.constants import constants as c
from AeViz.grid.grid import grid
from AeViz.units.aerray import aerray

def integrate_momenta(time, PNS_mass, hydro_v, nue_v, nua_v, nux_v):
    dt = np.zeros(time.shape[0])
    dt[1:] = time[1:] - time[:-1]
    dt[0] = dt[1]
    dt = dt * time.unit
    hydro = {key: hydro_v[key] / PNS_mass for key in hydro_v.keys()}
    [hd.set(log=False, limits=[0, 7.5e7], name='kick_vel_hydro_'+key,
            label=r'$v_\mathrm{kick h,'+key+r'}$') for (key, hd) in
            hydro.items()]
    nue = {key: np.cumsum(comp * dt) / PNS_mass for (key, comp) in nue_v.items()}
    [nu.set(log=False, limits=[0, 7.5e7], name='kick_vel_nue_'+key,
            label=r'$v_{\mathrm{kick} \nu_\mathrm{e},'+key+r'}$') for (key, nu)
            in nue.items()]
    nua = {key: np.cumsum(comp * dt) / PNS_mass for (key, comp) in nua_v.items()}
    [nu.set(log=False, limits=[0, 7.5e7], name='kick_vel_nua_'+key,
            label=r'$v_{\mathrm{kick} \overline{\nu}_\mathrm{e},'+key+r'}$') 
            for (key, nu) in nua.items()]
    nux ={key: np.cumsum(comp * dt) / PNS_mass for (key, comp) in nux_v.items()}
    [nu.set(log=False, limits=[0, 7.5e7], name='kick_vel_nux_'+key,
            label=r'$v_{\mathrm{kick} \nu_\mathrm{x},'+key+r'}$') for (key, nu)
            in nux.items()]
    vhydro_tot = np.sqrt(hydro['x'] ** 2 + hydro['y'] ** 2 + hydro['z'] ** 2)
    vhydro_tot.set(log=False, limits=[0, 7.5e7], name='kick_vel_hydro_tot',
                   label=r'$v_\mathrm{kick h,tot}$')
    vnue_tot = np.sqrt(nue['x'] ** 2 + nue['y'] ** 2 + nue['z'] ** 2)
    vnue_tot.set(log=False, limits=[0, 7.5e7], name='kick_vel_nue_tot',
                    label=r'$v_{\mathrm{kick} \nu_\mathrm{e},tot}$')
    vnua_tot = np.sqrt(nua['x'] ** 2 + nua['y'] ** 2 + nua['z'] ** 2)
    vnua_tot.set(log=False, limits=[0, 7.5e7], name='kick_vel_nua_tot',
                    label=r'$v_{\mathrm{kick} \overline{\nu}_\mathrm{e},tot}$')
    vnux_tot = np.sqrt(nux['x'] ** 2 + nux['y'] ** 2 + nux['z'] ** 2)
    vnux_tot.set(log=False, limits=[0, 7.5e7], name='kick_vel_nux_tot',
                    label=r'$v_{\mathrm{kick} \nu_\mathrm{x},tot}$')
    vnu_tot = np.sqrt((nue['x'] + nua['x'] + nux['x']) ** 2 + 
                      (nue['y'] + nua['y'] + nux['y']) ** 2 + 
                      (nue['z'] + nua['z'] + nux['z']) ** 2)
    vnu_tot.set(log=False, limits=[0, 7.5e7], name='kick_vel_nu_tot',
                    label=r'$v_{\mathrm{kick} \nu,tot}$')
    vtot = np.sqrt((hydro['x'] + nue['x'] + nua['x'] + nux['x']) ** 2 +
                   (hydro['y'] + nue['y'] + nua['y'] + nux['y']) ** 2 +
                   (hydro['z'] + nue['z'] + nua['z'] + nux['z']) ** 2)
    vtot.set(log=False, limits=[0, 7.5e7], name='kick_vel_tot',
                    label=r'$v_\mathrm{kick,tot}$')
    return create_series(time, hydro, nue, nua, nux, vhydro_tot, vnue_tot,
                         vnua_tot, vnux_tot, vnu_tot, vtot)
    
def velocity_kick(simulation, file_name, PNS_radius, gcells, dV, dOmega, r400):
    """
    Calculates the velocity kick of the PNS for one timestep.
    """
    if simulation.time(file_name) <= 0 or simulation.dim == 1:
        return [[0.*u.cm/u.s*u.M_sun, 0.*u.cm/u.s*u.M_sun, 0.*u.cm/u.s*u.M_sun],
                [0.*u.M_sun*u.cm/u.s**2, 0.*u.M_sun*u.cm/u.s**2, 0.*u.M_sun*u.cm/u.s**2],
                [0.*u.M_sun*u.cm/u.s**2, 0.*u.M_sun*u.cm/u.s**2, 0.*u.M_sun*u.cm/u.s**2],
                [0.*u.M_sun*u.cm/u.s**2, 0.*u.M_sun*u.cm/u.s**2, 0.*u.M_sun*u.cm/u.s**2]]
    
    outer_radius = 1e10 * u.cm
    if outer_radius >= simulation.cell.radius(simulation.ghost)[-1]:
        outer_radius = simulation.cell.radius(simulation.ghost)[-2]
    
    mask = (simulation.cell.radius(simulation.ghost) >= \
        simulation.ghost.remove_ghost_cells_radii(PNS_radius,
                                                    simulation.dim,
                                                    **gcells)[..., None]) & \
            (simulation.cell.radius(simulation.ghost) <= outer_radius)
    rho = simulation.rho(file_name) * dV
    gh = grid(simulation.dim, simulation.cell.radius(simulation.ghost),
                              simulation.cell.theta(simulation.ghost),
                              simulation.cell.phi(simulation.ghost))
    vx, vy, vz = gh.velocity_sph_to_cart(
        simulation.radial_velocity(file_name) * rho,
        simulation.theta_velocity(file_name) * rho,
        simulation.phi_velocity(file_name) * rho)
    nue_flux, nua_flux, nux_flux = simulation.neutrino_momenta_grey(file_name)
    nue_flux = nue_flux[..., r400, :] * dOmega[..., None]
    nua_flux = nua_flux[..., r400, :] * dOmega[..., None]
    nux_flux = nux_flux[..., r400, :] * dOmega[..., None]
    nue_flux = list(gh.velocity_sph_to_cart(
                                        nue_flux[..., 0][..., None],
                                        nue_flux[..., 1][..., None],
                                        nue_flux[..., 2][..., None]
                                       )
                    )
    nua_flux = list(gh.velocity_sph_to_cart(
                                        nua_flux[..., 0][..., None],
                                        nua_flux[..., 1][..., None],
                                        nua_flux[..., 2][..., None]
                                        )
                    )
    nux_flux = list(gh.velocity_sph_to_cart(
                                        nux_flux[..., 0][..., None],
                                        nux_flux[..., 1][..., None],
                                        nux_flux[..., 2][..., None]
                                        )
                    )
    nue_flux = [-np.sum(comp).to(u.M_sun/u.s**3) \
                * ((400*u.km) ** 2 / c.c) for comp in nue_flux]
    nua_flux = [-np.sum(comp).to(u.M_sun/u.s**3) \
                * ((400*u.km) ** 2 / c.c) for comp in nua_flux]
    nux_flux = [-np.sum(comp).to(u.M_sun/u.s**3) \
                * (4 * (400*u.km) ** 2 / c.c) for comp in nux_flux]
    
    vz = np.sum(vz[mask]).to(u.cm/u.s*u.M_sun)
    if simulation.dim == 2:
        vy = 0.0 * u.cm/u.s*u.M_sun
        vx = 0.0 * u.cm/u.s*u.M_sun
        for nu in [nue_flux, nua_flux, nux_flux]:
            nu[0] = 0.0 * u.M_sun*u.cm/u.s**2
            nu[1] = 0.0 * u.M_sun*u.cm/u.s**2
    else:
        vy = np.sum(vy[mask]).to(u.cm/u.s*u.M_sun)
        vx = np.sum(vx[mask]).to(u.cm/u.s*u.M_sun)
    
    return [-vx, -vy, -vz], nue_flux, nua_flux, nux_flux

def calculate_kick(simulation, save_checkpoints=True, no_new=False):
    PNSmass = simulation.PNS_mass_ene(comp='mass').data
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
                return integrate_momenta(time, PNSmass, hydro_v, nu_flux['nue'],
                                         nu_flux['nua'], nu_flux['nux'])
            else:
                start_point = 0
                hydro_v = 0
                nu_flux = 0
                processed_hdf = []
        elif processed_hdf[-1].decode("utf-8") == simulation.hdf_file_list[-1] \
            or no_new:
            return integrate_momenta(time, PNSmass, hydro_v, nu_flux['nue'],
                                         nu_flux['nua'], nu_flux['nux'])
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
    PNS_radius, pgcells = simulation.PNS_radius(rad='full')
    PNS_radius = PNS_radius.data
    ## Get the grid
    dV = simulation.cell.dVolume_integration(simulation.ghost)
    dOmega = simulation.cell.dOmega(simulation.ghost)
    r400 = np.argmax(simulation.cell.radius(simulation.ghost) >= \
                     (400 * u.km))
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
            time = np.concatenate((time, simulation.time(file)))
            hydro_v = {'x': np.concatenate((hydro_v['x'], hydro_v_f[0])),
                       'y': np.concatenate((hydro_v['y'], hydro_v_f[1])),
                       'z': np.concatenate((hydro_v['z'], hydro_v_f[2]))}
            nu_flux['nue'] = {'x': np.concatenate((nu_flux['nue']['x'], nue_flux_f[0])),
                              'y': np.concatenate((nu_flux['nue']['y'], nue_flux_f[1])),
                              'z': np.concatenate((nu_flux['nue']['z'], nue_flux_f[2]))}
            nu_flux['nua'] = {'x': np.concatenate((nu_flux['nua']['x'], nua_flux_f[0])),
                              'y': np.concatenate((nu_flux['nua']['y'], nua_flux_f[1])),
                              'z': np.concatenate((nu_flux['nua']['z'], nua_flux_f[2]))}
            nu_flux['nux'] = {'x': np.concatenate((nu_flux['nux']['x'], nux_flux_f[0])),
                              'y': np.concatenate((nu_flux['nux']['y'], nux_flux_f[1])),
                              'z': np.concatenate((nu_flux['nux']['z'], nux_flux_f[2]))}
        except Exception as e:
            print('Error: ', e)
            time = simulation.time(file)
            hydro_v = {'x': hydro_v_f[0],
                       'y': hydro_v_f[1],
                       'z': hydro_v_f[2]}
            nu_flux = {'nue': {'x': nue_flux_f[0],
                               'y': nue_flux_f[1],
                               'z': nue_flux_f[2]},
                       'nua': {'x': nua_flux_f[0],
                               'y': nua_flux_f[1],
                               'z': nua_flux_f[2]},
                       'nux': {'x': nux_flux_f[0],
                               'y': nux_flux_f[1],
                               'z': nux_flux_f[2]}
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
    return integrate_momenta(time, PNSmass, hydro_v, nu_flux['nue'],
                             nu_flux['nua'], nu_flux['nux'])

def read_kick(simulation):
    """
    Reads the kick file.
    """
    with h5py.File(os.path.join(simulation.storage_path, 
                                  'kick_velocity.h5'), 'r') as f:
        time = aerray(f['time'][...], u.s, 'time',
               r'$t-t_\mathrm{b}$', None, [-0.005, f['time'][-1]])
        hydro_v = {'x': aerray(f['hydro/x'][...], u.cm/u.s*u.M_sun),
                   'y': aerray(f['hydro/y'][...], u.cm/u.s*u.M_sun),
                   'z': aerray(f['hydro/z'][...], u.cm/u.s*u.M_sun)}
        nu_flux = {'nue': {'x': aerray(f['nu_flux/nue/x'][...],u.M_sun*u.cm/u.s**2),
                           'y': aerray(f['nu_flux/nue/y'][...],u.M_sun*u.cm/u.s**2),
                           'z': aerray(f['nu_flux/nue/z'][...],u.M_sun*u.cm/u.s**2)},
                   'nua': {'x': aerray(f['nu_flux/nua/x'][...],u.M_sun*u.cm/u.s**2),
                           'y': aerray(f['nu_flux/nua/y'][...],u.M_sun*u.cm/u.s**2),
                           'z': aerray(f['nu_flux/nua/z'][...],u.M_sun*u.cm/u.s**2)},
                   'nux': {'x': aerray(f['nu_flux/nux/x'][...],u.M_sun*u.cm/u.s**2),
                           'y': aerray(f['nu_flux/nux/y'][...],u.M_sun*u.cm/u.s**2),
                           'z': aerray(f['nu_flux/nux/z'][...],u.M_sun*u.cm/u.s**2)}}
        if 'processed' in f:
            processed = f['processed'][...]
        else:
            processed = None
    return time, hydro_v, nu_flux, processed
    