import numpy as np
import h5py, os
from AeViz.utils.utils import (check_existence, progressBar, checkpoints)
from AeViz.utils.files.file_utils import save_hdf, create_series
from AeViz.units import u
from AeViz.grid.grid import grid
from AeViz.units. aerray import aerray
from AeViz.units.constants import constants as c

def integrate_momenta(time, Lx, Ly, Lz):
    dt = np.zeros(time.shape[0])
    dt[1:] = time[1:] - time[:-1]
    dt[0] = dt[1]
    dt = dt * time.unit
    Lx = {key: np.cumsum(Lx[key] * dt) / c.c ** 2 
          for key in Lx.keys()}
    Ly = {key: np.cumsum(Ly[key] * dt) / c.c ** 2
          for key in Ly.keys()}
    Lz = {key: np.cumsum(Lz[key] * dt) / c.c ** 2
          for key in Lz.keys()}
    return Lx, Ly, Lz

def create_aerrays(time, Lx, Ly, Lz):
    nus = [r'$\nu_e$', r'$\overline{\nu}_e$', r'$\nu_x$']
    Lx, Ly, Lz = integrate_momenta(time, Lx, Ly, Lz)
    [Lx[key].set(name='PNS_Lx_'+key, limits=[-1e47, 1e47],
                 label=r'$L_\mathrm{'+nu+r',x}$', log=True) for (key, nu)
                 in zip(Lx.keys(), nus)]
    [Ly[key].set(name='PNS_Ly_'+key, limits=[-1e47, 1e47],
                 label=r'$L_\mathrm{'+nu+r',y}$', log=True) for (key, nu)
                 in zip(Ly.keys(), nus)]
    [Lz[key].set(name='PNS_Lz_'+key, limits=[0, 1e49],
                 label=r'$L_\mathrm{'+nu+r',z}$', log=True) for (key, nu)
                 in zip(Lz.keys(), nus)]
    Ltotx = Lx['nue'] + Lx['nua'] + Lx['nux']
    Ltotx.set(name='PNS_Ltotx_nu', limits=[-1e47, 1e47],
               label=r'$L_\mathrm{\nu,x}$', log=True)
    Ltoty = Ly['nue'] + Ly['nua'] + Ly['nux']
    Ltoty.set(name='PNS_Ltoty_nu', limits=[-1e47, 1e47],
                label=r'$L_\mathrm{\nu,y}$', log=True)
    Ltotz = Lz['nue'] + Lz['nua'] + Lz['nux']
    Ltotz.set(name='PNS_Ltotz_nu', limits=[0, 1e49],
                label=r'$L_\mathrm{\nu,z}$', log=True)
    Ltot = np.sqrt(Ltotx ** 2 + Ltoty ** 2 + Ltotz ** 2)
    Ltot.set(name='PNS_Ltot_nu', limits=[0, 1e49],
             label=r'$L_\mathrm{\nu}$', log=True)
    Lx = {key: Lx[key].to(u.erg*u.s) for key in Lx.keys()}
    Ly = {key: Ly[key].to(u.erg*u.s) for key in Ly.keys()}
    Lz = {key: Lz[key].to(u.erg*u.s) for key in Lz.keys()}
    return create_series(time, Lx, Ly, Lz, Ltotx, Ltoty, Ltotz, Ltot)

def PNS_angular_momentum_neutrinos(simulation, file_name, PNS_radius,
                                   PNS_avg_radius, indices, dOmega, grid,
                                   radius, gcells):
    """
    Function to calculate the angular momentum of the PNS due to neutrinos.
    input:
        simulation: Simulation object
        file_name: string
        PNS_radius: array of floats, 1 or 2 dimensional
        PNS_avg_radius: float
        indices: tupleof array of integers or meshgrid
        dOmega: array of floats
        grid: grid object
        gcells: dictionary
    output:
        Lx, Ly, Lz: array of floats        
    """
    
    Flux = simulation.neutrino_momenta_grey(file_name)
    Flux = [fl * dOmega[..., None, None] for fl in Flux] #nue, nua, nux
    Fx = {}
    Fy = {}
    Fz = {}
    Fx['nue'], Fy['nue'], Fz['nue'] = grid.spherical_to_cartesian(Flux[0][..., 0],
                                                                  Flux[0][..., 1],
                                                                  Flux[0][..., 2])
    Fx['nua'], Fy['nua'], Fz['nua'] = grid.spherical_to_cartesian(Flux[1][..., 0],
                                                                  Flux[1][..., 1],
                                                                  Flux[1][..., 2])
    Fx['nux'], Fy['nux'], Fz['nux'] = grid.spherical_to_cartesian(Flux[2][..., 0],
                                                                  Flux[2][..., 1],
                                                                  Flux[2][..., 2])
    PNS_surface = simulation.ghost.remove_ghost_cells_radii(PNS_radius,
                                                            simulation.dim,
                                                            **gcells)
    Omega = simulation.omega(file_name) * dOmega[..., None] / dOmega.sum()
    surface_indices = np.argmax(radius >= PNS_surface[..., None], axis=-1)
    if simulation.dim == 2:
        Omega = np.nansum(Omega[indices[0], surface_indices])
        Lx = {key: -np.sum(Fx[key][indices[0], surface_indices] * 
                           PNS_surface ** 2, axis=0) * PNS_avg_radius ** 2 *
                           Omega for key in Fx.keys()}
        Ly = {key: -np.sum(Fy[key][indices[0], surface_indices] * 
                            PNS_surface ** 2, axis=0) * PNS_avg_radius ** 2 *
                            Omega for key in Fy.keys()}
        Lz = {key: -np.sum(Fz[key][indices[0], surface_indices] * 
                            PNS_surface ** 2, axis=0) * PNS_avg_radius ** 2 * 
                            Omega for key in Fz.keys()}
    else:
        Omega = np.nansum(Omega[indices[0], indices[1], surface_indices])
        Lx = {key: -np.sum(Fx[indices[0], indices[1], surface_indices] * 
                            PNS_surface ** 2, axis=(0, 1)) * 
                            PNS_avg_radius ** 2 * Omega for key in Fx.keys()}
        Ly = {key: -np.sum(Fy[indices[0], indices[1], surface_indices] * 
                            PNS_surface[..., None] ** 2, axis=(0, 1)) * 
                            PNS_avg_radius ** 2 * Omega for key in Fy.keys()}
        Lz = {key: -np.sum(Fz[indices[0], indices[1], surface_indices] * 
                        PNS_surface[..., None] ** 2, axis=(0, 1)) * 
                        PNS_avg_radius ** 2 * Omega for key in Fx.keys()}
    
    return Lx, Ly, Lz

def calculate_angular_mom_PNS_nu(simulation, save_checkpoints=True):
    if simulation.dim == 1:
        print('No angular momentum calculation for 1D simulations.')
        return None
    if check_existence(simulation, 'PNS_angular_momentum_nu.h5'):
        time, Lx, Ly, Lz, processed_hdf = \
            read_angular_mom_PNS_nu(simulation)
        if processed_hdf[-1].decode("utf-8") == simulation.hdf_file_list[-1]:
            return create_aerrays(time, Lx, Ly, Lz)
        else:
            start_point = len(processed_hdf)
            processed_hdf = [ff.decode("utf-8") for ff in processed_hdf]
            print('Checkpoint found for the angular momentum file, starting' \
                  ' from checkpoint.\nPlease wait...')
    else:
        start_point = 0
        processed_hdf = []
        print('No checkpoint found for the angular momentum file, starting' \
              ' from the beginning.\nPlease wait...')
    if (checkpoints[simulation.dim] == False) or (not save_checkpoints):
        checkpoint = len(simulation.hdf_file_list)
    else:
        checkpoint = checkpoints[simulation.dim]
    PNS_r, g_cells = simulation.PNS_radius(rad='full')
    PNS_r = PNS_r.data
    av_r = simulation.PNS_radius(rad='avg').data
    r = simulation.cell.radius(simulation.ghost)
    dOmega = simulation.cell.dOmega(simulation.ghost)
    gr = grid(simulation.dim, r, simulation.cell.theta(simulation.ghost),
              simulation.cell.phi(simulation.ghost))
    while r.ndim < simulation.dim:
        r = r[None, :]
    
    if simulation.dim == 2:
        indices = (np.arange(len(simulation.cell.theta(simulation.ghost))), )
    else:
        indices = np.meshgrid(
            np.arange(len(simulation.cell.phi(simulation.ghost))),
            np.arange(len(simulation.cell.theta(simulation.ghost))),
            indexing='ij')
    findex = start_point
    check_index = 0
    progress_index = 0
    total_points = len(simulation.hdf_file_list) - start_point
    for file in simulation.hdf_file_list[start_point:]:
        progressBar(progress_index, total_points, suffix='Computing neutrino '
                    'angular momentum component...')
        Lx_comp, Ly_comp, Lz_comp = PNS_angular_momentum_neutrinos(
                                                    simulation, file,
                                                    PNS_r[..., findex],
                                                    av_r[findex],
                                                    indices, dOmega, gr,
                                                    r, g_cells)
        try:
            time = np.concatenate((time, simulation.time(file)))
            Lx = {'nue': np.concatenate((Lx['nue'], Lx_comp['nue'])),
                  'nua': np.concatenate((Lx['nua'], Lx_comp['nua'])),
                  'nux': np.concatenate((Lx['nux'], Lx_comp['nux']))}
            Ly = {'nue': np.concatenate((Ly['nue'], Ly_comp['nue'])),
                  'nua': np.concatenate((Ly['nua'], Ly_comp['nua'])),
                  'nux': np.concatenate((Ly['nux'], Ly_comp['nux']))}
            Lz = {'nue': np.concatenate((Lz['nue'], Lz_comp['nue'])),
                  'nua': np.concatenate((Lz['nua'], Lz_comp['nua'])),
                  'nux': np.concatenate((Lz['nux'], Lz_comp['nux']))}
        except Exception as e:
            time = simulation.time(file)
            Lx = {'nue': Lx_comp['nue'],
                  'nua': Lx_comp['nua'],
                  'nux': Lx_comp['nux']}
            Ly = {'nue': Ly_comp['nue'],
                  'nua': Ly_comp['nua'],
                  'nux': Ly_comp['nux']}
            Lz = {'nue': Lz_comp['nue'],
                  'nua': Lz_comp['nua'],
                  'nux': Lz_comp['nux']}
        processed_hdf.append(file)
        if (check_index >= checkpoint) and save_checkpoints:
            print('Checkpoint reached, saving...\n')
            save_hdf(os.path.join(simulation.storage_path, 
                                  'PNS_angular_momentum_nu.h5'),
                     ['time', 'Lx', 'Ly', 'Lz', 'processed'],
                     [time, Lx, Ly, Lz, processed_hdf])
            check_index = 0
        check_index += 1
        progress_index += 1
        findex += 1
    save_hdf(os.path.join(simulation.storage_path, 
                        'PNS_angular_momentum_nu.h5'),
            ['time', 'Lx', 'Ly', 'Lz', 'processed'],
            [time, Lx, Ly, Lz, processed_hdf])
    return create_aerrays(time, Lx, Ly, Lz)

def read_angular_mom_PNS_nu(simulation):
    """
    Reads the angular momentum file.
    """
    with h5py.File(os.path.join(simulation.storage_path, 
                                  'PNS_angular_momentum_nu.h5'), 'r') as f:
        time = aerray(f['time'][...], u.s, 'time',
               r'$t-t_\mathrm{b}$', None, [-0.005, f['time'][-1]])
        Lx = {'nue': aerray(f['Lx/nue'][...], u.erg*u.cm**2/u.s**2),
              'nua': aerray(f['Lx/nua'][...], u.erg*u.cm**2/u.s**2),
              'nux': aerray(f['Lx/nux'][...], u.erg*u.cm**2/u.s**2)}
        Ly = {'nue': aerray(f['Ly/nue'][...], u.erg*u.cm**2/u.s**2),
              'nua': aerray(f['Ly/nua'][...], u.erg*u.cm**2/u.s**2),
              'nux': aerray(f['Ly/nux'][...], u.erg*u.cm**2/u.s**2)}
        Lz = {'nue': aerray(f['Lz/nue'][...], u.erg*u.cm**2/u.s**2),
              'nua': aerray(f['Lz/nua'][...], u.erg*u.cm**2/u.s**2),
              'nux': aerray(f['Lz/nux'][...], u.erg*u.cm**2/u.s**2)}
        processed = f['processed'][...]
           
    return time, Lx, Ly, Lz, processed

    