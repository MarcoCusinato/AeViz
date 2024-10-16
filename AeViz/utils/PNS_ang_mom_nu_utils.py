import numpy as np
import h5py, os
from AeViz.utils.utils import (check_existence, progressBar, checkpoints)
from AeViz.utils.file_utils import save_hdf
from AeViz.units.units import units
from AeViz.grid.grid import grid

u = units()

def integrate_momenta(time, Lx, Ly, Lz):
    dt = np.zeros(time.shape[0])
    dt[1:] = time[1:] - time[:-1]
    dt[0] = dt[1]
    Lx = {key: np.cumsum(Lx[key] * dt) for key in Lx.keys()}
    Ly = {key: np.cumsum(Ly[key] * dt) for key in Ly.keys()}
    Lz = {key: np.cumsum(Lz[key] * dt) for key in Lz.keys()}
    return Lx, Ly, Lz

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
    if simulation.dim == 2:
        Flux = simulation.neutrino_momenta_grey(file_name) * \
               dOmega[..., None, None]
    else:
        Flux = simulation.neutrino_momenta_grey(file_name) * \
               dOmega[..., None, None, None]
    Fx, Fy, Fz = grid.spherical_to_cartesian(Flux[..., 0], Flux[..., 1],
                                           Flux[..., 2], add_back=1)
    del Flux
    PNS_surface = simulation.ghost.remove_ghost_cells_radii(PNS_radius,
                                                            simulation.dim,
                                                            **gcells)
    Omega = simulation.omega(file_name) * dOmega[..., None] / dOmega.sum()
    surface_indices = np.argmax(radius >= PNS_surface[..., None], axis=-1)
    if simulation.dim == 2:
        Omega = np.nansum(Omega[indices[0], surface_indices])
        Lx = -np.sum(Fx[indices[0], surface_indices, :] * \
                     PNS_surface[..., None] ** 2, axis=0) * \
                        PNS_avg_radius ** 2 * Omega
        Ly = -np.sum(Fy[indices[0], surface_indices, :] * \
                     PNS_surface[..., None] ** 2, axis=0) * \
                        PNS_avg_radius ** 2 * Omega
        Lz = -np.sum(Fz[indices[0], surface_indices, :] * \
                     PNS_surface[..., None] ** 2, axis=0) * \
                        PNS_avg_radius ** 2 * Omega
    else:
        Omega = np.nansum(Omega[indices[0], indices[1], surface_indices])
        Lx = -np.sum(Fx[indices[0], indices[1], surface_indices,:] * \
                    PNS_surface[..., None] ** 2, axis=(0, 1)) * \
                        PNS_avg_radius ** 2 * Omega
        Ly = -np.sum(Fy[indices[0], indices[1], surface_indices,:] * \
                    PNS_surface[..., None] ** 2, axis=(0, 1)) * \
                        PNS_avg_radius ** 2 * Omega
        Lz = -np.sum(Fz[indices[0], indices[1], surface_indices, :] * \
                    PNS_surface[..., None] ** 2, axis=(0, 1)) * \
                        PNS_avg_radius ** 2 * Omega
    
    return Lx, Ly, Lz

def calculate_angular_mom_PNS_nu(simulation, save_checkpoints=True):
    if simulation.dim == 1:
        print('No angular momentum calculation for 1D simulations.')
        return None
    if check_existence(simulation, 'PNS_angular_momentum_nu.h5'):
        time, Lx, Ly, Lz, processed_hdf = \
            read_angular_mom_PNS_nu(simulation)
        if processed_hdf[-1].decode("utf-8") == simulation.hdf_file_list[-1]:
            Lx, Ly, Lz = integrate_momenta(time, Lx, Ly, Lz)
            Ltotx = Lx['nue'] + Lx['nua'] + Lx['nux']
            Ltoty = Ly['nue'] + Ly['nua'] + Ly['nux']
            Ltotz = Lz['nue'] + Lz['nua'] + Lz['nux']
            Ltot = np.sqrt(Ltotx ** 2 + Ltoty ** 2 + Ltotz ** 2)
            return time, Lx, Ly, Lz, Ltotx, Ltoty, Ltotz, Ltot
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
    _, PNS_r, _, _, av_r, g_cells = simulation.PNS_radius()
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
        progressBar(progress_index, total_points, suffix='Computing...')
        Lx_comp, Ly_comp, Lz_comp = PNS_angular_momentum_neutrinos(
                                                    simulation, file,
                                                    PNS_r[..., findex],
                                                    av_r[findex],
                                                    indices, dOmega, gr,
                                                    r, g_cells)
        try:
            time = np.append(time, simulation.time(file))
            Lx = {'nue': np.append(Lx['nue'], Lx[..., 0]),
                  'nua': np.append(Lx['nua'], Lx[..., 1]),
                  'nux': np.append(Lx['nux'], Lx[..., 2])}
            Ly = {'nue': np.append(Ly['nue'], Ly[..., 0]),
                  'nua': np.append(Ly['nua'], Ly[..., 1]),
                  'nux': np.append(Ly['nux'], Ly[..., 2])}
            Lz = {'nue': np.append(Lz['nue'], Lz[..., 0]),
                  'nua': np.append(Lz['nua'], Lz[..., 1]),
                  'nux': np.append(Lz['nux'], Lz[..., 2])}
        except:
            time = np.array([simulation.time(file)])
            Lx = {'nue': np.array([Lx_comp[..., 0]]),
                  'nua': np.array([Lx_comp[..., 1]]),
                  'nux': np.array([Lx_comp[..., 2]])}
            Ly = {'nue': np.array([Ly_comp[..., 0]]),
                  'nua': np.array([Ly_comp[..., 1]]),
                  'nux': np.array([Ly_comp[..., 2]])}
            Lz = {'nue': np.array([Lz_comp[..., 0]]),
                  'nua': np.array([Lz_comp[..., 1]]),
                  'nux': np.array([Lz_comp[..., 2]])}
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
    Lx, Ly, Lz = integrate_momenta(time, Lx, Ly, Lz)
    Ltotx = Lx['nue'] + Lx['nua'] + Lx['nux']
    Ltoty = Ly['nue'] + Ly['nua'] + Ly['nux']
    Ltotz = Lz['nue'] + Lz['nua'] + Lz['nux']
    Ltot = np.sqrt(Ltotx ** 2 + Ltoty ** 2 + Ltotz ** 2)
    return time, Lx, Ly, Lz, Ltotx, Ltoty, Ltotz, Ltot

def read_angular_mom_PNS_nu(simulation):
    """
    Reads the angular momentum file.
    """
    with h5py.File(os.path.join(simulation.storage_path, 
                                  'PNS_angular_momentum_nu.h5'), 'r') as f:
        time = f['time'][...]
        Lx = {'nue': f['Lx/nue'][...],
              'nua': f['Lx/nua'][...],
              'nux': f['Lx/nux'][...]}
        Ly = {'nue': f['Ly/nue'][...],
              'nua': f['Ly/nua'][...],
              'nux': f['Ly/nux'][...]}
        Lz = {'nue': f['Lz/nue'][...],
              'nua': f['Lz/nua'][...],
              'nux': f['Lz/nux'][...]}
        processed = f['processed'][...]
        
        
    return time, Lx, Ly, Lz, processed

    