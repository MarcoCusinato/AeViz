from AeViz.utils.utils import (check_existence, progressBar, checkpoints)
from AeViz.simulation.simulation import Simulation
from AeViz.utils.file_utils import save_hdf
from AeViz.grid.grid import grid
import os, h5py
import numpy as np

def PNS_mass_energy(simulation, file_name, PNS_radius, gcells, dV, grid, gr):
    """
    Calculates the mass, kinetic, magnetic, rotational, gravitational
    and convective energy of the PNS for one timestep.
    """
    if simulation.dim == 1:
        mask = (simulation.cell.radius(simulation.ghost) <= \
            PNS_radius)
    else:
        mask = (simulation.cell.radius(simulation.ghost) <= \
            simulation.ghost.remove_ghost_cells_radii(PNS_radius,
                                                      simulation.dim,
                                                  **gcells)[..., None])
    rho = simulation.rho(file_name)[mask] * dV[mask]
    if simulation.dim == 1:
        Lx, Ly, Lz = 0, 0, 0
        L_tot = 0
    else:
        vr = simulation.radial_velocity(file_name)
        vt = simulation.theta_velocity(file_name)
        vp = simulation.phi_velocity(file_name)
        vx, vy, vz = gr.velocity_sph_to_cart(vr, vt, vp)
        masked_grid = [i[mask] if not type(i) == int else i for i in grid] 
        Lx, Ly, Lz = (rho * (vz[mask] * masked_grid[1] - vy[mask] * \
            masked_grid[2])).sum(), \
            (rho * (vx[mask] * masked_grid[2] - vz[mask] * \
                masked_grid[0])).sum(), \
            (rho * (vy[mask] * masked_grid[0] - vx[mask] * \
                masked_grid[1])).sum()
        L_tot = np.sqrt(Lx ** 2 + Ly ** 2 + Lz ** 2)
    return  Lx, Ly, Lz, L_tot

def calculate_masses_energies(simulation, save_checkpoints=True):
    if check_existence(simulation, 'MOM.h5'):
        PNS_me  = \
            read_masses_energies(simulation)
        if len(simulation.hdf_file_list) == len(PNS_me[0]):
            return PNS_me
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
    _, PNS_radius, _, _, _, pgcells = simulation.PNS_radius()
    ## Get the grid
    if simulation.dim == 1:
        gr = grid(1, simulation.cell.radius(simulation.ghost))
        X, Y, Z = (gr.cartesian_grid(), 0, 0)
    elif simulation.dim == 2:
        gr = grid(2, simulation.cell.radius(simulation.ghost), 
                  simulation.cell.theta(simulation.ghost))
        X, Z = gr.cartesian_grid()
        Y = 0
    else:
        gr = grid(3, simulation.cell.radius(simulation.ghost), 
                  simulation.cell.theta(simulation.ghost),
                  simulation.cell.phi(simulation.ghost))
        X, Y, Z = gr.cartesian_grid()

    ## Get the volume elements
    dV = simulation.cell.dVolume_integration(simulation.ghost)
    ## Indices
    findex = start_point
    check_index = 0
    progress_index = 0
    total_points = len(simulation.hdf_file_list) - start_point
    for file in simulation.hdf_file_list[start_point:]:
        progressBar(progress_index, total_points, suffix='Computing...')
        PNS_data = PNS_mass_energy(simulation, file, PNS_radius[..., findex],
                                   pgcells, dV, (X, Y, Z), gr)
        try:
            time = np.concatenate((time, simulation.time(file)))
            
            ## PNS
            PNS_me['L']['Lx'] = np.concatenate((PNS_me['L']['Lx'],
                                                np.array([PNS_data[0]])))
            PNS_me['L']['Ly'] = np.concatenate((PNS_me['L']['Ly'],
                                                np.array([PNS_data[1]])))
            PNS_me['L']['Lz'] = np.concatenate((PNS_me['L']['Lz'],
                                                np.array([PNS_data[2]])))
            PNS_me['L']['L_tot'] = np.concatenate((PNS_me['L']['L_tot'],
                                                np.array([PNS_data[3]])))
            
        except:
            time = simulation.time(file)           
            PNS_me = {
                'L': {
                    'Lx': np.array([PNS_data[0]]),
                    'Ly': np.array([PNS_data[1]]),
                    'Lz': np.array([PNS_data[2]]),
                    'L_tot': np.array([PNS_data[3]])
                }
            }
        if (check_index >= checkpoint) and save_checkpoints:
            print('Checkpoint reached, saving...\n')
            save_hdf(os.path.join(simulation.storage_path,
                    'MOM.h5'),
                     ['PNS'],
                     [PNS_me])
            
            check_index = 0
        check_index += 1
        progress_index += 1
        findex += 1
    print('Computation completed, saving...')
    save_hdf(os.path.join(simulation.storage_path, 'MOM.h5'),
                    ['PNS'],
                    [ PNS_me])
    return None
        

def read_masses_energies(simulation):
    masses_energies_data = h5py.File(os.path.join(simulation.storage_path, 
                                            'MOM.h5'), 'r')
    data = [
        {
            'L': {
                'Lx': masses_energies_data['PNS/L/Lx'][...],
                'Ly': masses_energies_data['PNS/L/Ly'][...],
                'Lz': masses_energies_data['PNS/L/Lz'][...],
                'L_tot': masses_energies_data['PNS/L/L_tot'][...]
            }
        }
    ]
    return data


sims = ['A05-1', 'A08-1', 'A17-1', 'A20-1', 'A26-1', 'A30-1', 'A39-1']
for sim in sims:
    print('Computing masses and energies for simulation: {}'.format(sim))
    simulation = Simulation(sim)
    calculate_masses_energies(simulation)
    print('Masses and energies for simulation {} computed and saved.'.format(sim))
    del simulation