
from __future__ import annotations
from AeViz.utils.utils import (check_existence, progressBar, checkpoints)
from AeViz.simulation.simulation import Simulation
import os
from AeViz.utils.physics.PNS_postprocessing import (declare_PNS_dictionary,
                                                    compute_PNS_postprocessing,
                                                    PNS_quantity_metadata)
from AeViz.units import u
from AeViz.utils.files.file_utils import load_dataset, save_merge_dictionary_hdf
import numpy as np

def supernova_postprocessing(simulation: Simulation,
                             save_checkpoints: bool = True) -> None:
    if (checkpoints[simulation.dim] == False) or (not save_checkpoints):
        checkpoint = len(simulation.hdf_file_list)
    else:
        checkpoint = checkpoints[simulation.dim]
    
    dV = simulation.cell.dVolume_integration(simulation.ghost)
    radius = simulation.cell.radius(simulation.ghost)
    grid = simulation.grid.cartesian_grid()
    dOmega = simulation.cell.dOmega(simulation.ghost)
    file_list = simulation.hdf_file_list
    pns_processed = []
    ## unpack the grid
    if simulation.dim == 1:
        X = grid * radius.unit
        rperp2 = radius ** 2
    elif simulation.dim == 2:
        rperp2 = (radius[None, :] *
                  np.sin(simulation.cell.theta(simulation.ghost))[:, None]) ** 2
        X, Y, Z = grid[0] * radius.unit, grid[1] * radius.unit, 0 * radius.unit
    elif simulation.dim == 3:
        rperp2 = (radius[None, None, :] * 
                  np.sin(simulation.cell.theta(simulation.ghost))[None, :, None]) ** 2
        X, Y, Z = grid[0] * radius.unit, X[1] * radius.unit, grid[3] * radius.unit
    if check_existence(simulation, 'PNS_postprocessing.h5'):
        pns_processed = load_dataset(simulation.storage_path, 
                                     'PNS_postprocessing.h5',
                                     'processed_hdf', False)
    start_index = len(pns_processed)
    file_list = file_list[start_index:]
    tot_points = len(file_list)
    PNS_dict = declare_PNS_dictionary(simulation.dim, simulation.magdim)
    pns_metadata = PNS_quantity_metadata()
    pns_time = []
    check_index = 0
    for findex, file_name in enumerate(file_list):
        progressBar(findex, tot_points, 'Computing PNS postprocessing...')
        ## Compute the angular momentum
        dmass = simulation.rho(file_name) * dV
        time = simulation.time(file_name)
        if simulation.dim == 1:
            jx, jy, jz = None, None, None
            jtot = None
        else:
            vr = simulation.radial_velocity(file_name)
            vt = simulation.theta_velocity(file_name)
            vp = simulation.phi_velocity(file_name)
            vx, vy, vz = simulation.grid.spherical_to_cartesian(dmass * vr,
                                                                dmass * vt,
                                                                dmass * vp)
            jx, jy, jz = (vz * X - vy * Z), (vx * Z - vz * X), (vy * X - vx * Y)
            jtot = np.sqrt(jx ** 2 + jy ** 2 + jz ** 2)
            shell_jx = np.nansum(jx, axis=range(simulation.dim - 1))
            shell_jy = np.nansum(jy, axis=range(simulation.dim - 1))
            shell_jz = np.nansum(jz, axis=range(simulation.dim - 1))
            shell_jtot = np.nansum(jtot, axis=range(simulation.dim - 1))
        ## Compute inertia moment
        if simulation.dim < 3:
            Inertia = dmass * rperp2
        else:
            nx = shell_jx / shell_jtot
            ny = shell_jy / shell_jtot
            nz = shell_jz / shell_jtot
            mm = np.isnan(nx) | np.isnan(ny) | np.isnan(nz)
            nx[mm] = 0.
            ny[mm] = 0.
            nz[mm] = 1.
            n = np.array([nx, ny, nz])
            X_perp = np.linalg.norm((grid - np.einsum('ir,iptr->ptr', n, grid)[None, ...]
                                      * n[:, None, None, :]), axis=0) + 1e-12
            X_perp2 = np.where((n[2] >= .99)[None, None, :], rperp2, X_perp ** 2)
            Inertia = dmass * X_perp2

        compute_PNS_postprocessing(PNS_dict, simulation, file_name,
                                   dV, dOmega,
                                   jx, jy, jz)
        pns_time.append(time)
        pns_processed.append(file_name)
        if (check_index >= checkpoint) and save_checkpoints:
            print('Checkpoint reached, saving files...\n')
            save_merge_dictionary_hdf(simulation, (pns_time, PNS_dict, pns_processed),
                                      simulation.storage_path,
                                      'PNS_postprocessing.h5', **pns_metadata)
            check_index = 0
            PNS_dict = declare_PNS_dictionary(simulation.dim, simulation.magdim)
            pns_time = []
        else:
            check_index += 1