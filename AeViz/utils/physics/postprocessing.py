
from __future__ import annotations
from AeViz.utils.utils import (check_existence, progressBar, checkpoints)
from AeViz.simulation.simulation import Simulation
import os
from AeViz.utils.physics.PNS_postprocessing import (declare_PNS_dictionary,
                                                    compute_PNS_postprocessing)
from AeViz.units import u
from AeViz.utils.files.file_utils import load_dataset
import numpy as np

def supernova_postprocessing(simulation: Simulation,
                             save_checkpoints: bool = True) -> None:
    if (checkpoints[simulation.dim] == False) or (not save_checkpoints):
        checkpoint = len(simulation.hdf_file_list)
    else:
        checkpoint = checkpoints[simulation.dim]
    
    dV = simulation.cell.dVolume_integration(simulation.ghost)
    radius = simulation.cell.radius(simulation.ghost)
    X = simulation.grid.cartesian_grid()
    dOmega = simulation.cell.dOmega(simulation.ghost)
    file_list = simulation.hdf_file_list
    pns_processed = []
    ## unpack the grid
    if simulation.dim == 1:
        X = X * radius.unit
        rperp2 = radius ** 2
    elif simulation.dim == 2:
        rperp2 = (radius[None, :] *
                  np.sin(simulation.cell.theta)[:, None]) ** 2
        X, Y, Z = X[0] * radius.unit, X[1] * radius.unit, 0 * radius.unit
    elif simulation.dim == 3:
        rperp2 = (radius[None, None, :] * 
                  np.sin(simulation.cell.theta)[None, :, None]) ** 2
        X, Y, Z = X[0] * radius.unit, X[1] * radius.unit, X[3] * radius.unit
    if check_existence(simulation, 'PNS_postprocessing.h5'):
        pns_processed = load_dataset(simulation.storage_path, 
                                     'PNS_postprocessing.h5',
                                     'processed_hdf', False)
    start_index = len(pns_processed)
    file_list = file_list[start_index:]
    tot_points = len(file_list)
    for findex, file_name in enumerate(file_list):
        progressBar(findex, tot_points, 'Computing PNS postprocessing...')
        ## Compute the angular momentum
        rho = simulation.rho(file_name) * dV
        if simulation.dim == 1:
            jx, jy, jz = None, None, None
            jtot = None
        else:
            vr = simulation.radial_velocity(file_name)
            vt = simulation.theta_velocity(file_name)
            vp = simulation.phi_velocity(file_name)
            vx, vy, vz = simulation.grid.spherical_to_cartesian(rho * vr,
                                                                rho * vt,
                                                                rho * vp)
            jx, jy, jz = (vz * X - vy * Z), (vx * Z - vz * X), (vy * X - vx * Y)
            jtot = np.sqrt(jx ** 2 + jy ** 2 + jz ** 2)
            shell_jx = np.nansum(jx, axis=range(simulation.dim - 1))
            shell_jy = np.nansum(jy, axis=range(simulation.dim - 1))
            shell_jz = np.nansum(jz, axis=range(simulation.dim - 1))
            shell_jtot = np.nansum(jtot, axis=range(simulation.dim - 1))
        ## Compute inertia moment
        if simulation.dim < 3:
            Inertia = rho * rperp2
        else:
            nx = jx / jtot
            ny = jy / jtot
            nz = jz / jtot
            
        
    
    
    