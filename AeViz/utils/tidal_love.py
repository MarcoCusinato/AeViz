
from AeViz.units import u
from AeViz.utils.file_utils import save_hdf
from AeViz.utils.math_utils import function_average
from AeViz.utils.utils import check_existence, progressBar, checkpoints
from scipy.interpolate import Akima1DInterpolator
from scipy.integrate import solve_ivp
import h5py
import numpy as np
import os

"""
Everything here is taken from the following paper:
https://doi.org/10.1103/PhysRevD.81.123016
In particular, equation (12) for the linearized metric, equation (14)
for the Love number.
"""

def linearized_metric(r, y, p, m, cs, rho):
    dH_dr, H = y
    if r == 0:
        return [y[0], 0]
    
    coeff = 1 / (1 - 2 * m(r) / r)

    d2H_dr2 = 2 * coeff * H * \
        (-2 * np.pi * (5 * rho(r) + 9 * p(r) + (rho(r) + p(r)) / cs(r) ** 2) + \
         3 / r ** 2 + 2 * coeff * (m(r) / r ** 2 + 4 * np.pi * r * p(r)) ** 4) + \
         2 * dH_dr / r * coeff * (-1 + m(r) / r + 2 * np.pi * r ** 2 * (rho(r) - p(r)))
    return [d2H_dr2, dH_dr]

def Love_number(xi, y):
    """
    Compute the love number for a given compacness and
    R H'(R)/H(R) = y
    """
    XI = (1 - 2 * xi)
    kappa2 = 8 / 5 * xi ** 5 * XI ** 2 * (2 + 2 * xi * (y- 1) - y) / \
          (2 * xi * (6 - 3 * y + 3 * xi * (5 * y - 8)) + 4 * xi ** 3 * \
           (13 - 11 * y + xi * (3 * y - 2) + 2 * xi ** 2 * (1 + y)) + \
            3 * XI ** 2 * (2 - y+ 2 * xi * (y- 1)) * (np.log(XI)))
    return kappa2

def tidal_deformability(kappa2, xi):
    """
    Compute the tidal deformability from the love number
    """
    return 2 / 3 * kappa2 / xi ** 5

def solve_tidal_love(xi, pres, mass, soundspeed, dens, radius):
    """
    Solve the tidal Love number for a given compactness xi
    """
    pres = Akima1DInterpolator(radius, pres, method='akima', extrapolate=True)
    dens = Akima1DInterpolator(radius, dens, method='akima', extrapolate=True)
    mass = Akima1DInterpolator(radius, mass, method='akima', extrapolate=True)
    soundspeed = Akima1DInterpolator(radius, soundspeed, method='akima',
                                     extrapolate=True)
    solution = solve_ivp(linearized_metric, [radius[0], radius[-1]],
                         [2 * radius[0], radius[0] ** 2],
                         method='LSODA', t_eval=radius,
                         args=[pres, mass, soundspeed, dens], max_step=10000)
    y = solution.t[-1] * solution.y[0][-1] / solution.y[1][-1]
    kappa2 = Love_number(xi, y)
    tidal_d = tidal_deformability(kappa2, xi)
    return kappa2, tidal_d

def solve_tidal_love_profile(simulation, save_checkpoints=True):
    """
    Derives and saves the tidal love number and tidal deformability for
    the PNS and PNS core.
    """
    if check_existence(simulation, 'tidal.h5'):
        time, pns, core, processed_hdf = \
            read_tidal(simulation)
        if processed_hdf[-1].decode("utf-8") == simulation.hdf_file_list[-1]:
            return time, pns, core
        else:
            start_point = len(processed_hdf)
            processed_hdf = [ff.decode("utf-8") for ff in processed_hdf]
            print('Checkpoint found for the tidal file, starting' \
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
    t, _, _, _, PNS_radius, _ = simulation.PNS_radius()
    _, _, _, _, core_radius, _ = simulation.PNS_radius()
    ## Get the calculated profiles
    _, radius, rho_prof = simulation.radial_profile('rho')
    _, _, pgas_prof = simulation.radial_profile('gas_pressure')
    ## Convert to cactus units
    radius = u.convert_to_cactus_lenght(radius)
    rho_prof = u.convert_to_cactus_rho(rho_prof)
    pgas_prof = u.convert_to_cactus_pressure(pgas_prof)
    PNS_radius = u.convert_to_cactus_lenght(PNS_radius)
    core_radius = u.convert_to_cactus_lenght(core_radius)
    dOmega = simulation.cell.dOmega(simulation.ghost)
    dV = simulation.cell.dVolume_integration(simulation.ghost)
    ## Indices
    findex = start_point
    check_index = 0
    progress_index = 0
    total_points = len(simulation.hdf_file_list) - start_point
    for file in simulation.hdf_file_list[start_point:]:
        progressBar(progress_index, total_points, suffix='Computing...')
        if t[findex] < 0:
            tidal_core, tidal_pns = 0, 0
            love_core, love_pns = 0, 0
        else:
            ##compute the speed of sound profile
            csound = u.convert_from_cactus_velocity(
                function_average(simulation.soundspeed(file), simulation.dim,
                                'Omega', dOmega))
            ##compute the mass profile
            mass = u.convert_to_solar_masses(np.cumsum(
                        np.sum(simulation.rho(file) * dV,
                            axis=tuple(range(simulation.dim - 1)))))
            ## compute the core and PNS radius indices
            core_index = np.argmax(radius > core_radius[findex])
            pns_index = np.argmax(radius > PNS_radius[findex])
            ## compute the tidal love number and tidal deformability
            love_core, tidal_core = solve_tidal_love(
                mass[core_index] / core_radius[findex],
                pgas_prof[:core_index, findex], mass[:core_index, findex],
                csound[:core_index, findex], rho_prof[:core_index, findex],
                radius[:core_index, findex])
            love_pns, tidal_pns = solve_tidal_love(
                mass[pns_index] / PNS_radius[findex],
                pgas_prof[:pns_index, findex], mass[:pns_index, findex],
                csound[:pns_index, findex], rho_prof[:pns_index, findex],
                radius[:pns_index, findex])
        ## Save the results into dictionaries
        try:
            time = np.concatenate((time, [t[findex]]))
            pns['kappa2'] = np.concatenate((pns['kappa2'], [love_pns]))
            pns['lambda'] = np.concatenate((pns['lambda'], [tidal_pns]))
            core['kappa2'] = np.concatenate((core['kappa2'], [love_core]))
            core['lambda'] = np.concatenate((core['lambda'], [tidal_core]))
        except:
            time = np.array([t[findex]])
            pns = {'kappa2': np.array([love_pns]),
                   'lambda': np.array([tidal_pns])}
            core = {'kappa2': np.array([love_core]),
                    'lambda': np.array([tidal_core])}
        processed_hdf.append(file)
        if (check_index >= checkpoint) and save_checkpoints:
            print('Checkpoint reached, saving...\n')
            save_hdf(os.path.join(simulation.storage_path, 'tidal.h5'),
                     ['time', 'PNS', 'PNS_core', 'processed'],
                     [time, pns, core, processed_hdf])
            
            check_index = 0
        check_index += 1
        progress_index += 1
        findex += 1
    print('Computation completed, saving...')
    save_hdf(os.path.join(simulation.storage_path, 'tidal.h5'),
                     ['time', 'PNS', 'PNS_core', 'processed'], 
                     [time, pns, core, processed_hdf])
    return time, pns, core

def read_tidal(simulation):
    """
    Read the tidal file
    """
    data_file = h5py.File(os.path.join(simulation.storage_path, 'tidal.h5'), 'r')
    data = [
        data_file['time'][...],
        {
            'kappa2': data_file['PNS']['kappa2'][...],
            'lambda': data_file['PNS']['lambda'][...]
        },
        {
            'kappa2': data_file['PNS_core']['kappa2'][...],
            'lambda': data_file['PNS_core']['lambda'][...]
        },
        data_file['processed'][...]
    ]
    return data