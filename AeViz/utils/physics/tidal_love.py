
from AeViz.units import u
from AeViz.units.aeseries import aerray, aeseries
from AeViz.units.constants import constants as c
from AeViz.utils.files.file_utils import save_hdf, create_series
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

def to_cactus_len(x):
    """
    Convert the radius to cactus units
    """
    return x / c.G / c.Msol * c.c ** 2

def to_cactus_time(x):
    """
    Convert the time to cactus units
    """
    return x / c.G / c.Msol * c.c ** 3

def to_cactus_dens(x):
    """
    Convert the density to cactus units
    """
    return x / c.Msol * c.G ** 3  * c.Msol ** 3 / c.c ** 6

def to_cactus_vel(x):
    """
    Convert the velocity to cactus units
    """
    return x / c.c

def to_cactus_pres(x):
    """
    Convert the pressure to cactus units
    """
    return to_cactus_dens(x) / c.c ** 2

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
        if processed_hdf[-1].decode("utf-8") == simulation.hdf_file_list[-1] or \
            simulation.no_new:
            time = aerray(time, u.s, name='time', label=r'$t-t_\mathrm{b}$',
                          limits=[-0.05, time[-1]])
            lambda_pns = aerray(pns['lambda'], u.dimensionless_unscaled,
                                name='lambda_pns',  limits=[0, 10000], log=True,
                                label=r'$\Lambda_\mathrm{PNS}$')
            kappa_pns = aerray(pns['kappa2'], u.dimensionless_unscaled,
                               name='kappa_pns', limits=[0, 0.002],
                               label=r'$\kappa_2^\mathrm{PNS}$')
            lambda_core = aerray(core['lambda'], u.dimensionless_unscaled,
                                name='lambda_core', limits=[0, 10000], log=True,
                                label=r'$\Lambda_\mathrm{core}$')
            kappa_core = aerray(core['kappa2'], u.dimensionless_unscaled,
                               name='kappa_core', limits=[0, 0.002],
                               label=r'$\kappa_2^\mathrm{core}$')
            return create_series(time, lambda_pns, kappa_pns, lambda_core,
                                  kappa_core)
        else:
            start_point = len(processed_hdf)
            processed_hdf = [ff.decode("utf-8") for ff in processed_hdf]
            print('Checkpoint found for the tidal file, starting' \
                  ' from checkpoint.\nPlease wait...')
    else:
        start_point = 0
        processed_hdf = []
        print('No checkpoint found for the tidal deformablity file, starting' \
              ' from the beginning.\nPlease wait...')
    if (checkpoints[simulation.dim] == False) or (not save_checkpoints):
        checkpoint = len(simulation.hdf_file_list)
    else:
        checkpoint = checkpoints[simulation.dim]
    ## Get the radii
    PNS_radius = simulation.PNS_radius(rad='avg')
    core_radius = simulation.PNS_radius(rad='avg')
    ## Get the calculated profiles
    rho_prof = simulation.radial_profile('rho')
    pgas_prof = simulation.radial_profile('gas_pressure')
    ## Convert to cactus units
    t = PNS_radius.time
    radius = to_cactus_len(rho_prof.radius)
    rho_prof = to_cactus_dens(rho_prof.data)
    pgas_prof = to_cactus_pres(pgas_prof.data)
    PNS_radius = to_cactus_len(PNS_radius.data)
    core_radius = to_cactus_len(core_radius.data)
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
            csound = to_cactus_vel(
                function_average(simulation.soundspeed(file), simulation.dim,
                                'Omega', dOmega))
            ##compute the mass profile
            mass = np.cumsum(np.sum(simulation.rho(file) * dV,
                            axis=tuple(range(simulation.dim - 1))).to(u.M_sun)) / u.M_sun
            ## compute the core and PNS radius indices
            core_index = np.argmax(radius > core_radius[findex])
            pns_index = np.argmax(radius > PNS_radius[findex])
            ## compute the tidal love number and tidal deformability
            love_core, tidal_core = solve_tidal_love(
                mass[core_index] / core_radius[findex],
                pgas_prof[:core_index, findex], mass[:core_index],
                csound[:core_index], rho_prof[:core_index, findex],
                radius[:core_index])
            love_pns, tidal_pns = solve_tidal_love(
                mass[pns_index] / PNS_radius[findex],
                pgas_prof[:pns_index, findex], mass[:pns_index],
                csound[:pns_index], rho_prof[:pns_index, findex],
                radius[:pns_index])
        ## Save the results into dictionaries
        try:
            time = np.concatenate((time, [t[findex]]))
            pns['kappa2'] = np.concatenate((pns['kappa2'], [love_pns]))
            pns['lambda'] = np.concatenate((pns['lambda'], [tidal_pns]))
            core['kappa2'] = np.concatenate((core['kappa2'], [love_core]))
            core['lambda'] = np.concatenate((core['lambda'], [tidal_core]))
        except Exception as e:
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
    time = aerray(time, u.s, name='time', label=r'$t-t_\mathrm{b}$',
                  limits=[-0.05, time[-1]])
    lambda_pns = aerray(pns['lambda'], u.dimensionless_unscaled,
                        name='lambda_pns',  limits=[0, 10000], log=True,
                        label=r'$\Lambda_\mathrm{PNS}$')
    kappa_pns = aerray(pns['kappa2'], u.dimensionless_unscaled,
                       name='kappa_pns', limits=[0, 0.002],
                        label=r'$\kappa_2^\mathrm{PNS}$')
    lambda_core = aerray(core['lambda'], u.dimensionless_unscaled,
                        name='lambda_core', limits=[0, 10000], log=True,
                        label=r'$\Lambda_\mathrm{core}$')
    kappa_core = aerray(core['kappa2'], u.dimensionless_unscaled,
                        name='kappa_core', limits=[0, 0.002],
                        label=r'$\kappa_2^\mathrm{core}$')
    return create_series(time, lambda_pns, kappa_pns, lambda_core,
                            kappa_core)

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
