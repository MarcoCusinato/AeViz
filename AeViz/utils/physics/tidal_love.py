
from __future__ import annotations
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
from AeViz.simulation import Simulation

"""
Everything here is taken from the following paper:
https://doi.org/10.1103/PhysRevD.81.123016
In particular, equation (12) for the linearized metric, equation (14)
for the Love number.
"""

def to_cactus_len(x: aerray) -> aerray:
    """
    Convert the radius to cactus units

    Parameters
    ----------
    x : aerray
        radius

    Returns
    -------
    aerray
    """
    return x / c.G / c.Msol * c.c ** 2

def to_cactus_time(x: aerray) -> aerray:
    """
    Convert the time to cactus units

    Parameters
    ----------
    x : aerray
        time

    Returns
    -------
    aerray
    """
    return x / c.G / c.Msol * c.c ** 3

def to_cactus_dens(x: aerray) -> aerray:
    """
    Convert the density to cactus units

    Parameters
    ----------
    x : aerray
        density

    Returns
    -------
    aerray
    """
    return x / c.Msol * c.G ** 3  * c.Msol ** 3 / c.c ** 6

def to_cactus_vel(x: aerray) -> aerray:
    """
    Convert the velocity to cactus units

    Parameters
    ----------
    x : aerray
        velocity

    Returns
    -------
    aerray
    """
    return x / c.c

def to_cactus_pres(x: aerray) -> aerray:
    """
    Convert the pressure to cactus units

    Parameters
    ----------
    x : aerray
        pressure

    Returns
    -------
    aerray
    """
    return to_cactus_dens(x) / c.c ** 2

def linearized_metric(r: aerray,
                      y: list[aerray],
                      p: aerray,
                      m: aerray,
                      cs: aerray,
                      rho: aerray) -> list[aerray]:
    """
    Computes the linearised metric as in equation 12 of Hinderer et al,
    2010, https://doi.org/10.1103/PhysRevD.81.123016.
    
    .. math::
        \frac{d\beta}{dr}
            =
            2\left(1-\frac{2m_r}{r}\right)^{-1} H
            \left\{
            -2\pi \left[5\epsilon + 9p + f(\epsilon + p)\right]
            + \frac{3}{r^2}
            + 2\left(1-\frac{2m_r}{r}\right)^{-1}
            \left(\frac{m_r}{r^2}+4\pi r p\right)^2
            \right\}
            +\frac{2\beta}{r}
            \left(1-\frac{2m_r}{r}\right)^{-1}
            \left\{
            -1+\frac{m_r}{r}+2\pi r^2(\epsilon-p)
            \right\}.

    Parameters
    ----------
    r : aerray
        radius
    y : list[aerray]
        .. math::
            \frac{dH}{dr} \qquad H
    p : aerray
        pressure
    m : aerray
        mass
    cs : aerray
        soundspeed
    rho : aerray
        density

    Returns
    -------
    list[aerray]
        Returns the first and second dericative of H
    """
    dH_dr, H = y
    if r == 0:
        return [y[0], 0]
    
    coeff = 1 / (1 - 2 * m(r) / r)

    d2H_dr2 = 2 * coeff * H * \
        (-2 * np.pi * (5 * rho(r) + 9 * p(r) + (rho(r) + p(r)) / cs(r) ** 2) +
         3 / r ** 2 + 2 * coeff * 
         (m(r) / r ** 2 + 4 * np.pi * r * p(r)) ** 4) + \
         2 * dH_dr / r * coeff * (-1 + m(r) / r + 2 * np.pi * r ** 2 * 
                                  (rho(r) - p(r)))
    return [d2H_dr2, dH_dr]

def Love_number(xi: aerray,
                y: aerray) -> aerray:
    """
    Compute the love number for a given compacness and
    R H'(R)/H(R) = y, as of equation 14 of Hinderer et al 2010.
    .. math::
        k_2 = \frac{8C^5}{5(1-2C)^2}
            \left[
            2 + 2C(y-1) - y
            \left\{
            2C\left[6 - 3y + 3C(5y-8)\right]
            + 4C^3\left[13 - 11y + C(3y-2) + 2C^2(1+y)\right]
            + 3(1-2C)^2\left[2 - y + 2C(y-1)\right]\ln(1-2C)
            \right\}
            \right]^{-1}

    Parameters
    ----------
    xi : aerray
        compactness
    y : aerray
        .. math::
            R \frac{dH}{dr}\frac{1}{H}

    Returns
    -------
    aerray
        returns the tidal love number
    """
    XI = (1 - 2 * xi)
    kappa2 = 8 / 5 * xi ** 5 * XI ** 2 * (2 + 2 * xi * (y- 1) - y) / \
          (2 * xi * (6 - 3 * y + 3 * xi * (5 * y - 8)) + 4 * xi ** 3 * \
           (13 - 11 * y + xi * (3 * y - 2) + 2 * xi ** 2 * (1 + y)) + \
            3 * XI ** 2 * (2 - y+ 2 * xi * (y- 1)) * (np.log(XI)))
    return kappa2

def tidal_deformability(kappa2: aerray,
                        xi: aerray) -> aerray:
    """
    Compute the tidal deformability from the love number

    Parameters
    ----------
    kappa2 : aerray
        tidal love number
    xi : aerray
        compactness

    Returns
    -------
    aerray
        tidal deformability of the star
    """
    return 2 / 3 * kappa2 / xi ** 5

def solve_tidal_love(xi: aerray,
                     pres: aerray,
                     mass: aerray,
                     soundspeed: aerray,
                     dens: aerray,
                     radius: aerray) -> list[aerray]:
    """
    Computes the tidal Love number for a given compactness xi

    Parameters
    ----------
    xi : aerray
        compactness
    pres : aerray
        pressure
    mass : aerray
        enclosed mass
    soundspeed : aerray
        speed of sound of the matter
    dens : aerray
        density
    radius : aerray
        radius

    Returns
    -------
    list[aerray]
        list contining the tidal love number and tidal deformability
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

def solve_tidal_love_profile(simulation: Simulation,
                             save_checkpoints: bool = True) -> list[aeseries]:
    """
    Derives and saves the tidal love number and tidal deformability for
    the PNS and PNS core.

    Parameters
    ----------
    simulation : Simulation
        simulation from which to compute the tidal deformability and
        love number from
    save_checkpoints : bool, optional
        If true will save checpoints every several timesteps,
        by default True

    Returns
    -------
    list[aeseries]
        contains the evolution of the PNS and PNS core tidal number.
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
    dr = simulation.cell.dr_integration(simulation.ghost)
    ## Get the calculated profiles
    rho_prof = simulation.radial_profile('rho')
    pgas_prof = simulation.radial_profile('gas_pressure')
    cs_prof = simulation.radial_profile('soundspeed')
    ## Convert to cactus units
    t = PNS_radius.time
    radius = to_cactus_len(rho_prof.radius)
    mass_prof = np.cumsum(rho_prof.data * 4 * np.pi * dr[None, :],
                          axis=0).to(u.Msun) / u.Msun
    rho_prof = to_cactus_dens(rho_prof.data)
    pgas_prof = to_cactus_pres(pgas_prof.data)
    cs_prof = to_cactus_vel(cs_prof.data)
    PNS_radius = to_cactus_len(PNS_radius.data)
    core_radius = to_cactus_len(core_radius.data)
    ## Indices
    findex = start_point
    check_index = 0
    progress_index = 0
    total_points = len(t) - start_point
    flist = simulation.hdf_file_list
    for findex in range(start_point, len(t)):
        progressBar(progress_index - start_point, total_points,
                    suffix='Computing tidal quantities...')
        if t[findex] < 0:
            tidal_core, tidal_pns = 0, 0
            love_core, love_pns = 0, 0
        else:
            ##compute the speed of sound profile
            csound = cs_prof[:, findex]
            ##compute the mass profile
            mass = mass_prof[:, findex]
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
        processed_hdf.append(flist[findex])
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

def read_tidal(simulation: Simulation) -> np.ndarray:
    """
    Reads the tidal and love number data from and hdf file.

    Parameters
    ----------
    simulation : Simulation
        simulation object from which the data where computed.

    Returns
    -------
    np.ndarray
        tidal data
    """
    data_file = h5py.File(os.path.join(simulation.storage_path, 'tidal.h5'),
                          'r')
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
