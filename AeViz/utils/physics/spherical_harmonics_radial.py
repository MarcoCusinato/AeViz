from __future__ import annotations
from AeViz.spherical_harmonics.spherical_harmonics import SphericalHarmonics
import numpy as np
import scipy.special as sp
from AeViz.utils.utils import (check_existence, progressBar, checkpoints)
from AeViz.utils.files.file_utils import save_hdf
import os, h5py
from AeViz.units import u
from typing import Literal
from AeViz.units import aerray, aeseries
from AeViz.simulation import Simulation

def get_radius_indices(simulation, r, radius):
    if isinstance(r, str):
        rs = r.split('-')
        if len(rs) == 1:
            rr = rs[0]
            rt = 'avg'
            rc = None
        elif len(rs) == 2:
            rr = rs[0]
            rt = rs[1] if rs[1] in ['avg', 'max', 'min'] else 'avg'
            rc = rs[1] if rs[1] not in ['avg', 'max', 'min', 'full'] else None
        elif len(rs) == 3:
            rr = rs[0]
            rt = rs[1] if rs[1] in ['avg', 'max', 'min'] else rs[2] \
                if rs[2] in ['avg', 'max', 'min'] else 'avg'
            rc = rs[1] if rs[1] not in ['avg', 'max', 'min', 'full'] else \
                rs[2] if rs[2] not in ['avg', 'max', 'min', 'full'] else None
        r = getattr(simulation, rr)(rad=rt) if rc is None else \
            getattr(simulation, rr)(rad=rt, comp=rc)
        rindex = np.argmax(radius[:, None] >= r.data[None, :], axis=0)
    else:
        if r >= radius[-1]:
            rindex = -1
        elif r <= radius[0]:
            rindex = 0
        else:
            rindex = np.argmax(radius >= r)
    return rindex

def Harmonics_decomposition_rho(simulation: Simulation,
                                file_name: str,
                                theta: aerray,
                                phi: aerray,
                                dOmega: aerray,
                                SpH: SphericalHarmonics,
                                lmax: int = 4) -> aerray:
    """
    Computes the spherical harmonics decomposition of the density up to
    a value of l equal to lmax for all ms. 

    Parameters
    ----------
    simulation : Simulation
        smiulation of which to compute the spherical harmonics
        decomposition
    file_name : str
        name of the timestep containing the local output.
    theta : aerray
        grid in the theta direction.
    phi : aerray
        grid in the phi direction
    dOmega : aerray
        solid angle element
    SpH : SphericalHarmonics
        class containing the methods to compute the spherical harmonics
    lmax : int, optional
        maximum l to consider, by default 4

    Returns
    -------
    aerray
        contains the decomposition in spherical harmonics of the
        timestep
    """
    rho = simulation.rho(file_name)
    out_array = np.zeros((int(sp.factorial(lmax)) + 1, rho.shape[-1]))
    harm_index = 0
    for l in range( lmax + 1 ):
        for m in range( -l, l + 1 ):
            Ylm = SpH.Ylm_norm(m, l, theta, phi)
            out_array[harm_index, :] = np.sum( rho * Ylm[..., None] * 
                                              dOmega[..., None],
                                        axis=tuple(range(simulation.dim-1)))
            harm_index += 1
    return out_array

def Harmonics_decomposition_rho_msum(simulation: Simulation,
                                     file_name: str,
                                     theta: aerray,
                                     phi: aerray,
                                     dOmega: aerray,
                                     SpH: SphericalHarmonics,
                                     lmax: int = 40):
    """
    Computes the spherical harmonics decomposition of the density up to
    a value of l equal to lmax summed for all ms. 

    Parameters
    ----------
    simulation : Simulation
        smiulation of which to compute the spherical harmonics
        decomposition
    file_name : str
        name of the timestep containing the local output.
    theta : aerray
        grid in the theta direction.
    phi : aerray
        grid in the phi direction
    dOmega : aerray
        solid angle element
    SpH : SphericalHarmonics
        class containing the methods to compute the spherical harmonics
    lmax : int, optional
        maximum l to consider, by default 40

    Returns
    -------
    aerray
        contains the decomposition in spherical harmonics of the
        timestep
    """
    rho = simulation.rho(file_name)
    out_array = np.zeros((lmax+1, rho.shape[-1]))
    harm_index = 0
    for l in range( lmax + 1 ):
        for m in range( -l, l + 1 ):
            Ylm = SpH.Ylm_norm(m, l, theta, phi)
            out_array[harm_index, :] += np.sum( rho * Ylm[..., None] * 
                                               dOmega[..., None],
                                        axis=tuple(range(simulation.dim-1))) ** 2
        harm_index += 1
    return np.sqrt(out_array)
   
def calculate_rho_decomposition(simulation: Simulation,
                                save_checkpoints: bool = True,
                                msum: bool = False,
                                no_new:bool = False) -> bool:
    """
    Computes the density decomposition in spherical harmonics of all
    timestep of a simulation.

    Parameters
    ----------
    simulation : Simulation
        the simultion object to consider
    save_checkpoints : bool, optional
        it will save the result every several timestep depending on the
        simulatin dimensionality, by default True
    msum : bool, optional
        sums over the azimuthal number, by default False
    no_new : bool, optional
        if to compute the rest of the decomposition, by default False

    Returns
    -------
    bool
        True if the computation is finished or it does not have to be
        resumed.
    """
    if msum:
        lmax = 40
        fname = 'rho_decomposition_SpH_msum.h5'
    else:
        lmax = 4
        fname = 'rho_decomposition_SpH.h5'
    if check_existence(simulation, fname):
        time, decomposition, processed_hdf = read_rho_decomposition(simulation, 
                                                                    lmax, msum)
        if processed_hdf[-1].decode("utf-8") == simulation.hdf_file_list[-1] \
            or no_new:
            return True
        else:
            start_point = len(processed_hdf)
            processed_hdf = [ff.decode("utf-8") for ff in processed_hdf]
            print('Checkpoint found for the harmonics decomposition file, ' \
                  'starting from the beginning.\nPlease wait...')
    else:
        start_point = 0
        processed_hdf = []
        print('No checkpoint found for the harmonics decomposition file, ' \
              'starting from the beginning.\nPlease wait...')
    if (checkpoints[simulation.dim] == False) or (not save_checkpoints):
        checkpoint = len(simulation.hdf_file_list)
    else:
        checkpoint = checkpoints[simulation.dim]
    
    ## Set up the spherical harmonics
    SpH = SphericalHarmonics()
    ## Get the angular component
    dOmega = simulation.cell.dOmega(simulation.ghost)
    theta = simulation.cell.theta(simulation.ghost)
    phi = simulation.cell.phi(simulation.ghost)
    findex = start_point
    check_index = 0
    progress_index = 0
    total_points = len(simulation.hdf_file_list) - start_point
    for file in simulation.hdf_file_list[start_point:]:
        progressBar(progress_index, total_points,
                    suffix='Computing spherical harmonics...')
        if msum:
            in_data = (simulation, file, theta,
                                                       phi, dOmega, SpH)
        else:
            in_data = Harmonics_decomposition_rho(simulation, file, theta, phi,
                                                  dOmega, SpH)
        try:
            time = np.concatenate((time, simulation.time(file)))
            decomposition = np.concatenate((decomposition, in_data[..., None]),
                                           axis=-1)
        except Exception as e:
            time = simulation.time(file)
            decomposition = in_data[..., None]
        processed_hdf.append(file)
        if (check_index >= checkpoint) and save_checkpoints:
            print('Checkpoint reached, saving...\n')
            save_decomposition(simulation, decomposition, time, processed_hdf,
                               lmax, msum)
            
            check_index = 0
        check_index += 1
        progress_index += 1
        findex += 1
    print('Computation completed, saving...')
    save_decomposition(simulation, decomposition, time, processed_hdf, lmax,
                       msum)
    return True

def save_decomposition(simulation: Simulation,
                       decomposition: aerray,
                       time: aerray,
                       processed_hdf: list[str],
                       lmax: int,
                       msum: bool) -> None:
    """
    Saves the decomposition in spherical harmonics in a hdf file.

    Parameters
    ----------
    simulation : Simulation
        simulation from which the spherical harmonics decomposition has
        been computed
    decomposition : aerray
        result of the computation
    time : aerray
        time series of the results
    processed_hdf : list[str]
        list of timestep considered.
    lmax : int
        maximum l considered in the computation
    msum : bool
        if the dcomposition has been done summing over m
    """
    keys = ['time']
    quantity = [time]
    if msum:
        file_name = 'rho_decomposition_SpH_msum.h5'
    else:
        file_name = 'rho_decomposition_SpH.h5'
    dec_index = 0
    for l in range(lmax + 1):
        if msum:
            keys.append('rho_l' + str(l))
            quantity.append(decomposition[dec_index, ...])
            dec_index += 1
        else:
            for m in range(-l, l + 1):
                keys.append('rho_l' + str(l) + 'm' + str(m))
                quantity.append(decomposition[dec_index, ...])
                dec_index += 1
    keys.append('processed')
    quantity.append(processed_hdf)
    save_hdf(os.path.join(simulation.storage_path, file_name),
                keys, quantity)
    
def read_rho_decomposition(simulation: Simulation,
                           lmax: int,
                           msum: bool) -> list:
    """
    Reads the density decomposed in spherical hgarmonics.

    Parameters
    ----------
    simulation : Simulation
        simulation from which the spherical harmonics decomposition has
        been computed
    lmax : int
        maximum l considered in the computation
    msum : bool
        if the dcomposition has been done summing over m

    Returns
    -------
    list
        list of the computation results containing the time, 
        decomposition, and procedd timestep.
    """
    if msum:
        fname = 'rho_decomposition_SpH_msum.h5'
        data_dim = lmax + 1
    else:
        fname = 'rho_decomposition_SpH.h5'
        data_dim = int(sp.factorial(lmax)) + 1
    decomposition_data = h5py.File(os.path.join(simulation.storage_path, 
                                            fname), 'r')
    data = [
        (decomposition_data['time'][...] * u.s)
    ]

    dec_data = np.zeros((data_dim,
                         len(simulation.cell.radius(simulation.ghost)),
                         len(decomposition_data['time'][...])))
    dec_index = 0
    for l in range(lmax + 1):
        if msum:
            key = 'rho_l' + str(l)
            dec_data[dec_index, ...] = decomposition_data[key][...]
            dec_index += 1
        else:
            for m in range(-l, l + 1):
                key = 'rho_l' + str(l) + 'm' + str(m)
                dec_data[dec_index, ...] = decomposition_data[key][...]
                dec_index += 1
    data.append(dec_data)
    data.append(decomposition_data['processed'][...])
    decomposition_data.close()
    return data
    
def get_sph_profile(simulation: Simulation,
                    l: int,
                    m: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    """
    Reads from an hdf file the decomposition

    Parameters
    ----------
    simulation : Simulation
        simulation from which the spherical harmonics decomposition has
        been computed
    l : int
        the l number to consider
    m : int | None, optional
        the m number to consider, if None returns the sum over them,
        by default None

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        list of time and profile
    """
    if m is None:
        fname = 'rho_decomposition_SpH_msum.h5'
        key = 'rho_l' + str(l)
    else:
        fname = 'rho_decomposition_SpH.h5'
        key = 'rho_l' + str(l) + 'm' + str(m)
    decomposition_data = h5py.File(os.path.join(simulation.storage_path, 
                                            fname), 'r')
    data = decomposition_data[key][...]
    time = decomposition_data['time'][...]
    decomposition_data.close()
    return time, data

def get_sph_profiles_r(simulation: Simulation,
                       l: int,
                       m: int = None,
                       zero_norm: bool = True,
                       rhomin: aerray | None = None,
                       rhomax: aerray | None = None,
                       r: aerray | None = None,
                       mode: Literal['mass', 'radius'] = 'radius') -> \
                           tuple[np.ndarray, np.ndarray]:
    """
    Returns the spherical profile or the value of a certain radius or 
    the average over a certain region enclosed by two values of density.

    Parameters
    ----------
    simulation : Simulation
        simulation from which the spherical harmonics decomposition has
        been computed 
    l : int
        the l number to consider
    m : int | None, optional
        the m number to consider, if None returns the sum over them,
        by default None
    zero_norm : bool, optional
        normalised by a_{00}, by default True
    rhomin : aerray | None, optional
        minimum value of density to consider, by default None
    rhomax : aerray | None, optional
        maximum value of density to consider, by default None
    r : aerray | None, optional
        radius at which to extract the spherical harmonics,
        by default None
    mode : Literal['mass', 'radius'], optional
        integration mode, either by mass or radius, by default 'radius'

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        time and profile.
    """
    rr = [rhomin, rhomax, r]
    assert rr.count(None) < 3, "Please provide at least one of the three " \
        "arguments: rhomin, rhomax, r"
    if m is None:
        time, r00 = get_sph_profile(simulation, 0)
    else:
        time, r00 = get_sph_profile(simulation, 0, 0)
    _, rlm = get_sph_profile(simulation, l, m)
    if zero_norm:
        rlm /= r00
    if r is not None:
        radius = simulation.cell.radius(simulation.ghost)
        rindex = get_radius_indices(simulation, r, radius)
        try:
            return time, rlm[rindex, np.arange(len(rindex))]
        except:
            return time, rlm[rindex, ...]
    else:
        rho = simulation.radial_profile('rho').data.value
        if rhomin is None:
            rhomin = 0
        if rhomax is None:
            rhomax = rho.max()
        mask = (rho >= rhomin) & (rho <= rhomax)        
        rlm[~mask] = np.nan
        ## Average over the selected region
        if mode == 'radius':
            dr = simulation.cell.dr_integration(simulation.ghost).value[:, None] * \
                np.ones(rlm.shape)
            rlm = rlm * dr
            dr[~mask] = np.nan
            rlm = np.nansum(rlm, axis=0) / np.nansum(dr, axis=0)
            rlm = np.nan_to_num(rlm)
        elif mode == 'mass':
            rlm = rlm * rho
            rho[~mask] = np.nan
            rlm = np.nansum(rlm, axis=0) / np.nansum(rho, axis=0)
            rlm = np.nan_to_num(rlm)
        return time, rlm
    
def get_data_for_barcode(simulation: Simulation,
                         lmax: int | None = None,
                         lmin: int | None = None,
                         rhomin: aerray | None = None,
                         rhomax:aerray | None = None,
                         r: aerray | None = None,
                         msum: bool = False,
                         zero_norm: bool = True,
                         mode: Literal['mass', 'radius'] = 'radius') -> \
                             tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns the evolution of the l and m over time for a specific radius
    or an average over two values of density.

    Parameters
    ----------
    simulation : Simulation
        simulation from which the spherical harmonics decomposition has
        been computed 
    lmax : int | None, optional
        maximum value of l to consider, if None default values are 
        considered, by default None
    lmin : int | None, optional
        minimum value of l to consider, if None default values are 
        considered, by default None
    rhomin : aerray | None, optional
        minimum value of density to consider, by default None
    rhomax : aerray | None, optional
        maximum value of density to consider, by default None
    r : aerray | None, optional
        radius at which to extract the spherical harmonics,
        by default None
        _description_, by default None
    msum : bool, optional
        if to return only the value of the ls or also the ms, by default
        False
    zero_norm : bool, optional
        normalised by a_{00}, by default True
    mode : Literal['mass', 'radius'], optional
        integration mode, either by mass or radius, by default 'radius'

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        contains, time, l numers, the profiles
    """

    if lmax is None and msum:
        lmax = 40
    elif lmax is None:
        lmax = 4
    if lmin is None:
        lmin = 0
    if not msum:
        Yscale = np.arange(lmin, int(sp.factorial(lmax)) + 1)
    else:
        Yscale = np.arange(lmin, lmax + 1)
    if msum:
        for l in range(lmin, lmax + 1):
            time, rlm = get_sph_profiles_r(simulation, l=l, m=None,
                                           zero_norm=zero_norm, rhomin=rhomin,
                                           rhomax=rhomax, r=r, mode=mode)
            if l == lmin:
                data = rlm[None, ...]
            else:
                data = np.concatenate((data, rlm[None, ...]), axis=0)
    else:
        for l in range(lmin, lmax + 1):
            for m in range(-l, l + 1):
                time, rlm = get_sph_profiles_r(simulation, l=l, m=m,
                                               zero_norm=zero_norm,
                                               rhomin=rhomin, rhomax=rhomax,
                                               r=r, mode=mode)
                if l == lmin and m == -l:
                    data = rlm[None, ...]
                else:
                    data = np.concatenate((data, rlm[None, ...]), axis=0)
    return time, Yscale, data
    
def Fourier_amplitude(simulation: Simulation,
                      save_checkpoints: bool = True,
                      no_new: bool = False) -> bool:
    """
    Computes the fourier amplitude for the first 11 ms in a 3D simulation

    Parameters
    ----------
    simulation : Simulation
        simulation from which the Fourier harmonics decomposition has
        been computed
    save_checkpoints : bool, optional
        it will save the result every several timestep depending on the
        simulatin dimensionality, by default True
    no_new : bool, optional
        if to compute the rest of the decomposition, by default False

    Returns
    -------
    bool
        True if the computation is finished or it does not have to be
        resumed.
    """
    if check_existence(simulation, 'rho_fourier.h5'):
        time, rhom_series, processed_hdf = read_rho_fourier(simulation)
        if processed_hdf[-1].decode("utf-8") == simulation.hdf_file_list[-1] \
            or no_new:
            return True
        else:
            start_point = len(processed_hdf)
            processed_hdf = [ff.decode("utf-8") for ff in processed_hdf]
            print('Checkpoint found for the Fourier coefficients file, ' \
                  'starting from the checkpoint.\nPlease wait...')
    else:
        start_point = 0
        processed_hdf = []
        print('No checkpoint found for the Fourier coefficients file, ' \
              'starting from the beginning.\nPlease wait...')
    if (checkpoints[simulation.dim] == False) or (not save_checkpoints):
        checkpoint = len(simulation.hdf_file_list)
    else:
        checkpoint = checkpoints[simulation.dim]
    dtheta = simulation.cell.dtheta_integration(simulation.ghost).value
    N_theta = len(dtheta) // 2
    dtheta = dtheta[None, N_theta-2:N_theta+2, None]
    theta_norm = dtheta.sum()
    phi = simulation.cell.phi(simulation.ghost).value[:, None, None]
    dphi = simulation.cell.dphi(simulation.ghost).value[:, None, None]
    findex = start_point
    check_index = 0
    progress_index = 0
    ## Compute all the stuff we can just one time
    mexp = {}
    for m in range(11):
        mexp[m] = np.exp(1.j * m * phi) * dphi / theta_norm * dtheta
    total_points = len(simulation.hdf_file_list) - start_point
    for file in simulation.hdf_file_list[start_point:]:
        progressBar(progress_index, total_points,
                    suffix='Computing Fourier coefficients')
        rho = simulation.rho(file).value[:, N_theta-2:N_theta+2, :]
        ## Compute the radial m coefficients
        rhom = {}
        for m in range(11):
            rhom[m] = np.sum(mexp[m] * rho, axis=(0, 1))
        tm = simulation.time(file)
        try:
            time = np.concatenate((time, tm))
            for m in range(11):
                rhom_series[m] = np.concatenate((rhom_series[m],
                                                 rhom[m][..., None]), axis=-1)
        except Exception as e:
            print(e)
            time = tm
            rhom_series = {}
            for m in range(11):
                rhom_series[m] = rhom[m][..., None]

        processed_hdf.append(file)
        if (check_index >= checkpoint and save_checkpoints):
            print('Checkpoint reached, saving...\n')
            save_hdf(os.path.join(simulation.storage_path, 'rho_fourier.h5'),
                     ['time', 'Pm', 'processed'],
                     [time, rhom_series, processed_hdf])
            check_index = 0
        check_index += 1
        progress_index += 1
    print('Computation complete, saving...\n')
    save_hdf(os.path.join(simulation.storage_path, 'rho_fourier.h5'),
                     ['time', 'Pm', 'processed'],
                     [time, rhom_series, processed_hdf])
    return True

def read_rho_fourier(simulation: Simulation) -> list[aerray]:
    """
    Reads the computed Fourier coefficients from a hdf file.

    Parameters
    ----------
    simulation : Simulation
        simulation from which the Fourier decomposition has been
        computed 

    Returns
    -------
    list[aerray]
        list of data
    """
    fourier_data = h5py.File(os.path.join(simulation.storage_path, 
                                                'rho_fourier.h5'), 'r')
    data = [
        (fourier_data['time'][...] * u.s)
    ]

    Pm = {}
    for m in range(11):
        Pm[m] = fourier_data[f'Pm/{m}'][...]
    data.append(Pm)
    data.append(fourier_data['processed'][...])
    fourier_data.close()
    return data

def get_rho_fourier(simulation: Simulation,
                    m: int,
                    mode: Literal['phase', 'amplitude'] = 'amplitude',
                    r: aerray | None = None,
                    zero_norm: bool = True) -> aeseries:
    """
    Returns the evolution of the Fourier coefficient with azimuthal
    number m at a certain radius r

    Parameters
    ----------
    simulation : Simulation
        simulation from which the Fourier decomposition has been
        computed 
    m : int
        azimuthal number
    mode : Literal['phase', 'amplitude'], optional
        if the amplitude of the coefficient or its phase has to be
        returned, by default 'amplitude'
    r : aerray | None, optional
        radius at whcih to extract the coefficient. If None returns the
        profile, by default None
    zero_norm : bool, optional
        if the coeffiecient has to be normalised by a_0, by default True

    Returns
    -------
    aeseries
        series containing the time, data and radius (if a profile is
        returned).
    """
    time, Pms, _ = read_rho_fourier(simulation)
    radius = simulation.cell.radius(simulation.ghost)
    time.set(name='time', label=r'$t-t_\mathrm{b}$', cmap=None, log=False,
             limits=[-0.005, time.value.max()])
    P0 = np.abs(Pms[0])
    Pm = Pms[m]
    if mode == 'phase':
        outdata = aerray(np.angle(Pm), u.radian, name=f'phase_{m}',
                    label=f'$\\phi_{m}$', cmap='rainbow', limits=[-np.pi, np.pi])
    elif mode == 'amplitude':
        if zero_norm:
            outdata = aerray(np.abs(Pm)/P0, u.dimensionless_unscaled, name=f'ampl_{m}',
                    label=r'$\tilde{P}_{%d}/\tilde{P}_{0}$' % m, cmap='cividis',
                    limits=[(np.abs(Pm)/P0).min() * 1.1, (np.abs(Pm)/P0).max() * 0.9])
        
        else:
            outdata = aerray(np.abs(Pm), u.dimensionless_unscaled, name=f'ampl_{m}',
                    label=r'$\\tilde{P}'+f'_{m}$', cmap='cividis',
                    limits=[np.abs(Pm).min() * 1.1, np.abs(Pm).max() * 0.9])
    if r is not None:
        rindex = np.argmax(radius >= r)
        outdata = outdata[rindex, ...]
        return aeseries(outdata, time=time)
    else:
        return aeseries(outdata, time=time, radius=radius)
        
    