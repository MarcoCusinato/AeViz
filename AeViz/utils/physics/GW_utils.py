from __future__ import annotations
import numpy as np
from AeViz.utils.math_utils import IDL_derivative, gradient
from AeViz.utils.files.string_utils import merge_strings
from numpy.fft import fft, fftfreq
import os, h5py
from AeViz.utils.utils import check_existence, progressBar, checkpoints
from AeViz.utils.files.file_utils import save_hdf, create_series
from AeViz.spherical_harmonics.spherical_harmonics import SphericalHarmonics
from AeViz.units.aeseries import aeseries, aerray
from AeViz.units import u
from AeViz.units.constants import constants as c
from AeViz.GWs import GWstrain
from AeViz.simulation import Simulation
from typing import Literal


## ---------------------------------------------------------------------
## GW strain
## ---------------------------------------------------------------------

def GW_strain(sim_dim: int,
              column_change:int | None,
              data: np.ndarray,
              index: int | None,
              ref: int,
              distance: aerray | None,
              tob: aerray) -> GWstrain:
    """
    computes the gws from a log file of an Aenus CCSN simulation

    Parameters
    ----------
    sim_dim : int
        dimension of the simulation
    column_change : int | None
        row where the number of columns increases or reduces
    data : np.ndarray
        array containing time and quadrupolar tensor components
    index : int | None
        index indicating the end of the tail to compute the GW displacement
    ref : int
        sampling of the GW signal
    distance : aerray | None
        distance at which to simulate the waveform
    tob : aerray
        time of bounce

    Returns
    -------
    GWstrain
        _description_
    """
    assert sim_dim in [1, 2, 3], "Simulation MUST be 1, 2 or 3D."
    if distance:
        if not isinstance(distance, aerray):
            distance = distance * u.cm
    if sim_dim == 1:
        return GW_strain_1D(data)
    elif sim_dim == 2:
        return GWstrain(
            correct_zero(2, GW_strain_2D(data[::ref, :]), index),
            u.s,
            u.cm,
            distance,
            False,
            tob
        )
    else:
        GWs, comps = GW_strain_3D(data)
        if column_change is not None:    
            GWs[:column_change, 1] = \
                GW_strain_2D(data)[:column_change, 1]
            GWs[:column_change-1, 2:] = np.zeros((column_change-1, 3))
            GWs = match_remap(remove_3D_spikes(GWs, column_change),
                              column_change)
            GWs = correct_zero(3, GWs[::ref, :], index)
            comps = comps[::ref, :]
        return GWstrain(
            GWs,
            u.s,
            u.cm,
            distance,
            False,
            tob,
            **{
              'txx': comps[:, 0],
              'tyy': comps[:, 1],
              'tzz': comps[:, 2],
              'txy': comps[:, 3],
              'txz': comps[:, 4],
              'tyz': comps[:, 5]
            }
        )

def GW_strain_1D(data):
    print("No GW for you :'(")
    return None

def GW_strain_2D(data: np.ndarray) -> np.ndarray:
    """
    Computes the GW strain for a 2D simulation computed as the time
    derivative of the matter NE220.
        
    Parameters
    ----------
    data : np.ndarray
        File containing the evolution of the matter qudrupolar moment
        tensor
    Returns
    -------
    np.ndarray
        array containing in the first column the time and second the
        GW strain eveolution
    """
    const = -0.125 *  np.sqrt(15/np.pi)
    GWs = np.stack((data[:, 2], const * IDL_derivative(data[:,2], data[:,5])),
                    axis=-1)
    return GWs

def GW_strain_3D(data: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    r"""
    Computes the GW strain as the first time detivative of the
    quadrupolar tensor component as
    
    .. marh::
        h^{\mathrm{eqtr}}_{+} = 2 (h^{1}_{zz} - h^{1}_{yy}) \\
        h^{\mathrm{pole}}_{+} = 2 (h^{1}_{xx} - h^{1}_{yy}) \\
        h^{\mathrm{eqtr}}_{\times} = - 2 (h^{1}_{yz} + h^{1}_{zy}) \\
        h^{\mathrm{pole}}_{\times} = 2 (h^{1}_{xy} + h^{1}_{yx})  \\
    
    where
    .. math::
        h_{ij}= \frac{1}{D} \frac{\mathrm{d}}{\mathrm{d}t} t^1_{ij}
    
    It also returns the hij components
        

    Parameters
    ----------
    data : np.ndarray
        File containing the evolution of the matter qudrupolar moment
        tensor

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        contains the GW strains and the array of the second derivatives
        of the quadrupolar moment
    """
    hxx = IDL_derivative(data[:, 2], data[:, 9])
    hxy = IDL_derivative(data[:, 2], data[:, 10])
    hxz = IDL_derivative(data[:, 2], data[:, 11])
    hyy = IDL_derivative(data[:, 2], data[:, 13])
    hyx = IDL_derivative(data[:, 2], data[:, 12])
    hyz = IDL_derivative(data[:, 2], data[:, 14])
    hzz = IDL_derivative(data[:, 2], data[:, 17])
    hzy = IDL_derivative(data[:, 2], data[:, 16])
    
    GWs = np.stack((
        data[:, 2],
        2 * (hzz - hyy),
        2 * (hxx - hyy),
        -2 * (hyz + hzy),
        2 * (hxy + hyx)
    ), axis=-1)
    
    comp = np.stack((
        hxx, hyy, hzz,
        hxy, hxz, hyz
    ), axis=-1)
    return GWs, comp

def correct_zero(sim_dim: int,
                 GWs: np.ndarray,
                 index: int | None) -> np.ndarray:
    """
    Some GWs display aconstant displacement from 0 due to high negative
    velocities at the outer boundary of the computational domain.
    This function aims at solve this problem by simply shifting the
    strain.

    Parameters
    ----------
    sim_dim : int
        dimensionality of the simulation
    GWs : np.ndarray
        data oredered as time, h+eq, h+pol, hxeq, hxpol
    index : int | None
        the index upon which find the diplacement

    Returns
    -------
    np.ndarray
        corrected GW signal
    """
    
    if sim_dim == 1:
        pass
    else:
        if index is None:
            return GWs
        if sim_dim == 2:
            GWs[:, 1] -= GWs[:index, 1].mean()
        else:
            for igws in range(1, GWs.shape[1]):
                GWs[:, igws] -= GWs[:index, igws].mean()
        return GWs

def remove_3D_spikes(GWs: np.ndarray,
                     index: int) -> np.ndarray:
    """
    Around the 2D-3D matching there are spikes at least in some cases,
    let's try to remove them.

    Parameters
    ----------
    GWs : np.ndarray
        the GW strain with all polarisations and los
    index : int
        index of remapping

    Returns
    -------
    np.ndarray
        returns the GWs without spikes (hopefully)
    """
    for igws in range(1, GWs.shape[1]):
        for i in range(index - 50, index + 50):
            if GWs[i, igws] >= 1e3 or GWs[i, igws] <= -1e3:
                GWs[i, igws] = GWs[i-1, igws]            
    return GWs

def match_remap(GWs: np.ndarray,
                index: int) -> np.ndarray:
    """
    In some cases the pre bounce 2d evolution and the post 2d evolution
    of the GW signal are not matched because of high velocities at the
    outer boundary of the computational domain. We  try to match the two.

    Parameters
    ----------
    GWs : np.ndarray
        array containing the polarisations and the time
    index : int
        remap index

    Returns
    -------
    np.ndarray
        matched polarisations and time
    """
    for igws in range(1, GWs.shape[1]):
        diff = np.mean(GWs[index-8:index, igws]) - \
            np.mean(GWs[index:index+8, igws])
        GWs[index:, igws] += diff
    return GWs

def universal_modes_relation(PNS_mass, PNS_radius,
                             mode:Literal['2f_torres', '2p1_torres',
                                          '2p2_torres', '2p3_torres', 
                                          '2g1_torres', '2g2_torres',
                                          '2g3_torres'],
                             rhoC=None, pC=None):
    """
    Universal relations coming from:
        Torres-Forné+18 https://arxiv.org/pdf/1902.10048
        Sotani+21
    """
    modes = {
        '2f_torres':{'a': 0, 'b': 1.41e5, 'c': -4.23e6, 'd': 0,
                     'mexp': 0.5, 'rexp': 3/2, 'nm':'2f', 'lb':r'$^2f$'},
        '2p1_torres':{'a': 0, 'b': 2.205e5, 'c': 4.63e6, 'd': 0,
                     'mexp': 0.5, 'rexp': 3/2, 'nm':'2p1', 'lb':r'$^2p_1$'},
        '2p2_torres':{'a': 0, 'b': 4.02e5, 'c': 7.4e6, 'd': 0,
                     'mexp': 0.5, 'rexp': 3/2, 'nm':'2p2', 'lb':r'$^2p_2$'},
        '2p3_torres':{'a': 0, 'b': 6.21e5, 'c': -1.9e6, 'd': 0,
                     'mexp': 0.5, 'rexp': 3/2, 'nm':'2p3', 'lb':r'$^2p_3$'},
        '2g1_torres':{'a': 0, 'b': 8.67e5, 'c': -51.9e6, 'd':0,
                     'mexp': 1, 'rexp': 2, 'nm':'2g1', 'lb':r'$^2g_1$'},
        '2g2_torres':{'a': 0, 'b': 5.88e5, 'c': -86.2e6, 'd': 4.67e10,
                     'mexp': 1, 'rexp': 2, 'nm':'2g2', 'lb':r'$^2g_2$'},
        '2g3_torres':{'a': 905, 'b': -79.9, 'c': -11000, 'd': 0,
                     'mexp': 0.5, 'rexp': 3/2, 'nm':'2g3', 'lb':r'$^2g_3$'},          
    }
    md = modes[mode]
    x = PNS_mass.data.to(u.Msun).value ** md['mexp'] / \
        PNS_radius.data.to(u.km).value ** md['rexp']
    if mode == '2g3_torres':
        x = x * pC.data[0, :].value / rhoC.data[0, :].value ** 2.5
    val = md['a'] + md['b'] * x + md['c'] * x ** 2 + md['d'] * x ** 3
    frequency = aerray(val, u.Hz, md['nm'], md['lb'], None,
                       [val.min(), val.max()], False)
    return aeseries(frequency, time=PNS_mass.time.copy())

## ---------------------------------------------------------------------
## GW spectrogram
## ---------------------------------------------------------------------

def GWs_spectrogram(sim_dim: int,
                    GWs: list[GWstrain],
                    window_size: aerray,
                    scale_to: str,
                    **kwargs) -> aeseries | list[aeseries]:
    """
    Computes the spectrograms of the GW strain

    Parameters
    ----------
    sim_dim : int
        dimension of the simulation
    GWs : list[GWstrain]
        GW strains for the line of sights and polarisations
    window_size : aerray
        windon lenght used to perform the STFT
    scale_to : str
        computation of the STFT scaled to magnitud or psd

    Returns
    -------
    aeseries | list[aeseries]
        list or single aeseries containing the frequency, time and
        modulo of the STFT.
    """
    assert sim_dim in [1, 2, 3], "Simulation MUST be 1, 2 or 3D."
    if sim_dim == 1:
        print("And also no spectrogram for you :'(\nごめんなさい")
        return None
    elif sim_dim == 2:
        return GWs.stft(window_size=window_size, scale_to=scale_to, **kwargs)
    else:
        return [h.stft(window_size=window_size, scale_to=scale_to, **kwargs)
                for h in GWs]

## ---------------------------------------------------------------------
## SNR
## ---------------------------------------------------------------------

def compute_SNR(simulation: Simulation,
                detector: str,
                comp: Literal['eq', 'pol'],
                distance: aerray = (u.kpc *10),
                time_range: list[aerray] | None = None,
                mode: Literal['fft', 'energy'] = 'fft',
                evolution: bool = False,
                **kwargs) -> float | aeseries:
    """
    Computes the SNR at the given distance from a detrended, windowed
    and padded gw strain.

    Parameters
    ----------
    simulation : Simulation
        simulation object, allows to get the strains and the detectors
        psd
    detector : str
        name of the detector
    comp : Literal['eq', 'pol']
        line of sight to consider
    distance : aerray, optional
        distance from the source, by default (u.kpc *10)
    time_range : list[aerray] | None, optional
        portion of the waveform to consider, by default None
    mode : Literal['fft', 'energy'], optional
        type of charactyeristic strain to use, by default 'fft'
    evolution : bool, optional
        if the evolution of the SNR has to be returned, by default False
    **kwargs : specifically to be used for the windowing process or
               remove padding or detrending

    Returns
    -------
    tuple[float | aeseries, str]
        last value of the SNR or the cumulative integral plus the
        detector name
    """
    ## get the signal
    if simulation._Simulation__gws is None:
        simulation.GW_Amplitudes() 
    GW_strain = simulation._Simulation__gws.copy()
    if time_range is not None:
        istart = np.argmax(GW_strain.time >= time_range[0])
        istop = np.argmax(GW_strain.time >= time_range[1])
        if istop == 0:
            istop = len(GW_strain.time)
        GW_strain = GW_strain[istart:istop]
    GW_strain.set_distance(distance)
    ## set some default parameters for the kwargs
    kwargs.setdefault('regularise', True)
    kwargs.setdefault('dt', None)
    kwargs.setdefault('n', None)
    kwargs.setdefault('pad', True)
    kwargs.setdefault('pad_value', 0)
    kwargs.setdefault('pad_length', (1*u.s))
    kwargs.setdefault('apply_window', True)
    kwargs.setdefault('window_type', 'hann')
    
    kwargs.setdefault('detrend', True)
    ## Detrend the  GW strain
    if kwargs['detrend']:
        GW_strain.detrend()
    ## set up the fft configuration
    GW_strain.set_fft_config(**kwargs)
    ## get the detector and load it
    ASD = simulation.ASD(detector)
    GW_strain.load_asd(ASD)
    
    return (GW_strain.get_SNR(detector = ASD.data.name,
                              los = comp,
                              mode = mode,
                              evolution = evolution),
            detector)

## ---------------------------------------------------------------------
## GWs strain from the postprocessing
## ---------------------------------------------------------------------

def calculate_h(simulation, D=1, THETA=np.pi/2, PHI=0,
                save_checkpoints=True, **kwargs):
    """
    Calculates the h cross and x from postprocessing quantities for 
    every timestep of a simulation.
    Returns
        2D
            time: array of time step
            AE220: len(radius), len(time) array
            Full_strain: len(time) array
            PNS_nucleus_strain: len(time) array
            Convection_strain: len(time) array
            Oter_innercore_strain: len(time) array
        3D
            time
            [h_+, h_x]:  len(radius), len(time) array
            [h_+, h_x]_full: len(time)
            [h_+, h_x]_nucl: len(time)
            [h_+, h_x]_conv: len(time)
            [h_+, h_x]_out: len(time)
    """
    r1 = kwargs.setdefault("r1", None)
    r2 = kwargs.setdefault("r2", None)
    r3 = kwargs.setdefault("r3", None)
    apply_correction = kwargs.setdefault("apply_correction", True)
    radii = [r1, r2, r3]
    radii = [r for r in radii if r is not None]
    if simulation.dim == 1:
        print("No GWs for you :'(")
        return None
    elif simulation.dim == 2:
        return NE220_2D_timeseries(simulation, save_checkpoints, D, radii,
                                   apply_correction)
    elif simulation.dim == 3:
        return Qdot_timeseries(simulation, save_checkpoints, D, THETA, PHI, radii,
                               apply_correction)

## 2D

def NE220_2D_timeseries(simulation, save_checkpoints, D, radii,
                        apply_correction):
    """
    Calculates the NE220 from density and velocities for every timestep
    of a 2D simulation. It also calculates the full, nucleus, convection
    and outer core contributions to the strain.
    """
    if check_existence(simulation, 'NE220.h5'):
        time, NE220, full_NE220, nuc_NE220, conv_NE220, outer_NE220, \
            NE220_rad_corr, processed_hdf = \
            read_NE220(simulation)
        if processed_hdf is None:
            save_hdf(os.path.join(simulation.storage_path, 'NE220.h5'),
                     ['time', 'NE220', 'full_NE220', 'nucleus_NE220',
                      'convection_NE220', 'outer_NE220', 'NE220_corr', 'processed'],
                     [time, NE220, full_NE220, nuc_NE220, conv_NE220,
                      outer_NE220, NE220_rad_corr,
                      simulation.hdf_file_list[:len(time)]])
            time, NE220, full_NE220, nuc_NE220, conv_NE220, outer_NE220, \
                NE220_rad_corr, processed_hdf = \
            read_NE220(simulation)
        if len(processed_hdf) == 0:
            start_point = 0
            processed_hdf = []
            print("No checkpoint found. Starting from step 0")
        elif processed_hdf[-1].decode("utf-8") == simulation.hdf_file_list[-1] or \
            simulation.no_new:
            return calculate_strain_2D(simulation, D, time,
                                       simulation.cell.radius(simulation.ghost),
                                       NE220, full_NE220, nuc_NE220,
                                       conv_NE220, outer_NE220, NE220_rad_corr,
                                       radii)
        else:
            start_point = len(processed_hdf)
            processed_hdf = [ff.decode("utf-8") for ff in processed_hdf]
            print("Checkpoint found." \
                "Starting from step {}".format(start_point))
    else:
        start_point = 0
        processed_hdf = []
        print("No checkpoint found. Starting from step 0")
    checkpoint = checkpoints[simulation.dim]
    findex = start_point
    check_index = 0
    progress_index = 0
    dV = -simulation.cell.dVolume_integration(simulation.ghost)
    dOmega = simulation.cell.dOmega(simulation.ghost)
    ctheta = np.cos(simulation.cell.theta(simulation.ghost))[:, None]
    inner_rad, igcells = simulation.innercore_radius(rad='full')
    nuc_rad, ngcells = simulation.PNS_nucleus_radius(rad='full')

    inner_rad = inner_rad.data
    nuc_rad = nuc_rad.data

    for file in simulation.hdf_file_list[start_point:]:
        fNE220, ffull, finner, fnuc, fouter, corr = NE220_2D(simulation,
            file, dV, dOmega, ctheta, inner_rad[..., findex], igcells,
            nuc_rad[..., findex], ngcells)
        try:
            time = np.concatenate((time, simulation.time(file)))
            NE220 = np.concatenate((NE220, fNE220[..., None]), axis=-1)
            NE220_rad_corr = np.concatenate((NE220_rad_corr, corr[..., None]),
                                            axis=-1)
            full_NE220 = np.concatenate((full_NE220, ffull))
            nuc_NE220 = np.concatenate((nuc_NE220, fnuc))
            conv_NE220 = np.concatenate((conv_NE220, finner))
            outer_NE220 = np.concatenate((outer_NE220, fouter))
        except Exception as e:
            print(e)
            time = simulation.time(file)
            NE220 = fNE220[..., None]
            NE220_rad_corr = corr[..., None]
            full_NE220 = ffull
            nuc_NE220 = fnuc
            conv_NE220 = finner
            outer_NE220 = fouter
        processed_hdf.append(file)
        if save_checkpoints and check_index == checkpoint:
            save_hdf(os.path.join(simulation.storage_path, 'NE220.h5'),
                     ['time', 'NE220', 'full_NE220', 'nucleus_NE220',
                      'convection_NE220', 'outer_NE220', 'NE220_corr', 'processed'],
                     [time, NE220, full_NE220, nuc_NE220, conv_NE220,
                      outer_NE220, NE220_rad_corr, processed_hdf])
            print("Checkpoint reached, saving...")
            check_index = 0
        check_index += 1
        progress_index += 1
        findex += 1
        progressBar(progress_index, len(simulation.hdf_file_list))
    
    print("Computations done, saving...")
    save_hdf(os.path.join(simulation.storage_path, 'NE220.h5'),
                     ['time', 'NE220', 'full_NE220', 'nucleus_NE220',
                      'convection_NE220', 'outer_NE220', 'NE220_corr', 'processed'],
                     [time, NE220, full_NE220, nuc_NE220, conv_NE220,
                      outer_NE220, NE220_rad_corr, processed_hdf])
    return calculate_strain_2D(simulation, D, time, simulation.cell.radius(simulation.ghost),
                               NE220, full_NE220, nuc_NE220, conv_NE220,
                               outer_NE220, NE220_rad_corr, radii)

def Zha_correction_2D(dOmega, ctheta, r, rho, vr):
    """
    Computes the Zha correction into the strains computed on spherica shells
    """
    P2 = 0.5 * (3 * ctheta ** 2 - 1) * dOmega[:, None]
    
    prod = r[None, :] ** 4 * rho * vr * P2
    return prod

def Zha_surface_correction(r1_ind, r2_ind, Zha_corr):
    if r1_ind is None:
        r1 = 0 * Zha_corr.unit
    else:
        r1 = Zha_corr[np.arange(Zha_corr.shape[0]), r1_ind]
        r1 = r1.sum()
    if r2_ind is None:
        r2 = 0 * Zha_corr.unit
    else:
        r2 = Zha_corr[np.arange(Zha_corr.shape[0]), r2_ind]
        r2 = r2.sum()
    return r2-r1
    
def NE220_2D(simulation, file_name, dV, dOmega, ctheta, inner_rad, igcells,
          nuc_rad, ngcells):
    """
    Calculates the NE220 from density and velocities for for as single
    timestep. SAme process is employed in
    """
    radius = simulation.cell.radius(simulation.ghost)
    rho = simulation.rho(file_name)
    vr = simulation.radial_velocity(file_name)
    vt = simulation.theta_velocity(file_name)
    NE220 = dV * radius * rho * (vr * (3 * ctheta ** 2 - 1) - \
        3 * vt * ctheta * np.sqrt(1 - ctheta ** 2))
    mask_nuc = radius <= \
        simulation.ghost.remove_ghost_cells_radii(nuc_rad, simulation.dim,
                                                    **ngcells)[..., None]
    r_nuc_ind = np.argmax(radius[None, :] >= \
        simulation.ghost.remove_ghost_cells_radii(nuc_rad, simulation.dim,
                                                    **ngcells)[..., None], axis=-1)
    mask_inner = (radius <= \
        simulation.ghost.remove_ghost_cells_radii(inner_rad, simulation.dim, 
                                             **igcells)[..., None] + (2e6 * u.cm)) & \
        (np.logical_not(mask_nuc))
    r_inner_ind = np.argmax(radius[None, :] >= \
        simulation.ghost.remove_ghost_cells_radii(inner_rad, simulation.dim,
                                                    **ngcells)[..., None]+ (2e6 * u.cm), axis=-1)
    r_outer_ind = -np.ones(vt.shape[0], dtype=int)
    mask_outer = np.logical_not(mask_inner + mask_nuc)
    Zha_corr = Zha_correction_2D(dOmega, ctheta, radius, rho, vr)
    nuc_corr = Zha_surface_correction(None, r_nuc_ind, Zha_corr)
    inn_corr = Zha_surface_correction(r_nuc_ind+1, r_inner_ind, Zha_corr)
    out_corr = Zha_surface_correction(r_inner_ind+1, r_outer_ind, Zha_corr)
    return np.sum(NE220, axis=0), np.sum(NE220), np.sum(NE220 * mask_inner) + nuc_corr, \
        np.sum(NE220 * mask_nuc)+inn_corr, np.sum(NE220 * mask_outer)+out_corr, np.sum(Zha_corr, axis=0)
           
def read_NE220(simulation):
    """
    Reads the NE220 from a checkpoint file.
    """
    with h5py.File(os.path.join(simulation.storage_path, 'NE220.h5'), 'r') as data:
        if 'NE220_corr' not in data.keys():
            return 0, 0, 0, 0, 0, 0, 0, []
        time = data['time'][...] * u.s
        NE220 = data['NE220'][...] * u.cm ** 2 * u.g / u.s
        full_NE220 = data['full_NE220'][...] * u.cm ** 2 * u.g / u.s
        nuc_NE220 = data['nucleus_NE220'][...] * u.cm ** 2 * u.g / u.s
        conv_NE220 = data['convection_NE220'][...] * u.cm ** 2 * u.g / u.s
        outer_NE220 = data['outer_NE220'][...] * u.cm ** 2 * u.g / u.s
        correction = data['NE220_corr'][...] * u.cm ** 2 * u.g / u.s
        if 'processed' in data.keys():
            processed = data['processed'][...]
        else:
            processed = None
        
    return time, NE220, full_NE220, nuc_NE220, conv_NE220, outer_NE220, correction, processed

def calculate_strain_2D(simulation, D, time, radius, NE220, full_NE220, nuc_NE220,
                        conv_NE220, outer_NE220, corrections, radii):
    """
    Derives ancd fixes the constants of the strain.
    """
    const =  -0.125 *  np.sqrt(15/np.pi) * \
        (c.G * 8 * np.pi ** 0.5 / (np.sqrt( 15 ) * c.c ** 4))
    if D is not None:
        if not isinstance(D, aerray):
            D = D * u.cm
        const /= D
        add_lb = r''
    else:
        add_lb = r'$\mathcal{D}$'
        D = 1 * u.dimensionless_unscaled
    time.set(name='time', label=r'$t-t_\mathrm{b}$',
             cmap=None, limits=[-0.005, time[-1]])
    NE220_og = NE220.copy()
    NE220 = const * IDL_derivative(time, NE220)
    NE220.set(name='AE220', label=merge_strings(add_lb, r'$A^{E2}_{20}(r)$'),
              cmap='seismic', limits=[-3 / D.value, 3 / D.value])
    full_NE220 = const * IDL_derivative(time, full_NE220)
    full_NE220.set(name='full_NE220', label=merge_strings(add_lb, r'$h_{+,eq}$'),
                   cmap='seismic', limits=[-70 / D.value, 70 / D.value])
    out_strains = create_series(time, full_NE220)
    nuc_NE220 = const * IDL_derivative(time, nuc_NE220)
    nuc_NE220.set(name='nuc_NE220', cmap='seismic', limits=[-70 / D.value, 70 / D.value],
                  label=merge_strings(add_lb, r'$h_{+,\mathrm{eq,core}}$'))
    conv_NE220 = const * IDL_derivative(time, conv_NE220)
    conv_NE220.set(name='conv_NE220', cmap='seismic', limits=[-70 / D.value, 70 / D.value],
                   label=merge_strings(add_lb,r'$h_{+,\mathrm{eq,conv}}$'))
    outer_NE220 = const * IDL_derivative(time, outer_NE220)
    outer_NE220.set(name='outer_NE220', cmap='seismic', limits=[-70 / D.value, 70 / D.value],
                    label=merge_strings(add_lb, r'$h_{+,\mathrm{eq,outer}}$'))
    T, R = np.meshgrid(time, radius)
    ## Now we select the right radius
    if len(radii) == 1:
        if radii[0] == 'PNS_nucleus_radius-full':
            out_strains.extend(create_series(time, nuc_NE220))
        else:
            corr_r1, mask_r1 = get_correction_evolution(simulation, radii[0], radius, corrections)
            out_strains.extend(compute_partial_strain(NE220_og, time, None,
                                                      corr_r1, None, mask_r1,
                                                      const, D, add_lb,
                                                      r'$h_{+,\mathrm{eq,r1}}$'))
    elif len(radii) == 2:
        if radii == ['PNS_nucleus_radius-full', 'innercore_radius-full']:
            out_strains.extend(create_series(time, nuc_NE220, conv_NE220))
        else:
            corr_r1, mask_r1 = get_correction_evolution(simulation, radii[0], radius, corrections)
            out_strains.extend(compute_partial_strain(NE220_og, time, None,
                                                      corr_r1, None, mask_r1,
                                                      const, D, add_lb,
                                                      r'$h_{+,\mathrm{eq,r1}}$'))
            corr_r2, mask_r2 = get_correction_evolution(simulation, radii[1], radius, corrections)
            out_strains.extend(compute_partial_strain(NE220_og, time, corr_r1,
                                                      corr_r2, mask_r1, mask_r2,
                                                      const, D, add_lb,
                                                      r'$h_{+,\mathrm{eq,r2}}$'))
    elif len(radii) == 3:
        if radii[:2] == ['PNS_nucleus_radius-full', 'innercore_radius-full']:
            out_strains.extend(create_series(time, nuc_NE220, conv_NE220,
                                             outer_NE220))
        else:
            corr_r1, mask_r1 = get_correction_evolution(simulation, radii[0], radius, corrections)
            out_strains.extend(compute_partial_strain(NE220_og, time, None,
                                                      corr_r1, None, mask_r1,
                                                      const, D, add_lb,
                                                      r'$h_{+,\mathrm{eq,r1}}$'))
            corr_r2, mask_r2 = get_correction_evolution(simulation, radii[1], radius, corrections)
            out_strains.extend(compute_partial_strain(NE220_og, time, corr_r1,
                                                      corr_r2, mask_r1, mask_r2,
                                                      const, D, add_lb,
                                                      r'$h_{+,\mathrm{eq,r2}}$'))
            corr_r3, mask_r3 = get_correction_evolution(simulation, radii[2], radius, corrections)
            out_strains.extend(compute_partial_strain(NE220_og, time, corr_r2,
                                                      corr_r3, mask_r2, mask_r3,
                                                      const, D, add_lb,
                                                      r'$h_{+,\mathrm{eq,r3}}$'))
            
    return aeseries(NE220, time=T, radius=R), out_strains

def compute_partial_strain(NE220, time, corr1, corr2, mask1, mask2, const, D,
                           Dlab, label):
    if corr1 is None:
        corr = corr2
    else:
        corr = corr2 - corr1
    AE220 = NE220.copy()
    if mask1 is not None:
        if isinstance(mask1, np.int64):
            AE220[:mask1, :] = 0
        else:
            AE220 = NE220.copy()
            AE220[mask1] = 0
    if isinstance(mask2, np.int64):
        AE220[mask2:, :] = 0
    else:
        AE220[~mask2] = 0
    AE220 = np.sum(AE220, axis=0) - corr
    strain = const * IDL_derivative(time, AE220)
    strain.set(name='partial_strain', label=merge_strings(Dlab, label),
                   cmap='seismic', limits=[-70 / D.value, 70 / D.value])
    return create_series(time, strain)

def get_correction_evolution(simulation, r, radius, amplitude):
    ## get the correspondig radius
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
        rindex = np.argmax(radius[None, :] >= r.data[:, None], axis=-1)
        corr = amplitude[rindex, np.arange(amplitude.shape[1])]
        rindex = radius[:, None] * np.ones(r.data.shape)[None, :] <= r.data[None, :]
    elif isinstance(r, float):
        r = r * radius.unit
        rindex = np.argmax(radius >= r)
        corr = amplitude[rindex, :]
    elif isinstance(r, aerray):
        rindex = np.argmax(radius >= r)
        corr = amplitude[rindex, :]
    else:
        raise ValueError("Option not recognized")
    return corr, rindex
## 3D

def Qdot_timeseries(simulation, save_checkpoints, D, THETA, PHI, radii,
                    apply_correction):
    """
    Calculates the NE220 from density and velocities for every timestep
    of a 2D simulation. It also calculates the full, nucleus, convection
    and outer core contributions to the strain.
    """
    if check_existence(simulation, 'Qdot.h5'):
        time, Qdot_radial, Qdot_corr, Qdot_total, Qdot_inner, Qdot_nucleus, Qdot_outer, processed_hdf = \
            read_Qdot(simulation)
        if processed_hdf is None:
            save_hdf(os.path.join(simulation.storage_path, 'Qdot.h5'),
                     ['time', 'Qdot_total', 'Qdot_inner', 'Qdot_nucleus',
                      'Qdot_outer', 'Qdot_radial', 'Qdot_corr', 'processed'],
                     [time, Qdot_total, Qdot_inner, Qdot_nucleus, Qdot_outer,
                      Qdot_corr, Qdot_radial,
                      simulation.hdf_file_list[:len(time)]])
            time, Qdot_radial, Qdot_corr, Qdot_total, Qdot_inner, Qdot_nucleus, Qdot_outer, processed_hdf = \
            read_Qdot(simulation)
        if len(processed_hdf) == 0:
            start_point = 0
            processed_hdf = []
            print("No checkpoint found. Starting from step 0")
        elif processed_hdf[-1].decode("utf-8") == simulation.hdf_file_list[-1] or \
            simulation.no_new:
            return calculate_strain_3D(simulation, D, THETA, PHI, time,
                                       simulation.cell.radius(simulation.ghost),
                                       Qdot_radial, Qdot_total, Qdot_inner,
                                       Qdot_nucleus, Qdot_outer, Qdot_corr, radii,
                                       apply_correction)
        else:
            start_point = len(processed_hdf)
            processed_hdf = [ff.decode("utf-8") for ff in processed_hdf]
            print("Checkpoint found." \
                "Starting from step {}".format(start_point))
    else:
        start_point = 0
        processed_hdf = []
        print("No checkpoint found. Starting from step 0")
    checkpoint = checkpoints[simulation.dim]
    findex = start_point
    check_index = 0
    progress_index = 0
    dV = simulation.cell.dVolume_integration(simulation.ghost)
    dOmega = simulation.cell.dOmega(simulation.ghost)
    inner_rad, igcells = simulation.innercore_radius(rad='full')
    nuc_rad, ngcells = simulation.PNS_nucleus_radius(rad='full')

    inner_rad = inner_rad.data
    nuc_rad = nuc_rad.data
        
    grad, harm = get_spherical_harmonics(
                    simulation.cell.radius(simulation.ghost),
                    simulation.cell.theta(simulation.ghost),
                    simulation.cell.phi(simulation.ghost),
                    dOmega)
    
    for file in simulation.hdf_file_list[start_point:]:
        Qtot, Qinner, Qnuc, Qouter, Qradial, Qcorr = \
            calculate_Qdot(simulation, grad, harm,
                           file, dV, inner_rad[..., findex],
                           igcells, nuc_rad[..., findex], ngcells)
        try:
            time = np.concatenate((time, simulation.time(file)))
            Qdot_radial = np.concatenate((Qdot_radial, Qradial[..., None]),
                                         axis=-1)
            Qdot_total = np.concatenate((Qdot_total, Qtot[..., None]), axis=-1)
            Qdot_inner = np.concatenate((Qdot_inner, Qinner[..., None]),
                                        axis=-1)
            Qdot_nucleus = np.concatenate((Qdot_nucleus, Qnuc[..., None]),
                                          axis=-1)
            Qdot_outer = np.concatenate((Qdot_outer, Qouter[..., None]),
                                        axis=-1)
            Qdot_corr = np.concatenate((Qdot_corr, Qcorr[..., None]),
                                        axis=-1)
        except Exception as e:
            print(e)
            time = simulation.time(file)
            Qdot_radial = Qradial[..., None]
            Qdot_total = Qtot[..., None]
            Qdot_inner = Qinner[..., None]
            Qdot_nucleus = Qnuc[..., None]
            Qdot_outer = Qouter[..., None]
            Qdot_corr = Qcorr[..., None]
        processed_hdf.append(file)
            
        if save_checkpoints and check_index == checkpoint:
            print("Checkpoint reached. Saving...")
            save_hdf(os.path.join(simulation.storage_path, 'Qdot.h5'),
                     ['time', 'Qdot_total', 'Qdot_inner', 'Qdot_nucleus',
                      'Qdot_outer', 'Qdot_radial', 'Qdot_corr', 'processed'],
                     [time, Qdot_total, Qdot_inner, Qdot_nucleus, Qdot_outer,
                      Qdot_radial, Qdot_corr, processed_hdf])
            check_index = 0
        check_index += 1
        progress_index += 1
        findex += 1
        progressBar(progress_index, len(
            simulation.hdf_file_list[start_point:]))
    
    print("Computations done, saving...")
    save_hdf(os.path.join(simulation.storage_path, 'Qdot.h5'),
                     ['time', 'Qdot_total', 'Qdot_inner', 'Qdot_nucleus',
                      'Qdot_outer', 'Qdot_radial', 'Qdot_corr', 'processed'],
                     [time, Qdot_total, Qdot_inner, Qdot_nucleus, Qdot_outer,
                      Qdot_radial, Qdot_corr, processed_hdf])
    return calculate_strain_3D(simulation, D, THETA, PHI, time,
                               simulation.cell.radius(simulation.ghost),
                               Qdot_radial, Qdot_total, Qdot_inner,
                               Qdot_nucleus, Qdot_outer, Qdot_corr, radii,
                               apply_correction)

def Qdot_surface_correction(Qcorr, r1_ind, r2_ind):
    if r1_ind is None:
        r1 = 0 * Qcorr.unit
    else:
        r1 = Qcorr[np.arange(Qcorr.shape[0])[:, None],
                   np.arange(Qcorr.shape[1])[None, :], r1_ind]
        r1 = r1.sum()
    if r2_ind is None:
        r2 = 0 * Qcorr.unit
    else:
        r2 = Qcorr[np.arange(Qcorr.shape[0])[:, None],
                   np.arange(Qcorr.shape[1])[None, :], r2_ind]
        r2 = r2.sum()
    return r2-r1

def get_spherical_harmonics(radius, theta, phi, dOmega):
    """
    Calculates the gradient of the conjugate spherical harmonics times
    the radius for l = 2.
    Returns a list of arrays with dimension (len(phi), lewn(theta),
    len(radius), 3) containing the gradient of the conjugate spherical
    harmonics times the radius from m=-2 to m=2.
    """
    grd = []
    hr = []
    harmonics = SphericalHarmonics()
    for m in range(-2, 3):
        Y2m_r = harmonics.Ylm_conj(m, 2, theta, phi)[..., None] * \
            radius[None, None, :] ** 2
        grd.append(gradient(Y2m_r, radius, theta, phi, 'spherical'))
        hr.append(Y2m_r * radius[None, None, :] ** 2 * dOmega[..., None])    
    return grd, hr

def calculate_Qdot(simulation, gradY, Ylm, file_name, dV, 
                        inner_rad, igcells, nuc_rad, ngcells):
    """
    Calculates the Qdot for the different regions of the star.
    dot{Q} = dV ρ ∇(v Y*_2m)
    """
    radius = simulation.cell.radius(simulation.ghost)
    mask_nuc = radius <= \
        simulation.ghost.remove_ghost_cells_radii(nuc_rad, simulation.dim,
                                                    **ngcells)[..., None]
    mask_inner = (radius <= \
        simulation.ghost.remove_ghost_cells_radii(inner_rad, simulation.dim, 
                                             **igcells)[..., None] + (20*u.km)) & \
        (np.logical_not(mask_nuc))
    mask_outer = np.logical_not(mask_inner + mask_nuc)
    
    r_nuc_ind = np.argmax(radius[None, None, :] >= \
        simulation.ghost.remove_ghost_cells_radii(nuc_rad, simulation.dim,
                                                    **ngcells)[..., None],
                        axis=-1)
    r_inner_ind = np.argmax(radius[None, None, :] >= \
        simulation.ghost.remove_ghost_cells_radii(inner_rad, simulation.dim,
                                         **ngcells)[..., None] + (2e6 * u.cm),
                        axis=-1)
    r_outer_ind = -np.ones(dV.shape[:-1], dtype=int)
    
    rho = simulation.rho(file_name)
    v_r = simulation.radial_velocity(file_name)
    v_t = simulation.theta_velocity(file_name)
    v_p = simulation.phi_velocity(file_name)
    rho_vr = rho * v_r
    rho *= dV
    Qdot = (rho * (v_r * gradY[0][0, ...] + v_t * gradY[0][1, ...] + v_p \
        * gradY[0][2, ...]))
    Qcorr = Ylm[0] * rho_vr
    Qdot_tot = Qdot.sum()[..., None]
    Qdot_inner = (Qdot[mask_inner].sum() - \
        Qdot_surface_correction(Qcorr, None, r_nuc_ind))[..., None]
    Qdot_nuc = (Qdot[mask_nuc].sum() - \
        Qdot_surface_correction(Qcorr, r_nuc_ind, r_inner_ind))[..., None]
    Qdot_outer = (Qdot[mask_outer].sum() - \
        Qdot_surface_correction(Qcorr, r_inner_ind, r_outer_ind))[..., None]
    Qdot_radial = Qdot.sum(axis=(0,1))[..., None]
    Qdot_corr = Qcorr.sum(axis=(0,1))[..., None]
    for i in range(1, 5):
        Qdot = (rho * (v_r * gradY[i][0, ...] + v_t * \
            gradY[i][1, ...] + v_p * gradY[i][2, ...]))
        Qcorr = Ylm[i] * rho_vr
        Qdot_tot = np.concatenate((Qdot_tot, Qdot.sum()[..., None]), axis=-1)
        Qdot_inner = np.concatenate((Qdot_inner,
                                     (Qdot[mask_inner].sum() - \
                                         Qdot_surface_correction(Qcorr,
                                                                 None,
                                                                 r_nuc_ind))
                                     [..., None]), axis=-1)
        Qdot_nuc = np.concatenate((Qdot_nuc, (Qdot[mask_nuc].sum() - \
                                    Qdot_surface_correction(Qcorr, r_nuc_ind,
                                                    r_inner_ind))[..., None]),
                                  axis=-1)
        Qdot_outer = np.concatenate((Qdot_outer, (Qdot[mask_outer].sum() - \
                                    Qdot_surface_correction(Qcorr,
                                                            r_inner_ind,
                                                            r_outer_ind))
                                     [..., None]), axis=-1)
        Qdot_radial = np.concatenate((Qdot_radial, Qdot.sum(axis=(0, 1))
                                      [..., None]), axis=-1)
        Qdot_corr = np.concatenate((Qdot_corr, Qdot.sum(axis=(0, 1))
                                      [..., None]), axis=-1)
            
    return Qdot_tot, Qdot_inner, Qdot_nuc, Qdot_outer, Qdot_radial, Qdot_corr

def read_Qdot(simulation):
    """
    Reads the Qdot and masks from a checkpoint file.
    """
    with h5py.File(os.path.join(simulation.storage_path, 'Qdot.h5')) as data:
        if 'Qdot_corr' not in data.keys():
            return 0, 0, 0, 0, 0, 0, 0, []
        time = data['time'][...] * u.s
        Qdot_radial = data['Qdot_radial'][...] * u.g * u.cm ** 2 / u.s
        Qdot_total = data['Qdot_total'][...] * u.g * u.cm ** 2 / u.s
        Qdot_inner = data['Qdot_inner'][...] * u.g * u.cm ** 2 / u.s
        Qdot_nucleus = data['Qdot_nucleus'][...] * u.g * u.cm ** 2 / u.s
        Qdot_outer = data['Qdot_outer'][...] * u.g * u.cm ** 2 / u.s
        Qcorr = data['Qdot_corr'][...] * u.g * u.cm ** 2 / u.s
        if 'processed' in data.keys():
            processed_hdf = data['processed'][...]
        else:
            processed_hdf = None
    return time, Qdot_radial, Qcorr, Qdot_total, Qdot_inner, Qdot_nucleus, Qdot_outer, processed_hdf

def compute_partial_corrected_Qdotdot(Qdot, Ylm, time, corr1, corr2, mask1, mask2,
                                      apply_correction):
    Qdot_og = Qdot.copy()
    if corr1 is None:
        corr = corr2
    else:
        corr = corr2 - corr1
    if mask1 is not None:
        if isinstance(mask1, np.int64):
            Qdot_og[:mask1, :] = 0
        else:
            Qdot_og[mask1] = 0
    if isinstance(mask2, np.int64):
        Qdot_og[mask2:, :] = 0
    else:
        Qdot_og[~mask2] = 0
    if apply_correction:
        Qdot_og = np.sum(Qdot_og, axis=0) - corr
    else:
        Qdot_og = np.sum(Qdot_og, axis=0)
    return IDL_derivative(time, Qdot_og) * Ylm

def compute_partial_Qdotdot_3D(simulation, radii, radius, Qdot_radial, Ylm,
                               time, corrections, Qout, apply_correction):
    if len(radii) == 1:
        if radii[0] == 'PNS_nucleus_radius-full':
            return []
        else:
            corr_r1, mask_r1 = get_correction_evolution(simulation, radii[0],
                                                        radius, corrections)
            q0 = compute_partial_corrected_Qdotdot(Qdot_radial, Ylm, time,
                                                   None, corr_r1, None, mask_r1,
                                                   apply_correction)
            if len(Qout) == 0:
                Qout = [q0]
            else:
                Qout[0] += q0
    elif len(radii) == 2:
        if radii == ['PNS_nucleus_radius-full', 'innercore_radius-full']:
            return []
        else:
            corr_r1, mask_r1 = get_correction_evolution(simulation, radii[0],
                                                        radius, corrections)
            q0 = compute_partial_corrected_Qdotdot(Qdot_radial, Ylm, time,
                                                   None, corr_r1, None,
                                                   mask_r1, apply_correction)
            corr_r2, mask_r2 = get_correction_evolution(simulation, radii[1],
                                                        radius, corrections)
            q1 = compute_partial_corrected_Qdotdot(Qdot_radial, Ylm, time,
                                                   corr_r1, corr_r2,
                                                   mask_r1, mask_r2,
                                                   apply_correction)
            if len(Qout) == 0:
                Qout = [q0, q1]
            else:
                Qout[0] += q0
                Qout[1] += q1
    elif len(radii) == 3:
        if radii[:2] == ['PNS_nucleus_radius-full', 'innercore_radius-full']:
            return []
        else:
            corr_r1, mask_r1 = get_correction_evolution(simulation, radii[0],
                                                        radius, corrections)
            q0 = compute_partial_corrected_Qdotdot(Qdot_radial, Ylm, time,
                                                   None, corr_r1, None,
                                                   mask_r1, apply_correction)
            corr_r2, mask_r2 = get_correction_evolution(simulation, radii[1],
                                                        radius, corrections)
            q1 = compute_partial_corrected_Qdotdot(Qdot_radial, Ylm, time,
                                                   corr_r1, corr_r2,
                                                   mask_r1, mask_r2,
                                                   apply_correction)
            corr_r3, mask_r3 = get_correction_evolution(simulation, radii[2],
                                                        radius, corrections)
            q2 = compute_partial_corrected_Qdotdot(Qdot_radial, Ylm, time,
                                                   corr_r2, corr_r3,
                                                   mask_r2, mask_r3,
                                                   apply_correction)
            if len(Qout) == 0:
                Qout = [q0, q1, q2]
            else:
                Qout[0] += q0
                Qout[1] += q1
                Qout[2] += q2
    else:
        return []
    return Qout

def return_3D_strains(time, tot_st, cor_st, inn_st, out_st, oth_st, radii):
    out = create_series(time, tot_st)
    if len(radii) == 1:
        if radii[0] == 'PNS_nucleus_radius-full':
            out.extend(create_series(time, cor_st))
        else:
            out.extend(create_series(time, oth_st[0]))
    elif len(radii) == 2:
        if radii == ['PNS_nucleus_radius-full', 'innercore_radius-full']:
            out.extend(create_series(time, cor_st, inn_st))
        else:
            out.extend(create_series(time, oth_st)[0])
    elif len(radii) == 3:
        if radii[:2] == ['PNS_nucleus_radius-full', 'innercore_radius-full']:
            out.extend(create_series(time, cor_st, inn_st, out_st))
        else:
            out.extend(create_series(time, oth_st)[0])
    return out
    
def calculate_strain_3D(simulation, D, THETA, PHI, time, radius, Qdot_radial, Qdot_total,
                        Qdot_inner, Qdot_nucleus, Qdot_outer, corrections, radii,
                        apply_correction):
    if D is not None:
        if not isinstance(D, aerray):
            D = D * u.cm
        const /= D
        add_lb = r''
    else:
        add_lb = r'$\mathcal{D}$'
        D = 1 * u.dimensionless_unscaled
    harmonics = SphericalHarmonics()
    partialQ = []
    for m in range(5):
        Y22m = harmonics.spin_weighted_Ylm(-2, m-2, 2, THETA, PHI)
        partialQ = compute_partial_Qdotdot_3D(simulation, radii, radius,
                                              Qdot_radial[:, m, :], Y22m, time,
                                              corrections[:, m, :], partialQ,
                                              apply_correction)
        Qdot_radial[:, m, :] = IDL_derivative(time, Qdot_radial[:, m, :]) * \
            Y22m
        Qdot_total[m, :] = IDL_derivative(time, Qdot_total[m, :]) * Y22m
        Qdot_inner[m, :] = IDL_derivative(time, Qdot_inner[m, :]) * Y22m
        Qdot_nucleus[m, :] = IDL_derivative(time, Qdot_nucleus[m, :]) * Y22m
        Qdot_outer[m, :] = IDL_derivative(time, Qdot_outer[m, :]) * Y22m
    const = np.sqrt(2/3) * 8 * np.pi * c.G / (D * c.c ** 4 * 5) / u.s ## last u.s accounts for the derivative
    Qdot_radial = const * Qdot_radial.sum(axis=1)
    Qdot_total = const * Qdot_total.sum(axis=0)
    Qdot_inner = const * Qdot_inner.sum(axis=0)
    Qdot_nucleus = const * Qdot_nucleus.sum(axis=0)
    Qdot_outer = const * Qdot_outer.sum(axis=0)
    partialQ = [const * pQ for pQ in partialQ]
    time.set(name='time', label=r'$t-t_\mathrm{b}$', cmap=None,
             limits=[-0.005, time[-1]])
    ## Compute the polarizations
    hplus_radial = Qdot_radial.real
    hcross_radial = -Qdot_radial.imag
    hplus_tot = Qdot_total.real
    hcross_tot = -Qdot_total.imag
    hplus_nuc = Qdot_nucleus.real
    hcross_nuc = -Qdot_nucleus.imag
    hplus_inn = Qdot_inner.real
    hcross_inn = -Qdot_inner.imag
    hcross_out = -Qdot_outer.imag
    hplus_out = Qdot_outer.real
    partialQ_plus = [pQ.real for pQ in partialQ]
    partialQ_cross = [-pQ.imag for pQ in partialQ]
    
    ## Set the labels
    if np.isclose(THETA, np.pi, 0.05):
        hplus_radial.set(name='hpuls_radial_pol',
                         label=merge_strings(add_lb,
                                             r'$h_\mathrm{+,pol}(r)$'),
                         cmap='seismic',
                         limits=[-3 / D.value, 3 / D.value])
        hcross_radial.set(name='hcross_radial_pol',
                         label=merge_strings(add_lb,
                                             r'$h_\mathrm{\times,pol}(r)$'),
                         cmap='seismic',
                         limits=[-3 / D.value, 3 / D.value])
        hplus_tot.set(name='tot_hplus_pol',
                      cmap='seismic',
                      limits=[-70 / D.value, 70 / D.value],
                      label=merge_strings(add_lb, r'$h_{+,\mathrm{pol,tot}}$'))
        hcross_tot.set(name='tot_hcross_pol',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{pol,tot}}$'))
        hplus_nuc.set(name='nuc_hplus_pol',
                      cmap='seismic',
                      limits=[-70 / D.value, 70 / D.value],
                      label=merge_strings(add_lb, r'$h_{+,\mathrm{pol,core}}$'))
        hcross_nuc.set(name='nuc_hcross_pol',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{pol,core}}$'))
        hplus_inn.set(name='inn_hplus_pol',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{+,\mathrm{pol,conv}}$'))
        hcross_inn.set(name='inn_hcross_pol',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{pol,conv}}$'))
        hplus_out.set(name='out_hcross_pol',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{+,\mathrm{pol,out}}$'))
        hcross_out.set(name='out_hcross_pol',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{pol,out}}$'))
        [pQ.set(name=f'r{i}_hplus_pol',
                      cmap='seismic',
                      limits=[-70 / D.value, 70 / D.value],
                      label=merge_strings(add_lb, r'$h_{+,\mathrm{pol,r$',
                                          f'$_{i}$', r'$}}$'))
                      for i, pQ in enumerate(partialQ_plus)]
        [pQ.set(name=f'r{i}_hcross_pol',
                      cmap='seismic',
                      limits=[-70 / D.value, 70 / D.value],
                      label=merge_strings(add_lb, r'$h_{\times,\mathrm{pol,r$',
                                          f'$_{i}$', r'$}}$'))
                      for i, pQ in enumerate(partialQ_cross)]
    elif np.isclose(THETA, np.pi/2, 0.05):
        hplus_radial.set(name='hpuls_radial_eq',
                         label=merge_strings(add_lb,
                                             r'$h_\mathrm{+,eq}(r)$'),
                         cmap='seismic',
                         limits=[-3 / D.value, 3 / D.value])        
        hcross_radial.set(name='hcross_radial_eq',
                         label=merge_strings(add_lb,
                                             r'$h_\mathrm{\times,eq}(r)$'),
                         cmap='seismic',
                         limits=[-3 / D.value, 3 / D.value])
        hplus_tot.set(name='tot_hplus_eq',
                      cmap='seismic',
                      limits=[-70 / D.value, 70 / D.value],
                      label=merge_strings(add_lb, r'$h_{+,\mathrm{eq,tot}}$'))
        hcross_tot.set(name='tot_hcross_eq',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{eq,tot}}$'))
        hplus_nuc.set(name='nuc_hplus_eq',
                      cmap='seismic',
                      limits=[-70 / D.value, 70 / D.value],
                      label=merge_strings(add_lb, r'$h_{+,\mathrm{eq,core}}$'))
        hcross_nuc.set(name='nuc_hcross_eq',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{eq,core}}$'))
        hplus_inn.set(name='inn_hplus_eq',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{+,\mathrm{eq,conv}}$'))
        hcross_inn.set(name='inn_hcross_eq',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{eq,conv}}$'))
        hplus_out.set(name='out_hcross_eq',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{+,\mathrm{eq,out}}$'))
        hcross_out.set(name='out_hcross_eq',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{eq,out}}$'))
        [pQ.set(name=f'r{i}_hplus_eq',
              cmap='seismic',
              limits=[-70 / D.value, 70 / D.value],
              label=merge_strings(add_lb, r'$h_{+,\mathrm{eq,r$',
                                  f'$_{i}$', r'$}}$'))
              for i, pQ in enumerate(partialQ_plus)]
        [pQ.set(name=f'r{i}_hcross_eq',
              cmap='seismic',
              limits=[-70 / D.value, 70 / D.value],
              label=merge_strings(add_lb, r'$h_{\times,\mathrm{eq,r$',
                                  f'$_{i}$', r'$}}$'))
              for i, pQ in enumerate(partialQ_cross)]
    else:
        hplus_radial.set(name='hpuls_radial',
                         label=merge_strings(add_lb,
                                             r'$h_\mathrm{+}(r)$'),
                         cmap='seismic',
                         limits=[-3 / D.value, 3 / D.value])        
        hcross_radial.set(name='hcross_radial',
                         label=merge_strings(add_lb,
                                             r'$h_\mathrm{\times}(r)$'),
                         cmap='seismic',
                         limits=[-3 / D.value, 3 / D.value])
        hplus_tot.set(name='tot_hplus',
                      cmap='seismic',
                      limits=[-70 / D.value, 70 / D.value],
                      label=merge_strings(add_lb, r'$h_{+,\mathrm{tot}}$'))
        hcross_tot.set(name='tot_hcross',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{tot}}$'))
        hplus_nuc.set(name='nuc_hplus',
                      cmap='seismic',
                      limits=[-70 / D.value, 70 / D.value],
                      label=merge_strings(add_lb, r'$h_{+,\mathrm{core}}$'))
        hcross_nuc.set(name='nuc_hcross',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{core}}$'))
        hplus_inn.set(name='inn_hplus',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{+,\mathrm{conv}}$'))
        hcross_inn.set(name='inn_hcross',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{conv}}$'))
        hplus_out.set(name='out_hcross',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{+,\mathrm{out}}$'))
        hcross_out.set(name='out_hcross',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{out}}$'))
        [pQ.set(name=f'r{i}_hplus',
              cmap='seismic',
              limits=[-70 / D.value, 70 / D.value],
              label=merge_strings(add_lb, r'$h_{+,\mathrm{r$',
                                  f'$_{i}$', r'$}}$'))
              for i, pQ in enumerate(partialQ_plus)]
        [pQ.set(name=f'r{i}_hcross',
              cmap='seismic',
              limits=[-70 / D.value, 70 / D.value],
              label=merge_strings(add_lb, r'$h_{\times,\mathrm{r$',
                                  f'$_{i}$', r'$}}$'))
              for i, pQ in enumerate(partialQ_cross)]
    
    T, R = np.meshgrid(time, radius)
    return aeseries(data=hplus_radial, time=T.copy(), radius=R.copy()), \
           return_3D_strains(time, hplus_tot, hplus_nuc, hplus_inn, hplus_out,
                             partialQ_plus, radii), \
           aeseries(data=hcross_radial, time=T.copy(), radius=R.copy()), \
           return_3D_strains(time, hcross_tot, hcross_nuc, hcross_inn,
                             hcross_out, partialQ_cross, radii)
