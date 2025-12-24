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
from typing import Literal


## ---------------------------------------------------------------------
## GW strain
## ---------------------------------------------------------------------

def GW_strain(sim_dim, column_change, data, index, ref, distance,
              return_components=False):
    assert sim_dim in [1, 2, 3], "Simulation MUST be 1, 2 or 3D."
    if distance:
        if not isinstance(distance, aerray):
            distance = distance * u.cm
    if sim_dim == 1:
        return GW_strain_1D(data)
    elif sim_dim == 2:
        return correct_zero(2, GW_strain_2D(data[::ref, :], distance), index)
    else:
        GWs = GW_strain_3D(data, distance, return_components, ref)
        if return_components:
            return GWs
        if column_change is not None:    
            GWs[0].data[:column_change] = GW_strain_2D(data, distance).data[:column_change]
            for hh in range(1, len(GWs)):
                GWs[hh].data.value[:column_change-1] = np.zeros(column_change-1)
            GWs = match_remap(remove_3D_spikes(GWs, column_change, distance),
                              column_change)
        return correct_zero(3, refine_3D(GWs, ref), index)

def GW_strain_1D(data):
    print("No GW for you :'(")
    return None

def GW_strain_2D(data, distance):
    """
        GWs amplitudes calculated as the first partial time derivative
        of NE_220 and not with the second partial time derivative of 
        ME_220. Moreover, we only consider matter contribution to the 
        amplitude, being the main contribution to it.
    """
    const = -0.125 *  np.sqrt(15/np.pi)
    GWs = aeseries(
        aerray(const * IDL_derivative(data[:,2], data[:,5]), u.cm, 'h+eq',
               r'$\mathcal{D}h_{+,eq}$', 'Spectral_r', [-150, 150]),
        time=aerray(data[:, 2], u.s, 'time', r'$t$', None, [0, data[-1, 2]])
    )
    if distance:
        lb, cm, lm, nm = GWs.data.label, GWs.data.cmap, \
            GWs.data.limits, GWs.data.name
        lm = [lm[0] / distance.to(GWs.data.unit).value,
              lm[1] / distance.to(GWs.data.unit).value]
        GWs /= distance
        GWs.data.set(label=lb.replace(r'\mathcal{D}', r''), name=nm, cmap=cm,
                     limits=lm)
    return GWs

def GW_strain_3D(data, distance, return_components=False, ref=None):
    if return_components and ref:
        data = data[::ref, :]
    time = aerray(data[:, 2], u.s, 'time', r'$t$', None, [0, data[-1, 2]])
    const = 1 #8 * np.sqrt(np.pi / 15)
    if return_components:
        return time, [IDL_derivative(time, (data[:, i] * u.cm * u.s)) 
                      for i in [9, 13, 17, 10, 11, 14]]
    hD_pl_p = 2. * ( data[:,9] - data[:,13] )
    hD_pl_e = 2. * ( data[:,17] - data[:,13] )
    hD_cr_p = 2. * ( data[:,10] + data[:,12] )
    hD_cr_e = 2. * ( - data[:,14] - data[:,16] )

    hD_pl_e = aeseries(
        aerray(const * IDL_derivative( data[:,2], hD_pl_e ), u.cm, 'h+eq',
               r'$\mathcal{D}h_{+,eq}$', 'Spectral_r', [-150, 150]),
        time=time)
    hD_pl_p = aeseries(
        aerray(const * IDL_derivative( data[:,2], hD_pl_p ), u.cm, 'h+pol',
               r'$\mathcal{D}h_{+,pol}$', 'Spectral_r', [-150, 150]),
        time=time.copy())
    hD_cr_e = aeseries(
        aerray(const * IDL_derivative( data[:,2], hD_cr_e ), u.cm, 'hxeq',
               r'$\mathcal{D}h_{\times,eq}$', 'Spectral_r', [-150, 150]),
        time=time.copy())
    hD_cr_p = aeseries(
        aerray(const * IDL_derivative( data[:,2], hD_cr_p ), u.cm, 'hxpol',
               r'$\mathcal{D}h_{\times,pol}$', 'Spectral_r', [-150, 150]),
        time=time.copy())
    GWs = [hD_pl_e, hD_pl_p, hD_cr_e, hD_cr_p]
    if distance:
        for hh in range(len(GWs)):
            lb, cm, lm, nm = GWs[hh].data.label, GWs[hh].data.cmap, \
                GWs[hh].data.limits, GWs[hh].data.name
            lm = [lm[0] / distance.to(GWs[hh].data.unit).value,
              lm[1] / distance.to(GWs[hh].data.unit).value]
            GWs[hh] /= distance
            GWs[hh].data.set(label=lb.replace(r'\mathcal{D}', r''), name=nm, cmap=cm,
                        limits=lm)
    return GWs

def refine_3D(GWs, ref):
    if ref == 1:
        return GWs
    time = aerray(GWs[0].time.value[::ref], GWs[0].time.unit, GWs[0].time.name,
                  GWs[0].time.label, GWs[0].time.cmap, GWs[0].time.limits,
                  GWs[0].time.log)
    for hh in range(len(GWs)):
        h = aerray(GWs[hh].data.value[::ref], GWs[hh].data.unit,
                   GWs[hh].data.name, GWs[hh].data.label, GWs[hh].data.cmap,
                   GWs[hh].data.limits, GWs[hh].data.log)
        GWs[hh] = aeseries(h, time = time.copy())
    return GWs

def correct_zero(sim_dim, GWs, index):
    if sim_dim == 1:
        pass
    else:
        if index is None:
            return GWs
        if sim_dim == 2:
            GWs -= GWs.data[:index].mean()
        else:
            for i in range(len(GWs)):
                GWs[i] -= GWs[i].data[:index].mean()
        return GWs

def remove_3D_spikes(GWs, index, distance):
    """
    Around the 2D-3D matching there are spikes at least in some cases,
    let's try to remove them.
    """
    if distance:
        check_GWs = [hh.data * distance for hh in GWs]
    else:
        check_GWs =  [hh.data for hh in GWs]
    for h in range(len(GWs)):
        for i in range(index-50, index+50):
            if check_GWs[h][i] >= 1e3 or check_GWs[h][i] <= -1e3:
                GWs[h].data[i] = GWs[h].data[i-1]            
    return GWs

def match_remap(GWs, index):
    for h in range(len(GWs)):
        diff = np.mean(GWs[h].data[index-8:index]) - \
            np.mean(GWs[h].data[index:index+8])
        GWs[h].data[index:] += diff
    return GWs
        
## ---------------------------------------------------------------------
## GW Energy
## ---------------------------------------------------------------------

def GWs_energy(GWs, sim_dim):
    assert sim_dim in [1, 2, 3], "Simulation MUST be 1, 2 or 3D."
    if sim_dim == 1:
        return GWs_energy_1D(GWs)
    elif sim_dim == 2:
        return GWs_energy_2D(GWs)
    else:
        return GWs_energy_3D(GWs)

def GWs_energy_1D(GWs):
    print("No GW for you :'(")
    return None

def GWs_energy_2D(GWs):
    """
    Calculates the energy of the GWs in 2D
    """
    const = c.c ** 3 / (c.G * 32 * np.pi) / (-0.125 *  np.sqrt(15/np.pi)) ** 2
    ene = aeseries(
        IDL_derivative(GWs.time, GWs.data) ** 2 * const,
        time = GWs.time
    )
    ene.data.set(name='EGWs', label=r'$L_\mathrm{GW}$', cmap='Spectral',
                 limits=[1e45, 1e48], log=True)
    return ene.to(u.erg / u.s)

def GWs_energy_3D(GWs):
    """
    Calculates the energy of the GWs in 3D
    Formulae from `10.1051/0004-6361:20078577` eq: (13)
    """
    const = 2 * c.c ** 3 / 15 / c.G
    time, series = GWs
    for i, h in enumerate(series):
        series[i] = IDL_derivative(time, h)
    hxx, hyy, hzz, hxy, hxz, hyz = series
    ene = hxx ** 2 + hyy ** 2 + hzz ** 2 - \
        (hxx * hyy + hxx * hzz + hyy * hzz) + \
            3 * (hxy **2 + hxz ** 2 + hyz ** 2)
    ene *= const
    ene.set(name='EGWs', label=r'$L_\mathrm{GW}$', cmap='Spectral',
            limits=[1e45, 1e48], log=True)
    ene = aeseries(np.where(ene <= (1e51*u.erg/u.s), ene, 0), time = time)
    return ene.to(u.erg / u.s)

def GWs_energy_per_frequency(GWs, sim_dim, time_range=None, windowing='hanning'):
    assert sim_dim in [1, 2, 3], "Simulation MUST be 1, 2 or 3D."
    if sim_dim == 1:
        return GWs_energy_per_frequency_1D(GWs)
    elif sim_dim == 2:
        return GWs_energy_per_frequency_2D(GWs, time_range, windowing)
    else:
        return GWs_energy_per_frequency_3D(GWs, time_range, windowing)

def GWs_energy_per_frequency_1D(GWs):
    print("No GW for you :'(")
    return None

def GWs_energy_per_frequency_2D(GWs, time_range=None, windowing='hanning'):
    const1 = 1/8 * np.sqrt(15 / np.pi)
    const2 = c.c ** 3 / (16 * np.pi * c.G)
    nm, lb, cm, lg = GWs.data.name, GWs.data.label, GWs.data.cmap, GWs.data.log  
    GWs.data /= const1
    GWs.data.set(name=nm, label=lb, cmap=cm, log=lg) 
    spectro = GWs.rfft(norm='forward', time_range=time_range,
                       windowing=windowing)
    dE_df = const2 * (2 * np.pi * spectro.frequency) ** 2 * \
        np.abs(spectro.data ** 2)
    dE_df.set(name='dE_df_h', label=r'$\frac{\mathrm{d}E}{\mathrm{d}f}$',
              log=True, cmap=None)
    return aeseries(
        dE_df,
        frequency=spectro.frequency.copy()
    )

def GWs_energy_per_frequency_3D(GWs, time_range=None, windowing='hanning'):
    """
    Computed for the two observers, polar and equatorial following the
    formulation in Kuroda et al 2014. Eq (45). `1304.4372`
    """
    const = np.pi * c.c ** 3 / (4 * c.G)
    names = ['dE_df_heq', 'dE_df_hpol']
    labels = [r'$\frac{\mathrm{d}E}{\mathrm{d}f}_\mathrm{eq}$',
              r'$\frac{\mathrm{d}E}{\mathrm{d}f}_\mathrm{pol}$']
    spectro_eq_pl = GWs[0].rfft(norm='forward', time_range=time_range,
                                windowing=windowing)
    spectro_eq_cr = GWs[2].rfft(norm='forward', time_range=time_range,
                                windowing=windowing)
    dE_df_eq = const * spectro_eq_pl.frequency ** 2 * (
        np.abs(spectro_eq_pl.data ** 2) + np.abs(spectro_eq_cr.data ** 2)
    )
    dE_df_eq.set(name=names[0], label=labels[0], log=True)
    spectro_pol_pl = GWs[1].rfft(norm='forward', time_range=time_range,
                                windowing=windowing)
    spectro_pol_cr = GWs[3].rfft(norm='forward', time_range=time_range,
                                windowing=windowing)
    dE_df_pol = const * spectro_pol_pl.frequency ** 2 * (
        np.abs(spectro_pol_pl.data ** 2) + np.abs(spectro_pol_cr.data ** 2)
    )
    dE_df_pol.set(name=names[1], label=labels[1], log=True)

    dE_df = [
        aeseries(dE_df_eq, frequency=spectro_eq_cr.frequency.copy()),
        aeseries(dE_df_pol, frequency=spectro_pol_cr.frequency.copy())
    ]
    return dE_df

## ---------------------------------------------------------------------
## GW characteristic strain
## ---------------------------------------------------------------------

"""
Taken from:
Flanagan+98, https://journals.aps.org/prd/pdf/10.1103/PhysRevD.57.4535
"""

def characteristic_strain(GWs, sim_dim, time_range=None, windowing='hanning',
                          distance=(10 * u.kpc), divide_by_frequency=True):
    assert sim_dim in [1, 2, 3], "Simulation MUST be 1, 2 or 3D."
    if sim_dim == 1:
        return characteristic_strain_1D(GWs)
    elif sim_dim == 2:
        return characteristic_strain_2D(GWs, time_range, windowing, distance,
                                        divide_by_frequency)
    else:
        return characteristic_strain_3D(GWs, time_range, windowing, distance,
                                        divide_by_frequency)
    
def characteristic_strain_1D(GWs):
    print("No GW for you :'(")
    return None

def characteristic_strain_2D(GWs, time_range, windowing, distance,
                             divide_by_frequency):
    const = 1 / (distance.to(GWs.data.unit) * c.c * np.pi) * \
        np.sqrt(2 * c.G / c.c)
    dE_df = GWs_energy_per_frequency_2D(GWs, time_range, windowing)
    hchar = const * np.sqrt(dE_df.data)
    to_unit = u.dimensionless_unscaled
    if divide_by_frequency:
        hchar /= np.sqrt(dE_df.frequency)
        hchar.set(name='hchar', label=r'$h_\mathrm{char,+,eq}/\sqrt{f}$', log=True,
              limits=[np.nanmin(hchar), np.nanmax(hchar)])
        to_unit = u.Hz ** -0.5
    else:
        hchar.set(name='hchar', label=r'$h_\mathrm{char,+,eq}$', log=True,
              limits=[np.nanmin(hchar), np.nanmax(hchar)])
    dE_df.frequency.set(limits=[1, 4000], log=True)
    return aeseries(hchar.to(to_unit), frequency=dE_df.frequency.copy())

def characteristic_strain_3D(GWs, time_range, windowing, distance,
                             divide_by_frequency):
    """
    Computes the GWs characteristic frequency from Kuroda et al 2014,
    Eq (44)
    """
    const = 2 / np.pi ** 2 * c.G / c.c ** 3 / distance ** 2
    dE_df = GWs_energy_per_frequency_3D(GWs, time_range, windowing)
    names = ['hchar_eq', 'hchar_pol']
    labels = [r'$h_\mathrm{char,eq}$',
              r'$h_\mathrm{char,pol}$']
    to_unit = u.dimensionless_unscaled
    if divide_by_frequency:
        labels = [merge_strings(lb, r'$/\sqrt{f}$') for lb in labels]
        to_unit = u.Hz ** -0.5
    hchar = []
    for i, dedf in enumerate(dE_df):
        hhchar = np.sqrt(dedf.data * const)
        if divide_by_frequency:
            hhchar /= np.sqrt(dedf.frequency)
        hhchar.set(name=names[i], label=labels[i], log=True,
              limits=[np.nanmin(hhchar), np.nanmax(hhchar)])
        dedf.frequency.set(limits=[1, 4000], log=True)
        hchar.append(aeseries(hhchar.to(to_unit),
                              frequency=dedf.frequency.copy()))
    return hchar

## ---------------------------------------------------------------------
## GW spectrogram
## ---------------------------------------------------------------------

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

def GWs_spectrogram(sim_dim, GWs, window_size, scale_to, **kwargs):
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
## GWs peaks
## ---------------------------------------------------------------------

## PEAKS

def GWs_peak_indices(GWs, peak, interval, min_time, max_time):
    """
    Function that finds the coordinates of the minimum and maximum peak
    of the equatorial + GW strain as well as the coordinates of the 
    points before and after the oscillation. Namely the latter are the 
    first intersection point with the x axis after the peak and the
    third one before it.
    Parameters:
        peak: which peak to find, can be the bounce peak, highest in
        an interval
        interval: interval (ms) in which the peak has to be found,
        if only one value is provided, that would be used as the 
        right hand side
    Returns:
        list containing left index, peak index, right index
    """
    zeros = np.where(GWs[:-1, 1] * GWs[1:, 1] < 0 )[0] + 1
    if peak == 'bounce':
        ## FIND the bounce time
        bounce_index = np.argmax(GWs[:, 0] >= 0)
        ## FIND the min in the 1.5 ms after the bounce
        index_after_bounce = np.argmax(GWs[:,0] \
                                        >= u.convert_to_s(min_time)) + 1
        x_min = np.argmin(GWs[bounce_index:index_after_bounce, 1]) \
            + bounce_index
        ##FIND MAX AFTER the min in 1.5 ms
        index2_after_bounce = index_after_bounce = \
            np.argmax(GWs[:,0] >= GWs[x_min, 0] + u.convert_to_s(max_time)) + 1
            
        x_max = np.argmax(GWs[x_min:index2_after_bounce, 1]) + x_min
    elif peak == 'highest':
        ## CUT the GWS
        if interval[0] is not None:
            start_index = np.argmax(GWs[:, 0] >= u.convert_to_s(interval[0]))
            GWs = GWs[start_index:, :]
        else:
            start_index = None
        if interval[1] is not None:
            GWs = GWs[:np.argmax(GWs[:, 0] >= u.convert_to_s(interval[1])), :]
        ## FIND the peak
        x_max = np.argmax(np.abs(GWs[:, 1]))
        if GWs[x_max, 1] < 0:
            GWs[:, 1] = -GWs[:, 1]
        ## FIND the min
        min_index = np.argmax(GWs[:, 0] >= (GWs[x_max, 0] - \
                                            u.convert_to_s(max_time)))
        x_min = np.argmin(GWs[ min_index:np.argmax(GWs[:, 0] >= \
                    ( GWs[x_max, 0] + u.convert_to_s(max_time) )), 1 ]) + \
                    min_index
        if start_index is not None:
            x_min += start_index
            x_max += start_index
    else:
        raise ValueError("Peak must be either 'bounce' or 'highest'.")
    ## Find the beginning and end of the peak
    if x_max > x_min:
        zeros_end_index = np.argmax(zeros>x_max)
        end_index = zeros[zeros_end_index]
        start_index = zeros[np.argmax(zeros>x_min) - 4]
    elif x_max < x_min:
        zeros_end_index = np.argmax(zeros>x_min)
        end_index = zeros[zeros_end_index]
        start_index = zeros[np.argmax(zeros>x_max) - 4]
    if start_index > x_max:
        start_index = 0
    return start_index, x_min, x_max, end_index

def GWs_max_peak(GWs, peak, interval, min_time, max_time):
    """
    Function that finds the maximum value of the equatorial + 
    polarization of the GW strain as well as the time at which it
    occurs.
    """
    indices = GWs_peak_indices(GWs, peak, interval, min_time, max_time)
    return GWs[indices[1], 0], GWs[indices[1], 1]

## FREQUENCIES

def GWs_fourier_transform(GWs, indices):
    """
    This function applies FFT to a small portion of the GW strain to
    find the domiunant frequency of a specific oscillation
    Returns
        positive frequency range
        tilde{h} * sqrt{freq}
    """
    dt = np.abs(GWs[1, 0] - GWs[0, 0])
    ## Cut the GWs signal
    strain = np.zeros(10000)
    ## Pad the strain with zeros to increase the resolution
    if (GWs[indices[0]:indices[-1], 1]).size < 11000:
        strain[10000 - (GWs[indices[0]:indices[-1], 1]).size:] = \
            GWs[indices[0]:indices[-1], 1]
    else:
        strain = GWs[indices[0]:indices[-1], 1]
    ## Find the frequencies
    freq = fftfreq(strain.size, dt)[:strain.size//2]
    dft = np.abs(fft(strain))[:strain.size//2] * np.sqrt(freq)
    return freq, dft

def GWs_frequency_peak_indices(frequency, htilde):
    """
    Function: finds the indices of the first and second frequency
    peak on a fourier transformed GW strain
    Return:
        list containing: first peak index, second peak index 
    """
    dhtilde_df = IDL_derivative(frequency, htilde)
    sign = dhtilde_df[1:] * dhtilde_df[:-1]
    ## Find extrema points
    extr = np.where(sign < 0)[0]
    ## Find the two peaks
    first_peak_index = 0
    second_peak_index = 0
    for index in extr:
        if htilde[index] > htilde[first_peak_index]:
            second_peak_index = first_peak_index
            first_peak_index = index
        elif htilde[index] > htilde[second_peak_index]:
            second_peak_index = index
    return [first_peak_index, second_peak_index]

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

def compute_SNR(simulation, detector, comp, distance=(u.kpc *10),
                time_range=None):
    """
    Computes the SNR as in Moore 2015. Please be careful, for this to be
    done correctly, the signal should be zero padded and the frequency
    range on which to integrate over should be carefully selected.
    However, the method proposed here is automatized and does not
    deviate from the true value by that much.
    """
    def interpolate_ASD(ASD, freq):
        new_ASD =  np.interp(freq, ASD.frequency, ASD.data)
        new_ASD = new_ASD.data
        new_ASD = aerray(new_ASD, (u.Hz**-0.5), '', r'$ASD_\mathrm{' + ASD.data.name + '}$')
        ASD = aeseries(
            new_ASD,
            frequency = freq.copy()
        )
        return ASD
    ## Compute the characteristic strain wich is dimensionless
    hchar = simulation.hchar(comp=comp, distance=distance,
                             time_range=time_range, divide_by_frequency=False)
    ## Load the selected ASD
    ASD = simulation.ASD(detector)
    ## Now we need the ASD and the characteristic strain to have the same
    ## frequency range. So we cut one or the other frequency series.
    if ASD.frequency[-1] > hchar.frequency[-1]:
        iend = np.argmax(ASD.frequency >= hchar.frequency[-1])
        ASD = ASD[:iend]
    else:
        iend = np.argmax(hchar.frequency >= ASD.frequency[-1])
        hchar[:iend]
    if ASD.frequency[0] < hchar.frequency[0]:
        istart = np.argmax(ASD.frequency >= hchar.frequency[0])
        ASD = ASD[istart:]
    else:
        istart = np.argmax(hchar.frequency >= ASD.frequency[0])
        hchar = hchar[istart:]
    
    ## Now we have to interpolate the ASD (since is smoother) to the
    ## frequency range of the characteristic strain
    ASD = interpolate_ASD(ASD, hchar.frequency)
    ## Now we proceed with the actual computation of the SNR
    ## We assume constant spacing in the frequency domain
    df = np.mean(np.diff(hchar.frequency))
    SNR = df * hchar.data ** 2 / (ASD.data ** 2 * ASD.frequency ** 2)
    ##Consistency check
    if SNR.unit != u.dimensionless_unscaled:
        raise ValueError("Units mismatch")
    return np.sqrt(np.sum(SNR.value)), detector


