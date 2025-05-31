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


## ---------------------------------------------------------------------
## GW strain
## ---------------------------------------------------------------------

def GW_strain(sim_dim, column_change, data, index, distance):
    assert sim_dim in [1, 2, 3], "Simulation MUST be 1, 2 or 3D."
    if distance:
        if not isinstance(distance, aerray):
            distance = distance * u.cm
    if sim_dim == 1:
        return GW_strain_1D(data)
    elif sim_dim == 2:
        return correct_zero(2, GW_strain_2D(data, distance), index)
    else:
        GWs = GW_strain_3D(data, distance)
        if column_change is not None:
            GWs[0].data[:column_change] = GW_strain_2D(data, distance).data[:column_change]
            for hh in range(1, len(GWs)):
                GWs[hh].data.value[:column_change] = np.zeros(column_change)
        return correct_zero(3, GWs, index)

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

def GW_strain_3D(data, distance):
    time = aerray(data[:, 2], u.s, 'time', r'$t$', None, [0, data[-1, 2]])
    const = 1 #8 * np.sqrt(np.pi / 15)
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
               r'$\mathcal{D}h_{x,eq}$', 'Spectral_r', [-150, 150]),
        time=time.copy())
    hD_cr_p = aeseries(
        aerray(const * IDL_derivative( data[:,2], hD_cr_p ), u.cm, 'hxpol',
               r'$\mathcal{D}h_{x,pol}$', 'Spectral_r', [-150, 150]),
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
    ene.data.set(name='EGWs', label=r'$E_\mathrm{GW}$', cmap='Spectral',
                 limits=[1e45, 1e48], log=True)
    return ene.to(u.erg / u.s)

def GWs_energy_3D(GWs):
    """
    Calculates the energy of the GWs in 3D
    """
    const = 8 / 15 * c.c ** 3 / (16 * np.pi * c.G)
    ene = aeseries(IDL_derivative(GWs[0].time, GWs[0].data) ** 2,
                   time=GWs[0].time)
    for h in GWs[1:]:
        ene += IDL_derivative(h.time, h.data) ** 2
    ene *= const
    ene.data.set(name='EGWs', label=r'$E_\mathrm{GW}$', cmap='Spectral',
            limits=[1e45, 1e48], log=True)
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
    dE_df = []
    const = c.c ** 3 / (16 * np.pi * c.G)
    names = ['dE_df_h+eq', 'dE_df_h+pol', 'dE_df_hxeq', 'dE_df_hxpol']
    labels = [r'$\frac{\mathrm{d}E}{\mathrm{d}f}_\mathrm{+,eq}$',
              r'$\frac{\mathrm{d}E}{\mathrm{d}f}_\mathrm{+,pol}$',
              r'$\frac{\mathrm{d}E}{\mathrm{d}f}_\mathrm{+,pol}$',
              r'$\frac{\mathrm{d}E}{\mathrm{d}f}_\mathrm{x,pol}$']
    for i, GW in enumerate(GWs):
        spectro = GWs.rfft(norm='forward', time_range=time_range,
                       windowing=windowing)
        dedf = const * (2 * np.pi * spectro.frequency) ** 2 * \
               np.abs(spectro.data ** 2)
        dedf.set(name=names[i], label=labels[i], log=True)
        dE_df.append(aeseries(dedf, frequency=spectro.frequency.copy()))
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
    if divide_by_frequency:
        hchar /= np.sqrt(dE_df.frequency)
        hchar.set(name='hchar', label=r'$h_\mathrm{char,+,eq}/\sqrt{f}$', log=True,
              limits=[np.nanmin(hchar), np.nanmax(hchar)])
    else:
        hchar.set(name='hchar', label=r'$h_\mathrm{char,+,eq}$', log=True,
              limits=[np.nanmin(hchar), np.nanmax(hchar)])
    dE_df.frequency.set(limits=[1, 4000], log=True)
    return aeseries(hchar.to((u.Hz ** -0.5)), frequency=dE_df.frequency.copy())

def characteristic_strain_3D(GWs, time_range, windowing, distance,
                             divide_by_frequency):
    const = 1 / (distance.to(GWs[0].data.unit) * c.c * np.pi) * np.sqrt(2 * c.G / c.c)
    dE_df = GWs_energy_per_frequency_3D(GWs, time_range, windowing)
    names = ['hchar+eq', 'hchar+pol', 'hcharxeq', 'hcharxpol']
    labels = [r'$h_\mathrm{char,+,eq}$',
              r'$h_\mathrm{char,+,pol}$',
              r'$h_\mathrm{char,x,eq}$',
              r'$h_\mathrm{char,x,pol}$']
    if divide_by_frequency:
        labels = [merge_strings(lb, r'$\sqrt{f}$') for lb in labels]
    hchar = []
    for i, dedf in enumerate(dE_df):
        hhchar = const * np.sqrt(dedf.data)
        if divide_by_frequency:
            hhchar /= np.sqrt(dedf.frequency)
        hhchar.set(name=names[i], label=labels[i], log=True,
              limits=[np.nanmin(hhchar), np.nanmax(hhchar)])
        dedf.frequency.set(limits=[1, 4000], log=True)
        hchar.append(aeseries(hhchar.to((u.Hz ** -0.5)),
                              frequency=dedf.frequency.copy()))
    return hchar
    
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
                save_checkpoints=True):
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
    if simulation.dim == 1:
        print("No GWs for you :'(")
        return None
    elif simulation.dim == 2:
        return NE220_2D_timeseries(simulation, save_checkpoints, D)
    elif simulation.dim == 3:
        return Qdot_timeseries(simulation, save_checkpoints, D, THETA, PHI)

## 2D

def NE220_2D_timeseries(simulation, save_checkpoints, D):
    """
    Calculates the NE220 from density and velocities for every timestep
    of a 2D simulation. It also calculates the full, nucleus, convection
    and outer core contributions to the strain.
    """
    if check_existence(simulation, 'NE220.h5'):
        time, NE220, full_NE220, nuc_NE220, conv_NE220, outer_NE220, processed_hdf = \
            read_NE220(simulation)
        if processed_hdf is None:
            save_hdf(os.path.join(simulation.storage_path, 'NE220.h5'),
                     ['time', 'NE220', 'full_NE220', 'nucleus_NE220',
                      'convection_NE220', 'outer_NE220', 'processed'],
                     [time, NE220, full_NE220, nuc_NE220, conv_NE220,
                      outer_NE220, simulation.hdf_file_list[:len(time)]])
            time, NE220, full_NE220, nuc_NE220, conv_NE220, outer_NE220, processed_hdf = \
            read_NE220(simulation)
        if processed_hdf[-1].decode("utf-8") == simulation.hdf_file_list[-1]:
            return calculate_strain_2D(D, time,
                                       simulation.cell.radius(simulation.ghost),
                                       NE220, full_NE220, nuc_NE220,
                                       conv_NE220, outer_NE220)
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
        fNE220, ffull, finner, fnuc, fouter = NE220_2D(simulation,
            file, dV, dOmega, ctheta, inner_rad[..., findex], igcells,
            nuc_rad[..., findex], ngcells)
        try:
            time = np.concatenate((time, simulation.time(file)))
            NE220 = np.concatenate((NE220, fNE220[..., None]), axis=-1)
            full_NE220 = np.concatenate((full_NE220, ffull))
            nuc_NE220 = np.concatenate((nuc_NE220, fnuc))
            conv_NE220 = np.concatenate((conv_NE220, finner))
            outer_NE220 = np.concatenate((outer_NE220, fouter))
        except Exception as e:
            print(e)
            time = simulation.time(file)
            NE220 = fNE220[..., None]
            full_NE220 = ffull
            nuc_NE220 = fnuc
            conv_NE220 = finner
            outer_NE220 = fouter
        processed_hdf.append(file)
        if save_checkpoints and check_index == checkpoint:
            save_hdf(os.path.join(simulation.storage_path, 'NE220.h5'),
                     ['time', 'NE220', 'full_NE220', 'nucleus_NE220',
                      'convection_NE220', 'outer_NE220', 'processed'],
                     [time, NE220, full_NE220, nuc_NE220, conv_NE220,
                      outer_NE220, processed_hdf])
            check_index = 0
        check_index += 1
        progress_index += 1
        findex += 1
        progressBar(progress_index, len(simulation.hdf_file_list))
    
    print("Computations done, saving...")
    save_hdf(os.path.join(simulation.storage_path, 'NE220.h5'),
                     ['time', 'NE220', 'full_NE220', 'nucleus_NE220',
                      'convection_NE220', 'outer_NE220', 'processed'],
                     [time, NE220, full_NE220, nuc_NE220, conv_NE220,
                      outer_NE220, processed_hdf])
    return calculate_strain_2D(D, time, simulation.cell.radius(simulation.ghost),
                               NE220, full_NE220, nuc_NE220, conv_NE220,
                               outer_NE220)

def Zha_correction_2D(dOmega, ctheta, r1_ind, r2_ind, r, rho, vr):
    """
    Computes the Zha correction into the strains computed on spherica shells
    """
    P2 = 0.5 * (3 * ctheta ** 2 - 1) * dOmega[:, None]
    
    prod = r[None, :] ** 4 * rho * vr * P2
    if r1_ind is None:
        r1 = 0 * prod.unit
    else:
        r1 = prod[np.arange(prod.shape[0]), r1_ind]
        r1 = r1.sum()
    if r2_ind is None:
        r2 = 0 * prod.unit
    else:
        r2 = prod[np.arange(prod.shape[0]), r2_ind]
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
    nuc_corr = Zha_correction_2D(dOmega, ctheta, None, r_nuc_ind, radius, rho, vr)
    inn_corr = Zha_correction_2D(dOmega, ctheta, r_nuc_ind+1, r_inner_ind, radius, rho, vr)
    out_corr = Zha_correction_2D(dOmega, ctheta, r_inner_ind+1, r_outer_ind, radius, rho, vr)
    return np.sum(NE220, axis=0), np.sum(NE220), np.sum(NE220 * mask_inner) + nuc_corr, \
        np.sum(NE220 * mask_nuc)+inn_corr, np.sum(NE220 * mask_outer)+out_corr
           
def read_NE220(simulation):
    """
    Reads the NE220 from a checkpoint file.
    """
    data = h5py.File(os.path.join(simulation.storage_path, 'NE220.h5'), 'r')
    time = data['time'][...] * u.s
    NE220 = data['NE220'][...] * u.cm ** 2 * u.g / u.s
    full_NE220 = data['full_NE220'][...] * u.cm ** 2 * u.g / u.s
    nuc_NE220 = data['nucleus_NE220'][...] * u.cm ** 2 * u.g / u.s
    conv_NE220 = data['convection_NE220'][...] * u.cm ** 2 * u.g / u.s
    outer_NE220 = data['outer_NE220'][...] * u.cm ** 2 * u.g / u.s
    if 'processed' in data.keys():
        processed = data['processed'][...]
    else:
        processed = None
    data.close()
    return time, NE220, full_NE220, nuc_NE220, conv_NE220, outer_NE220, processed

def calculate_strain_2D(D, time, radius, NE220, full_NE220, nuc_NE220,
                        conv_NE220, outer_NE220):
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
    NE220 = const * IDL_derivative(time, NE220)
    NE220.set(name='AE220', label=merge_strings(add_lb, r'$A^{E2}_{20}(r)$'),
              cmap='seismic', limits=[-3 / D.value, 3 / D.value])
    full_NE220 = const * IDL_derivative(time, full_NE220)
    full_NE220.set(name='full_NE220', label=merge_strings(add_lb, r'$h_{+,eq}$'),
                   cmap='seismic', limits=[-70 / D.value, 70 / D.value])
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
    return aeseries(NE220, time=T, radius=R),\
            create_series(time, full_NE220, nuc_NE220, conv_NE220, outer_NE220)

## 3D

def Qdot_timeseries(simulation, save_checkpoints, D, THETA, PHI):
    """
    Calculates the NE220 from density and velocities for every timestep
    of a 2D simulation. It also calculates the full, nucleus, convection
    and outer core contributions to the strain.
    """
    if check_existence(simulation, 'Qdot.h5'):
        time, Qdot_radial, Qdot_total, Qdot_inner, Qdot_nucleus, Qdot_outer, processed_hdf = \
            read_Qdot(simulation)
        if processed_hdf is None:
            save_hdf(os.path.join(simulation.storage_path, 'Qdot.h5'),
                     ['time', 'Qdot_total', 'Qdot_inner', 'Qdot_nucleus',
                      'Qdot_outer', 'Qdot_radial', 'processed'],
                     [time, Qdot_total, Qdot_inner, Qdot_nucleus, Qdot_outer,
                      Qdot_radial, simulation.hdf_file_list[:len(time)]])
            time, Qdot_radial, Qdot_total, Qdot_inner, Qdot_nucleus, Qdot_outer, processed_hdf = \
            read_Qdot(simulation)
        if processed_hdf[-1].decode("utf-8") == simulation.hdf_file_list[-1]:
            return calculate_strain_3D(D, THETA, PHI, time,
                                       simulation.cell.radius(simulation.ghost),
                                       Qdot_radial,
                                       Qdot_total, Qdot_inner, Qdot_nucleus,
                                       Qdot_outer)
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
    inner_rad, igcells = simulation.innercore_radius(rad='full')
    nuc_rad, ngcells = simulation.PNS_nucleus_radius(rad='full')

    inner_rad = inner_rad.data
    nuc_rad = nuc_rad.data
        
    grad = spherical_harmonics_gradient(
                                    simulation.cell.radius(simulation.ghost),
                                    simulation.cell.theta(simulation.ghost),
                                    simulation.cell.phi(simulation.ghost))
    
    for file in simulation.hdf_file_list[start_point:]:
        Qtot, Qinner, Qnuc, Qouter, Qradial = calculate_Qdot(simulation, 
                            grad, file, dV,
                            inner_rad[..., findex], igcells,
                            nuc_rad[..., findex], ngcells)
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
        except:
            time = simulation.time(file)
            Qdot_radial = Qradial[..., None]
            Qdot_total = Qtot[..., None]
            Qdot_inner = Qinner[..., None]
            Qdot_nucleus = Qnuc[..., None]
            Qdot_outer = Qouter[..., None]
        processed_hdf.append(file)
            
            
        if save_checkpoints and check_index == checkpoint:
            print("Checkpoint reached. Saving...")
            save_hdf(os.path.join(simulation.storage_path, 'Qdot.h5'),
                     ['time', 'Qdot_total', 'Qdot_inner', 'Qdot_nucleus',
                      'Qdot_outer', 'Qdot_radial', 'processed'],
                     [time, Qdot_total, Qdot_inner, Qdot_nucleus, Qdot_outer,
                      Qdot_radial, processed_hdf])
            check_index = 0
        check_index += 1
        progress_index += 1
        findex += 1
        progressBar(progress_index, len(
            simulation.hdf_file_list[start_point:]))
    
    print("Computations done, saving...")
    save_hdf(os.path.join(simulation.storage_path, 'Qdot.h5'),
                     ['time', 'Qdot_total', 'Qdot_inner', 'Qdot_nucleus',
                      'Qdot_outer', 'Qdot_radial', 'processed'],
                     [time, Qdot_total, Qdot_inner, Qdot_nucleus, Qdot_outer,
                      Qdot_radial, processed_hdf])
    return calculate_strain_3D(D, THETA, PHI, time,
                               simulation.cell.radius(simulation.ghost),
                               Qdot_radial, Qdot_total,
                               Qdot_inner, Qdot_nucleus, Qdot_outer)

def spherical_harmonics_gradient(radius, theta, phi):
    """
    Calculates the gradient of the conjugate spherical harmonics times
    the radius for l = 2.
    Returns a list of arrays with dimension (len(phi), lewn(theta),
    len(radius), 3) containing the gradient of the conjugate spherical
    harmonics times the radius from m=-2 to m=2.
    """
    grd = []
    harmonics = SphericalHarmonics()
    for m in range(-2, 3):
        Y2m_r = harmonics.Ylm_conj(m, 2, theta, phi)[..., None] * \
            radius[None, None, :] ** 2
        grd.append(gradient(Y2m_r, radius, theta, phi, 'spherical'))
    return grd

def calculate_Qdot(simulation, gradY, file_name, dV, 
                        inner_rad, igcells, nuc_rad, ngcells):
    """
    Calculates the Qdot for the different regions of the star.
    dot{Q} = dV ρ ∇(v Y*_2m)
    """
    mask_nuc = simulation.cell.radius(simulation.ghost) <= \
        simulation.ghost.remove_ghost_cells_radii(nuc_rad, simulation.dim,
                                                    **ngcells)[..., None]
    mask_inner = (simulation.cell.radius(simulation.ghost) <= \
        simulation.ghost.remove_ghost_cells_radii(inner_rad, simulation.dim, 
                                             **igcells)[..., None] + (20*u.km)) & \
        (np.logical_not(mask_nuc))
    mask_outer = np.logical_not(mask_inner + mask_nuc)
    
    rho = simulation.rho(file_name) * dV
    v_r = simulation.radial_velocity(file_name)
    v_t = simulation.theta_velocity(file_name)
    v_p = simulation.phi_velocity(file_name)
    Qdot = (rho * (v_r * gradY[0][0, ...] + v_t * gradY[0][1, ...] + v_p \
        * gradY[0][2, ...]))
    Qdot_tot = Qdot.sum()[..., None]
    Qdot_inner = Qdot[mask_inner].sum()[..., None]
    Qdot_nuc = Qdot[mask_nuc].sum()[..., None]
    Qdot_outer = Qdot[mask_outer].sum()[..., None]
    Qdot_radial = Qdot.sum(axis=(0,1))[..., None]
    for i in range(1, 5):
        Qdot = (rho * (v_r * gradY[i][0, ...] + v_t * \
            gradY[i][1, ...] + v_p * gradY[i][2, ...]))
        Qdot_tot = np.concatenate((Qdot_tot, Qdot.sum()[..., None]), axis=-1)
        Qdot_inner = np.concatenate((Qdot_inner, Qdot[mask_inner].sum()
                                     [..., None]), axis=-1)
        Qdot_nuc = np.concatenate((Qdot_nuc, Qdot[mask_nuc].sum()[..., None]),
                                  axis=-1)
        Qdot_outer = np.concatenate((Qdot_outer, Qdot[mask_outer].sum()
                                     [..., None]), axis=-1)
        Qdot_radial = np.concatenate((Qdot_radial, Qdot.sum(axis=(0, 1))
                                      [..., None]), axis=-1)
    
    return Qdot_tot, Qdot_inner, Qdot_nuc, Qdot_outer, Qdot_radial

def read_Qdot(simulation):
    """
    Reads the Qdot and masks from a checkpoint file.
    """
    with h5py.File(os.path.join(simulation.storage_path, 'Qdot.h5')) as data:
        time = data['time'][...] * u.s
        Qdot_radial = data['Qdot_radial'][...] * u.g * u.cm ** 2 / u.s
        Qdot_total = data['Qdot_total'][...] * u.g * u.cm ** 2 / u.s
        Qdot_inner = data['Qdot_inner'][...] * u.g * u.cm ** 2 / u.s
        Qdot_nucleus = data['Qdot_nucleus'][...] * u.g * u.cm ** 2 / u.s
        Qdot_outer = data['Qdot_outer'][...] * u.g * u.cm ** 2 / u.s
        if 'processed' in data.keys():
            processed_hdf = data['processed'][...]
        else:
            processed_hdf = None
    return time, Qdot_radial, Qdot_total, Qdot_inner, Qdot_nucleus, Qdot_outer, processed_hdf

def calculate_strain_3D(D, THETA, PHI, time, radius, Qdot_radial, Qdot_total,
                        Qdot_inner, Qdot_nucleus, Qdot_outer):
    if D is not None:
        if not isinstance(D, aerray):
            D = D * u.cm
        const /= D
        add_lb = r''
    else:
        add_lb = r'$\mathcal{D}$'
        D = 1 * u.dimensionless_unscaled
    harmonics = SphericalHarmonics()
    for m in range(5):
        Y22m = harmonics.spin_weighted_Ylm(-2, m-2, 2, THETA, PHI)
        Qdot_radial[:, m, :] = IDL_derivative(time, Qdot_radial[:, m, :]) * \
            Y22m
        Qdot_total[m, :]= IDL_derivative(time, Qdot_total[m, :]) * Y22m
        Qdot_inner[m, :]= IDL_derivative(time, Qdot_inner[m, :]) * Y22m
        Qdot_nucleus[m, :]= IDL_derivative(time, Qdot_nucleus[m, :]) * Y22m
        Qdot_outer[m, :]= IDL_derivative(time, Qdot_outer[m, :]) * Y22m
    const = np.sqrt(2/3) * 8 * np.pi * c.G / (D * c.c ** 4) / u.s ## last u.s accounts for the derivative
    Qdot_radial = const * Qdot_radial.sum(axis=1)
    Qdot_total = const * Qdot_total.sum(axis=0)
    Qdot_inner = const * Qdot_inner.sum(axis=0)
    Qdot_nucleus = const * Qdot_nucleus.sum(axis=0)
    Qdot_outer = const * Qdot_outer.sum(axis=0)
    
    time.set(name='time', label=r'$t-t_\mathrm{b}$', cmap=None,
             limits=[-0.005, time[-1]])
    if np.isclose(THETA, np.pi, 0.05):
        hplus_radial = Qdot_radial.real
        hplus_radial.set(name='hpuls_radial_pol',
                         label=merge_strings(add_lb,
                                             r'$h_\mathrm{+,pol}(r)$'),
                         cmap='seismic',
                         limits=[-3 / D.value, 3 / D.value])
        hcross_radial = -Qdot_radial.imag
        
        hcross_radial.set(name='hcross_radial_pol',
                         label=merge_strings(add_lb,
                                             r'$h_\mathrm{\times,pol}(r)$'),
                         cmap='seismic',
                         limits=[-3 / D.value, 3 / D.value])
        hplus_tot = Qdot_total.real
        hplus_tot.set(name='tot_hplus_pol',
                      cmap='seismic',
                      limits=[-70 / D.value, 70 / D.value],
                      label=merge_strings(add_lb, r'$h_{+,\mathrm{pol,tot}}$'))
        hcross_tot = -Qdot_total.imag
        hcross_tot.set(name='tot_hcross_pol',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{pol,tot}}$'))
        hplus_nuc = Qdot_nucleus.real
        hplus_nuc.set(name='nuc_hplus_pol',
                      cmap='seismic',
                      limits=[-70 / D.value, 70 / D.value],
                      label=merge_strings(add_lb, r'$h_{+,\mathrm{pol,core}}$'))
        hcross_nuc = -Qdot_nucleus.imag
        hcross_nuc.set(name='nuc_hcross_pol',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{pol,core}}$'))
        hplus_inn = Qdot_inner.real
        hplus_inn.set(name='inn_hplus_pol',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{+,\mathrm{pol,conv}}$'))
        hcross_inn = -Qdot_inner.imag
        hcross_inn.set(name='inn_hcross_pol',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{pol,conv}}$'))
        hplus_out = Qdot_outer.real
        hplus_out.set(name='out_hcross_pol',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{+,\mathrm{pol,out}}$'))
        hcross_out = -Qdot_outer.imag
        hcross_out.set(name='out_hcross_pol',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{pol,out}}$'))
    elif np.isclose(THETA, np.pi/2, 0.05):
        hplus_radial = Qdot_radial.real
        hplus_radial.set(name='hpuls_radial_eq',
                         label=merge_strings(add_lb,
                                             r'$h_\mathrm{+,eq}(r)$'),
                         cmap='seismic',
                         limits=[-3 / D.value, 3 / D.value])
        hcross_radial = -Qdot_radial.imag
        
        hcross_radial.set(name='hcross_radial_eq',
                         label=merge_strings(add_lb,
                                             r'$h_\mathrm{\times,eq}(r)$'),
                         cmap='seismic',
                         limits=[-3 / D.value, 3 / D.value])
        hplus_tot = Qdot_total.real
        hplus_tot.set(name='tot_hplus_eq',
                      cmap='seismic',
                      limits=[-70 / D.value, 70 / D.value],
                      label=merge_strings(add_lb, r'$h_{+,\mathrm{eq,tot}}$'))
        hcross_tot = -Qdot_total.imag
        hcross_tot.set(name='tot_hcross_eq',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{eq,tot}}$'))
        hplus_nuc = Qdot_nucleus.real
        hplus_nuc.set(name='nuc_hplus_eq',
                      cmap='seismic',
                      limits=[-70 / D.value, 70 / D.value],
                      label=merge_strings(add_lb, r'$h_{+,\mathrm{eq,core}}$'))
        hcross_nuc = -Qdot_nucleus.imag
        hcross_nuc.set(name='nuc_hcross_eq',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{eq,core}}$'))
        hplus_inn = Qdot_inner.real
        hplus_inn.set(name='inn_hplus_eq',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{+,\mathrm{eq,conv}}$'))
        hcross_inn = -Qdot_inner.imag
        hcross_inn.set(name='inn_hcross_eq',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{eq,conv}}$'))
        hplus_out = Qdot_outer.real
        hplus_out.set(name='out_hcross_eq',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{+,\mathrm{eq,out}}$'))
        hcross_out = -Qdot_outer.imag
        hcross_out.set(name='out_hcross_eq',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{eq,out}}$'))
    else:
        hplus_radial = Qdot_radial.real
        hplus_radial.set(name='hpuls_radial',
                         label=merge_strings(add_lb,
                                             r'$h_\mathrm{+}(r)$'),
                         cmap='seismic',
                         limits=[-3 / D.value, 3 / D.value])
        hcross_radial = -Qdot_radial.imag
        
        hcross_radial.set(name='hcross_radial',
                         label=merge_strings(add_lb,
                                             r'$h_\mathrm{\times}(r)$'),
                         cmap='seismic',
                         limits=[-3 / D.value, 3 / D.value])
        hplus_tot = Qdot_total.real
        hplus_tot.set(name='tot_hplus',
                      cmap='seismic',
                      limits=[-70 / D.value, 70 / D.value],
                      label=merge_strings(add_lb, r'$h_{+,\mathrm{tot}}$'))
        hcross_tot = -Qdot_total.imag
        hcross_tot.set(name='tot_hcross',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{tot}}$'))
        hplus_nuc = Qdot_nucleus.real
        hplus_nuc.set(name='nuc_hplus',
                      cmap='seismic',
                      limits=[-70 / D.value, 70 / D.value],
                      label=merge_strings(add_lb, r'$h_{+,\mathrm{core}}$'))
        hcross_nuc = -Qdot_nucleus.imag
        hcross_nuc.set(name='nuc_hcross',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{core}}$'))
        hplus_inn = Qdot_inner.real
        hplus_inn.set(name='inn_hplus',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{+,\mathrm{conv}}$'))
        hcross_inn = -Qdot_inner.imag
        hcross_inn.set(name='inn_hcross',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{conv}}$'))
        hplus_out = Qdot_outer.real
        hplus_out.set(name='out_hcross',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{+,\mathrm{out}}$'))
        hcross_out = -Qdot_outer.imag
        hcross_out.set(name='out_hcross',
                       cmap='seismic',
                       limits=[-70 / D.value, 70 / D.value],
                       label=merge_strings(add_lb, r'$h_{\times,\mathrm{out}}$'))
    T, R = np.meshgrid(time, radius)
    return aeseries(data=hplus_radial, time=T.copy(), radius=R.copy()), \
           create_series(time, hplus_tot, hplus_nuc, hplus_inn, hplus_out), \
           aeseries(data=hcross_radial, time=T.copy(), radius=R.copy()), \
           create_series(time, hcross_tot, hcross_nuc, hcross_inn, hcross_out)