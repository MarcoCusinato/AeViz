import numpy as np
from AeViz.units.units import units
from AeViz.utils.math_utils import IDL_derivative
from scipy.signal import stft
from numpy.fft import fft, fftfreq

u = units()

## ---------------------------------------------------------------------
## GW strain
## ---------------------------------------------------------------------

def GW_strain(sim_dim, column_change, data, index):
    assert sim_dim in [1, 2, 3], "Simulation MUST be 1, 2 or 3D."
    if sim_dim == 1:
        return GW_strain_1D(data)
    elif sim_dim == 2:
        return correct_zero(2, GW_strain_2D(data), index)
    else:
        GWs = GW_strain_3D(data)
        if column_change is not None:
            GWs[:column_change, 1] = GW_strain_2D(data)[:column_change, 1]
            GWs[:column_change,2:] = 0
        return correct_zero(3, GWs, index)

def GW_strain_1D(data):
    print("No GW for you :'(")
    return None

def GW_strain_2D(data):
    """
        GWs amplitudes calculated as the first partial time derivative
        of NE_220 and not with the second partial time derivative of 
        ME_220. Moreover, we only consider matter contribution to the 
        amplitude, being the main contribution to it.
    """
    const = -0.125 *  np.sqrt(15/np.pi)
    return np.stack((data[:, 2], const * IDL_derivative(data[:,2], data[:,5])),
              axis = 1)

def GW_strain_3D(data):
    const = 1 #8 * np.sqrt(np.pi / 15)
    hD_pl_p = 2. * ( data[:,9] - data[:,13] )
    hD_pl_e = 2. * ( data[:,17] - data[:,13] )
    hD_cr_p = 2. * ( data[:,10] + data[:,12] )
    hD_cr_e = 2. * ( - data[:,14] - data[:,16] )

    hD_pl_p = const * IDL_derivative( data[:,2], hD_pl_p )
    hD_cr_p = const * IDL_derivative( data[:,2], hD_cr_p )
    hD_pl_e = const * IDL_derivative( data[:,2], hD_pl_e )
    hD_cr_e = const * IDL_derivative( data[:,2], hD_cr_e )

    return np.stack((data[:,2], hD_pl_e, hD_pl_p, hD_cr_e, hD_cr_p), axis = -1)

def correct_zero(sim_dim, GWs, index):
    if sim_dim == 1:
        pass
    else:
        if index is None:
            return GWs
        for i in range(1, GWs.shape[1]):
            GWs[:, i] -= GWs[:index, i].mean()
        return GWs
    
## ---------------------------------------------------------------------
## GW spectrogram
## ---------------------------------------------------------------------


def GW_spectrogram(sim_dim, GWs, window_size):
    assert sim_dim in [1, 2, 3], "Simulation MUST be 1, 2 or 3D."
    if sim_dim == 1:
        return GW_spectrogram_1D(GWs)
    elif sim_dim == 2:
        return GW_spectrogram_2D(GWs, window_size)
    else:
        return GW_spectrogram_3D(GWs, window_size)

def GW_spectrogram_1D(GWs):
    print("And also no spectrogram for you :'(\nごめんなさい")
    return None

def GW_spectrogram_2D(GWs, window_size):
    """
        
    """
    fs = (GWs[1, 0] - GWs[0, 0]) ** (-1)
    frequency, time, Zxx = stft(GWs[:,1], fs = fs, nperseg=window_size)
    time = time / time[-1] * GWs[-1, 0]
    return time, frequency, np.abs(Zxx)

def GW_spectrogram_3D(GWs, window_size):
    fs = (GWs[1, 0] - GWs[0, 0]) ** (-1)
    freq_pl_e, time_pl_e, Zxx_pl_e = stft(GWs[:,1], fs=fs, nperseg=window_size)
    time_pl_e = time_pl_e / time_pl_e[-1] * GWs[-1, 0]
    freq_pl_p, time_pl_p, Zxx_pl_p = stft(GWs[:,2], fs=fs, nperseg=window_size)
    time_pl_p = time_pl_p / time_pl_p[-1] * GWs[-1, 0]
    freq_cr_e, time_cr_e, Zxx_cr_e = stft(GWs[:,3], fs=fs, nperseg=window_size)
    time_cr_e = time_cr_e / time_cr_e[-1] * GWs[-1, 0]
    freq_cr_p, time_cr_p, Zxx_cr_p = stft(GWs[:,4], fs=fs, nperseg=window_size)
    time_cr_p = time_cr_p / time_cr_p[-1] * GWs[-1, 0]
    
    return np.stack((time_pl_e, time_pl_p, time_cr_e, time_cr_p), axis=-1), \
           np.stack((freq_pl_e, freq_pl_p, freq_cr_e, freq_cr_p), axis=-1), \
           np.stack((np.abs(Zxx_pl_e), np.abs(Zxx_pl_p), np.abs(Zxx_cr_e),
                     np.abs(Zxx_cr_p)), axis=-1)

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
        x_max = np.argmax(GWs[:, 1])
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
        $\tilde{h} * \sqrt{freq}$
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

def GWs_frequency_peak_indices(frequency, htilde):
    """
    Method that find the indices of the first and second frequency peak
    on a fourier transformed GW strain
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
## TODO: fix with the new decorators

def calculate_AE220(simulation):
    """
    Calculates the AE220 from density and velocities for a
    2D simulation.
    ONLY 2D
    Returns
        time: arry of time step
        AE220: len(radius), len(time) array
    """
    radius = simulation.cell.radius(simulation.ghost)
    dV = -simulation.cell.dVolume_integration(simulation.ghost)
    costheta = np.cos(simulation.cell.theta(simulation.ghost))
    file_list = simulation.file_list_hdf()
    inner_radius, time, gh_inner = simulation.get_innercore_radius(innercore_radius=True, ret_time=True,
                                                    ghost_cells=True)
    convective_radius, gh_conv = simulation.get_convective_radius(convective_radius=True,
                                                            ghost_cells=True)
    NE220_full = np.zeros((len(radius), len(time)))
    NE220_inner = np.zeros((len(radius), len(time)))
    NE220_convection = np.zeros((len(radius), len(time)))
    NE220_outer = np.zeros((len(radius), len(time)))
    radius = radius[None, ...]
    costheta = costheta[..., None]
    for index in range(len(file_list)):
        # Get all the quantities needed for the calculation
        try:
            data_h5 = simulation.open_h5(file_list[index])
        except:
            NE220_full[:, index] = NE220_full[:, index - 1]
            NE220_inner[:, index] = NE220_inner[:, index - 1]
            NE220_convection[:, index] = NE220_convection[:, index - 1]
            NE220_outer[:, index] = NE220_outer[:, index - 1]
            continue
        rho = simulation.rho(data_h5)
        vR = simulation.radial_velocity(data_h5)
        vT = simulation.theta_velocity(data_h5)
        simulation.close_h5(data_h5)
        # Calculate the NE220
        # NE220 of the full space
        NE220 = dV * radius * rho *  ( vR * \
            ( 3 * costheta ** 2 - 1 ) - 3 * vT * costheta * np.sqrt( 1 - costheta ** 2 ) )
        # Create masks
        mask_inner = radius <= simulation.ghost.remove_ghost_cells_radii(convective_radius[..., index],
                                                                        simulation.dim, **gh_conv)[..., None]
        mask_convection = (radius <= simulation.ghost.remove_ghost_cells_radii(inner_radius[..., index],
                                                                            simulation.dim, **gh_inner)[..., None] + 2e6) & \
                            ( np.logical_not(mask_inner) )
        mask_outer = np.logical_not(mask_inner + mask_convection)
        # Calculate the NE220s
        NE220_full[:, index] = np. sum( NE220, axis = 0 )
        NE220_inner[:, index] = np. sum( NE220 * mask_inner, axis = 0 )
        NE220_convection[:, index] = np. sum( NE220 * mask_convection, axis = 0 )
        NE220_outer[:, index] = np. sum( NE220 * mask_outer, axis = 0 )
    
    const = ( u.G * 16 * np.pi ** 0.5 / (np.sqrt( 15 ) * u.speed_light ** 4 ) )

    return time, IDL_derivative(time, NE220_full * const), calculate_strain(NE220_full, time), \
            calculate_strain(NE220_inner, time), calculate_strain(NE220_convection, time), \
            calculate_strain(NE220_outer, time)
            

def calculate_strain(NE220, time):
    """
    Calculates the strain from the NE220
    """
    # Put again the correct units
    NE220 *= ( u.G * 16 * np.pi ** 0.5 / (np.sqrt( 15 ) * u.speed_light ** 4 ) )
    # Calculate the strain
    strain = IDL_derivative(time, NE220)
    return strain.sum(axis=0)
