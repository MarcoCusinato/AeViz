from AeViz.simulation.methods import *
from AeViz.utils.GW_utils import (GW_strain, GWs_energy, calculate_h,
                                  GWs_spectrogram, GWs_peak_indices,
                                  GWs_fourier_transform,
                                  GWs_frequency_peak_indices)
from AeViz.utils.file_utils import load_file, find_column_changing_line
from typing import Literal

"""
Function to process gravitational waves data from a simulation in
spherical coordinates.
These functions are not meant to be used standalone, but rather to be
imported into the Simulation class.
"""

## -----------------------------------------------------------------
## GRAVIATIONAL WAVES DATA
## -----------------------------------------------------------------

@smooth
@derive
@subtract_tob
def GW_Amplitudes(self, distance=None, tob_corrected=True, 
                    zero_correction=True, lower_refinement=False,
                    **kwargs):
    """
    Params:
        distance: distance of the observer from the GW source
        zero_correction: shifts up or down the amplitude to make it
        centred with zero
        lower_refinement: on fine timestep refinements GWs are
        usually noisy, so we take 
        one every n points. n defined as the number to reach 0.1 ms
    Returns GWs amplitudes:
    1D:
        No GWs in spherical symmetry
    2D:
        Column 1: time
        Column 2: + ploarization
    3D:
        Column 1: time
        Column 2: + polarization equatorial plane
        Column 3: + polarization polar plane
        Column 4: x polarization equatorial plane
        Column 5: x polarization polar plane
    """
    data = load_file(self._Simulation__log_path, self._Simulation__grw_path)
    
    if lower_refinement:
        dt = data[1, 2] - data[0, 2]
        new_dt = dt
        n=1
        while new_dt < 5e-5:
            new_dt += dt
            n += 1
        data = data[::n, :]
    
    column_change = find_column_changing_line(self._Simulation__log_path,
                                                self._Simulation__grw_path)
    if zero_correction:
        index = np.argmax((data[:, 2] - self.tob) >= -0.01)
    else:
        index = None
    GWs = GW_strain(self.dim, column_change, data, index, distance)
    if GWs is None:
        return None
    return GWs


def GW_spectrogram(self, distance=None, window_size=aerray(10, u.ms),
                   tob_corrected=True,
                   scale_to:Literal['magnitude', 'psd']='magnitude',
                    **kwargs):
    """
    Parameters:
        distance: distance of the observer from the GW source
        tob_corrected: if the returned timeseries has to be 
        corrected for the tob
        window_size: size of the time window in which to perform the sft,
                    can be aerray or scalar. If scalar the uniit is the
                    one from time
        time_range: list of float or aerrays, crop the signal between
                    the two times
        check_spacing: if true check if the time array is equally spaced,
                       can cause issue for small timesteps
        windowing{'bartlett', 'blackman', 'hamming', 'hanning', 'kaiser'}:
                Window to apply to the signal, default is hann
        scale_to{'magnitude', 'psd'} default magnitude. Each STFT column
                represents either a 'magnitude' or a power spectral
                density ('psd') spectrum
    Returns:
        time: timeseries in s
        frequency: aray of the frequencies in Hz
        Zxx: magnitude
    In 3D simulations:
        h_pl_e, h_pl_p, h_cr_e, h_cr_p
    """
    GW_strain = self.GW_Amplitudes(distance=distance,
                                   tob_corrected=tob_corrected, **kwargs)
    return GWs_spectrogram(self.dim, GW_strain, window_size, scale_to, **kwargs)

def Deltah(self, peak:Literal['bounce', 'highest']='bounce',
            interval=[None, None], min_time=1.75, max_time=2, distance=1, 
            coordinates=False, tob_corrected=True):
    """
    Returns the Delta h of the gravitational wave strain as defined
    in Richers et al. 2017 (https://arxiv.org/pdf/1701.02752.pdf).
    Basically the difference between the first maximum and minimum 
    postbounce, the first peak that appears is not considered.
    In case the highest peak is selected the amplitude returned is 
    the maxima between left and right.
    min and max time are windows in which to search the low and high
    peaks in the strain
    Return:
        amplitude in cm
        if coordinates
            time, h of the highest and lowest peak
    """
    GWs = self.GW_Amplitudes(distance, tob_corrected)
    indices = GWs_peak_indices(GWs, peak, interval, min_time, max_time)
    Deltah = np.abs(GWs[indices[1], 1] - GWs[indices[2], 1])
    if coordinates:
        x = [GWs[indices[1], 0], GWs[indices[2], 0]]
        y = [GWs[indices[1], 1], GWs[indices[2], 1]]
        return Deltah, np.array(x), np.array(y)
    return Deltah

def GWs_peak_frequencies(self, peak:Literal['bounce', 'highest']='bounce',
                            min_time=1.75, max_time=2, interval=[None, None],
                            return_intensities=False, return_fourier=False):
    """
    Calculates the dominant and the second dominant frequency of a GWs peak
    Return:
        frequencies: dominat, second dominant
        if return intensities
            intensities: dominant, second dominant
        if return return fourier
            frequencies array
            htilde
    """
    GWs = self.GW_Amplitudes()
    indices = GWs_peak_indices(GWs, peak, interval, min_time, max_time)
    frequency, htilde = GWs_fourier_transform(GWs, indices)
    indices = GWs_frequency_peak_indices(frequency, htilde)
    return_list = [frequency[indices]]
    if return_intensities:
        return_list.append(htilde[indices])
    if return_fourier:
        return_list.append([frequency, htilde])
    if len(return_list) == 1:
        return return_list[0]
    return return_list

@smooth
def GWs_dE_dt(self, lower_refinement=False, tob_corrected=True):
    """
    Returns the energy carried away by the GWs in erg/s
    """
    GWs = self.GW_Amplitudes(tob_corrected=tob_corrected,
                                lower_refinement=lower_refinement)
    return GWs_energy(GWs, self.dim)

@smooth
@sum_tob
def hydro_strain(self, tob_corrected=True, D=None, theta=np.pi/2, phi=0,
            save_checkpoints=True):
    """
    Calculates the gravitational wave strain from the hydro for a
    simulation
    Returns
        2D
            radius
            time: array of time step
            AE220: len(radius), len(time) array
            full_strain: GWs strain from the full star
            nucleus_strain: GWs strain from the PNS nucleus
            convection_strain: GWs strain from the convection region
            outer_strain: GWs strain from the outer layers
        3D
            radius
            time: array of time step
            [h_+, h_x]:  len(radius), len(time) array
            [h_+, h_x]_full: len(time)
            [h_+, h_x]_nucl: len(time)
            [h_+, h_x]_conv: len(time)
            [h_+, h_x]_out: len(time)
    """
    GW_data = calculate_h(self, D, theta, phi, save_checkpoints)
    if GW_data is None:
        return None        
    return GW_data