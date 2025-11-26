from AeViz.simulation.methods import *
from AeViz.utils.physics.GW_utils import (GW_strain, GWs_energy, calculate_h,
                                  GWs_spectrogram, GWs_peak_indices,
                                  GWs_fourier_transform,
                                  GWs_frequency_peak_indices,
                                  characteristic_strain,
                                  GWs_energy_per_frequency,
                                  universal_modes_relation,
                                  get_spherical_harmonics)
from AeViz.utils.files.string_utils import merge_strings
from AeViz.utils.files.file_utils import (load_file, find_column_changing_line,
                                          load_asd)
from AeViz.spherical_harmonics.spherical_harmonics import SphericalHarmonics
from typing import Literal
from AeViz.simulation.methods.other_spherical_sym import radial_profile

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
    
    n = 1
    if lower_refinement:
        dt = data[1, 2] - data[0, 2]
        new_dt = dt
        n=1
        while new_dt < 5e-5:
            new_dt += dt
            n += 1
    
    if 'column_index' in kwargs:
        col = kwargs['column_index']
    else:
        col = None
    
    column_change = find_column_changing_line(self._Simulation__log_path,
                                                self._Simulation__grw_path,
                                                col)
    if zero_correction:
        index = np.argmax((data[:, 2] - self.tob) >= -0.01)
    else:
        index = None
    if 'return_components' in kwargs and self.dim == 3:
        GWs = GW_strain(self.dim, column_change, data, index, n, distance,
                        kwargs['return_components'])
    else:
        GWs = GW_strain(self.dim, column_change, data, index, n, distance)
    if GWs is None:
        return None
    if self.dim > 2:
        if 'comp' in kwargs:
            if kwargs['comp'] == 'all':
                pass
            elif kwargs['comp'] == 'h+eq':
                return GWs[0]
            elif kwargs['comp'] == 'h+pol':
                return GWs[1]
            elif kwargs['comp'] == 'hxeq':
                return GWs[2]
            elif kwargs['comp'] == 'hxpol':
                return GWs[3]
            else:
                raise TypeError("GW component not recognized")
    return GWs

@smooth
@derive
def GWs_dE_df(self, tob_corrected=True, time_range=None, windowing='hanning',
              **kwargs):
    kw = {'zero_correction': True,
          'lower_refinement': False}
    if 'lower_refinement' in kwargs:
        kw['lower_refinement'] = kwargs['lower_refinement']
    if 'zero_correction' in kwargs:
        kw['zero_correction'] = kwargs['zero_correction']
    GW_strain = self.GW_Amplitudes(tob_corrected=tob_corrected, comp='all',
                                   **kw)
    dedf = GWs_energy_per_frequency(GW_strain, self.dim, time_range, windowing)
    if self.dim > 2:
        if 'comp' in kwargs:
            if kwargs['comp'] == 'all':
                pass
            elif kwargs['comp'] == 'h+eq':
                return dedf[0]
            elif kwargs['comp'] == 'h+pol':
                return dedf[1]
            elif kwargs['comp'] == 'hxeq':
                return dedf[2]
            elif kwargs['comp'] == 'hxpol':
                return dedf[3]
            else:
                raise TypeError("GW component not recognized")
    return dedf

@smooth
@derive
def hchar(self, tob_corrected=True, time_range=None, windowing='hanning',
          distance=(10 * u.kpc), divide_by_frequency=True, **kwargs):
    kw = {'zero_correction': True,
          'lower_refinement': False}
    if 'lower_refinement' in kwargs:
        kw['lower_refinement'] = kwargs['lower_refinement']
    if 'zero_correction' in kwargs:
        kw['zero_correction'] = kwargs['zero_correction']
    GW_strain = self.GW_Amplitudes(tob_corrected=tob_corrected, comp='all',
                                   **kw)
    hchar = characteristic_strain(GW_strain, self.dim, time_range, windowing,
                                  distance, divide_by_frequency)
    if self.dim > 2:
        if 'comp' in kwargs:
            if kwargs['comp'] == 'all':
                pass
            elif kwargs['comp'] == 'h+eq':
                return hchar[0]
            elif kwargs['comp'] == 'h+pol':
                return hchar[1]
            elif kwargs['comp'] == 'hxeq':
                return hchar[2]
            elif kwargs['comp'] == 'hxpol':
                return hchar[3]
            else:
                raise TypeError("GW component not recognized")
    return hchar

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
@subtract_tob
def GWs_dE_dt(self, lower_refinement=False, tob_corrected=True, **kwargs):
    """
    Returns the energy carried away by the GWs in erg/s
    """
    if self.dim == 3:
        GWs = self.GW_Amplitudes(tob_corrected=False,
                                lower_refinement=lower_refinement,
                                return_components=True)
    else:
        GWs = self.GW_Amplitudes(tob_corrected=False,
                                lower_refinement=lower_refinement)
    
    return GWs_energy(GWs, self.dim)

@smooth
@sum_tob
def hydro_strain(self, tob_corrected=True, D=None, theta=np.pi/2, phi=0,
                 comp:Literal['h+eq', 'hxeq', 'h+pol', 'hxpol']=None,
                 save_checkpoints=True, **kwargs):
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
    if self.dim == 1:
        return None
    elif self.dim == 2:
        return calculate_h(self, D, theta, phi, save_checkpoints, **kwargs)
    elif self.dim == 3:
        if comp is None:
            return calculate_h(self, D, theta, phi, save_checkpoints, **kwargs)
        elif comp == 'h+eq':
            return calculate_h(self, D, np.pi/2, 0, save_checkpoints, **kwargs)[0:2]
        elif comp == 'hxeq':
            return calculate_h(self, D, np.pi/2, 0, save_checkpoints, **kwargs)[2:]
        elif comp == 'h+pol':
            return calculate_h(self, D, np.pi, 0, save_checkpoints, **kwargs)[:2]
        elif comp == 'hxpol':
            return calculate_h(self, D, np.pi, 0, save_checkpoints, **kwargs)[2:]

@get_grid
def hydro_strain_2D(self, file_name, **kwargs):
    ## find the two file indices
    if isinstance(file_name, int):
        file_1 = file_name
    elif isinstance(file_name, float) or isinstance(file_name, aerray):
        file_1 = self.hdf_file_list.index(self.find_file_from_time(file_name))
    else:
        file_1 = self.hdf_file_list.index(file_name)
    if file_1 == 0:
        file_0 = 0
    else:
        file_0 = file_1 - 1 
    kwargs.setdefault('comp', 'h+eq')
    ## set up constant
    if self.dim == 1:
        const = 1
    elif self.dim == 2:
        const =  -0.125 *  np.sqrt(15/np.pi) * \
        (c.G * 8 * np.pi ** 0.5 / (np.sqrt( 15 ) * c.c ** 4))
    else:
        const = np.sqrt(2/3) * 8 * np.pi * c.G / (c.c ** 4 * 5) / u.s
    kwargs.setdefault('D', None)
    if kwargs['D'] is not None:
        if not isinstance(kwargs['D'], aerray):
            kwargs['D'] = kwargs['D'] * u.cm
        const /= kwargs['D']
        add_lb = r''
    else:
        add_lb = r'$\mathcal{D}$'
        kwargs['D'] = 1 * u.dimensionless_unscaled
    
    if self.dim == 1:
        return None
    elif self.dim == 2:
        radius = self.cell.radius(self.ghost)
        dV = self.cell.dVolume_integration(self.ghost)
        ctheta = np.cos(self.cell.theta(self.ghost))[:, None]
        rho = self.rho(file_1)
        vr = self.radial_velocity(file_1)
        vt = self.theta_velocity(file_1)
        t1 = self.time(file_1)
        NE220_1 = dV * radius * rho * (vr * (3 * ctheta ** 2 - 1) - \
            3 * vt * ctheta * np.sqrt(1 - ctheta ** 2))
        rho = self.rho(file_0)
        vr = self.radial_velocity(file_0)
        vt = self.theta_velocity(file_0)
        t0 = self.time(file_0)
        NE220_0 = dV * radius * rho * (vr * (3 * ctheta ** 2 - 1) - \
            3 * vt * ctheta * np.sqrt(1 - ctheta ** 2))
        GWs = (NE220_1 - NE220_0) / (t1 - t0) * const
        kwargs['D'] = kwargs['D'].value if isinstance(kwargs['D'], aerray) else kwargs['D']
        GWs.set(name='AE220', label=merge_strings(add_lb, r'$A^{E2}_{20}(r)$'),
              cmap='seismic', limits=[-0.1 / kwargs['D'].value, 0.1 / kwargs['D'].value])
    else:
        if kwargs['comp'] in ['h+eq', 'hxeq']:
            THETA = np.pi / 2
        elif kwargs['comp'] in ['h+pol', 'hxpol']:
            THETA = np.pi
        dOmega = self.cell.dOmega(self.ghost)
        dV = self.cell.dVolume_integration(self.ghost)
        gradY, _ = get_spherical_harmonics(
                    self.cell.radius(self.ghost),
                    self.cell.theta(self.ghost),
                    self.cell.phi(self.ghost),
                    dOmega)
        rho = self.rho(file_1) * dV
        vr = self.radial_velocity(file_1)
        vt = self.theta_velocity(file_1)
        vp = self.phi_velocity(file_1)
        t1 = self.time(file_1)
        Qdot_1 = (rho * (vr * gradY[0][0, ...] + vt * gradY[0][1, ...] + vp \
            * gradY[0][2, ...]))[..., None]
        for i in range(1, 5):
            Qdot_1 = np.concatenate((Qdot_1, (rho * (vr * gradY[i][0, ...] + vt * \
                gradY[i][1, ...] + vp * gradY[i][2, ...]))[..., None]), axis=-1)
        rho = self.rho(file_0) * dV
        vr = self.radial_velocity(file_0)
        vt = self.theta_velocity(file_0)
        vp = self.phi_velocity(file_0)
        t0 = self.time(file_0)
        Qdot_0 = (rho * (vr * gradY[0][0, ...] + vt * gradY[0][1, ...] + vp \
            * gradY[0][2, ...]))[..., None]
        for i in range(1, 5):
            Qdot_0 = np.concatenate((Qdot_0, (rho * (vr * gradY[i][0, ...] + vt * \
                gradY[i][1, ...] + vp * gradY[i][2, ...]))[..., None]), axis=-1)
        dt = t1 - t0
        harmonics = SphericalHarmonics()
        for m in range(0, 5):
            Y22m = harmonics.spin_weighted_Ylm(-2, m-2, 2, THETA, 0)
            Qdot_1[..., m] = (Qdot_1[..., m] - Qdot_0[..., m]) / dt * Y22m
        GWs = const * Qdot_1.sum(axis=-1)
        kwargs['D'] = kwargs['D'].value if isinstance(kwargs['D'], aerray) else kwargs['D']
        if kwargs['comp'] == 'h+eq':
            GWs = GWs.real
            GWs.set(name='hpluseq',
                         label=merge_strings(add_lb,
                                             r'$h_\mathrm{+,eq}$'),
                         cmap='seismic',
                         limits=[-0.1 / kwargs['D'], 0.1 / kwargs['D']])
        elif kwargs['comp'] == 'hxeq':
            GWs = -GWs.imag
            GWs.set(name='htimeseq',
                         label=merge_strings(add_lb,
                                             r'$h_\mathrm{\times,eq}$'),
                         cmap='seismic',
                         limits=[-0.1 / kwargs['D'], 0.1 / kwargs['D']])
        elif kwargs['comp'] == 'h+pol':
            GWs = GWs.real
            GWs.set(name='hpluspol',
                         label=merge_strings(add_lb,
                                             r'$h_\mathrm{+,pol}$'),
                         cmap='seismic',
                         limits=[-0.1 / kwargs['D'], 0.1 / kwargs['D']])
        elif kwargs['comp'] == 'hxpol':
            GWs = -GWs.imag
            GWs.set(name='hcrosspol',
                         label=merge_strings(add_lb,
                                             r'$h_\mathrm{\times,pol}$'),
                         cmap='seismic',
                         limits=[-0.1 / kwargs['D'], 0.1 / kwargs['D']])
    return GWs

def ASD(self, detector, **kwargs):
    """
    Return the theoretical ASD for the selected detector
    """
    return load_asd(self.utils_path, detector)

@smooth
@sum_tob
def modes_universal_relations(self,
                              mode: Literal['2f_torres', '2p1_torres',
                                            '2p2_torres', '2p3_torres',
                                            '2g1_torres', '2g3_torres'],
                              tob_corrected=True, **kwargs):
    radius = self.PNS_radius(rad='avg')
    mass = self.PNS_mass_ene(comp='mass')
    if mode == '2g3_torres':
        rhoC = self.radial_profile('rho')
        pC = self.radial_profile('gas_pressure')
    else:
        rhoC = None
        pC = None
    return universal_modes_relation(mass, radius, mode, rhoC=rhoC, pC=pC)