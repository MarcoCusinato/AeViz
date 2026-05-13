from AeViz.simulation.methods import *
from AeViz.utils.physics.GW_utils import (GW_strain, calculate_h,
                                  GWs_spectrogram,
                                  universal_modes_relation,
                                  get_spherical_harmonics)
from AeViz.utils.files.string_utils import merge_strings
from AeViz.utils.files.file_utils import (load_file, find_column_changing_line,
                                          load_asd)
from AeViz.spherical_harmonics.spherical_harmonics import SphericalHarmonics
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
@sum_tob
def GW_Amplitudes(self, distance: aerray | None = None,
                  tob_corrected: bool = True, 
                  zero_correction: bool = True,
                  lower_refinement: bool = False,
                  **kwargs) -> aeseries | list[aeseries]:
    """
    Extract the GW strains from a simulation
    
    Parameters
    ----------
    distance : aerray | None, optional
        distance of the observer from the GW source, by default None
    tob_corrected : bool, optional
        centres the time to the core bounce time, by default True
    zero_correction : bool, optional
        sets the beginning of the evolution to zero, by default True
    lower_refinement : bool, optional
        decimates the sampling to 5e-5 s, by default False
    kwargs:
        components to extract

    Returns
    -------
    aeseries | list[aeseries]
        aeseries containing the selected polarisation or the list of all
        polarisations
    """
    if not self._Simulation__grw_path in self._Simulation__loaded_files:
        self._Simulation__loaded_files[self._Simulation__grw_path] = \
            load_file(self._Simulation__log_path, self._Simulation__grw_path)
        self._Simulation__gws_param['col_change'] = \
            find_column_changing_line(self._Simulation__log_path,
                                      self._Simulation__grw_path)
        if 'column_index' in kwargs:
            self._Simulation__gws_param['col'] = kwargs['column_index']
        else:
            self._Simulation__gws_param['col'] = 0
        self._Simulation__gws_param['lower_refinement'] = lower_refinement
        self._Simulation__gws_param['zero_correction'] = zero_correction
    
    low_ref = self._Simulation__gws_param['lower_refinement']
    ze_corr = self._Simulation__gws_param['zero_correction']
    if 'column_index' in kwargs:
        col = kwargs['column_index']
    else:
        col = self._Simulation__gws_param['col']
    ccol = self._Simulation__gws_param['col']
    
    if (self._Simulation__gws is None or low_ref != lower_refinement or
        ze_corr != zero_correction or col != ccol):
        ## get the data, to avoid spikes at the end we neglect the last points
        data = self._Simulation__loaded_files[self._Simulation__grw_path][:-50, :]
        column_change_list = self._Simulation__gws_param['col_change']
        ## update the parameters
        self._Simulation__gws_param['col'] = col
        self._Simulation__gws_param['lower_refinement'] = lower_refinement
        self._Simulation__gws_param['zero_correction'] = zero_correction
    
        if len(column_change_list) > 1:
            column_change = column_change_list[col]
        elif len(column_change_list) == 1:
            column_change = column_change_list[0]
        else:
            column_change = None
        if column_change < 3:
            column_change = None
    
        n = 1
        if lower_refinement:
            dt = data[1, 2] - data[0, 2]
            new_dt = dt
            n=1
            while new_dt < 5e-5:
                new_dt += dt
                n += 1
    
        if zero_correction:
            index = np.argmax((data[:, 2] - self.tob) >= -0.01)
        else:
            index = None
        self._Simulation__gws = GW_strain(self.dim,
                                          column_change,
                                          data,
                                          index,
                                          n,
                                          distance,
                                          self.tob)
    GWs = self._Simulation__gws
    if not 'comp' in kwargs:
        kwargs['comp'] = 'all'
    if distance:
        if not isinstance(distance, aerray):
            distance *= GWs.hple.unit
        GWs.set_distance(distance)
        return GWs.get_polarisation(kwargs['comp'], True)
    else:
        return GWs.get_polarisation(kwargs['comp'], False)

@smooth
@derive
def GWs_dE_df(self,
              comp:Literal['eq', 'h+eq', 'hxeq',
                           'pol' 'h+pol', 'hxpol'] = 'eq',
              time_range: list | None = None,
              window_type:str | None = 'hanning',
              **kwargs) -> aeseries:
    """
    Computes and returns the energy spectra over per frequency for
    a chosen line of sight

    Parameters
    ----------
    comp : Literal['eq', 'h+eq', 'hxeq', 'pol' 'h+pol', 'hxpol'], optional
        line of sight. Can also be given the polarisation and it changes
        automatically to the los, by default 'eq'
    time_range : list | None, optional
        if given cuts considers  only that interval, by default None
    window_type : str, optional
        if not None applies a window to the signal before
        performing a FFT, by default 'hanning'
    **kwargs

    Returns
    -------
    aeseries
        energy spectra and frequency for that los

    """
    if self._Simulation__gws is None:
        self.GW_Amplitudes()
    GWs = self._Simulation__gws
    kwargs['apply_window'] = True if window_type is not None else False
    kwargs['window_type'] = window_type
    
    GWs.set_fft_config(**kwargs)
    if 'eq' in comp:
        los = 'eq'
    else:
        los = 'pol'
    if time_range:
        istart = np.argmax(GWs.time >= time_range[0])
        istop = np.argmax(GWs.time >= time_range[1])
        return GWs[istart:istop].get_dEdf(los=los)
    return GWs.get_dEdf(los=los)

@smooth
@derive
def hchar(self,
          comp:Literal['eq', 'h+eq', 'hxeq',
                       'pol' 'h+pol', 'hxpol'] = 'eq',
          time_range: list | None = None,
          window_type:str | None = 'hanning',
          distance: aerray = (10 * u.kpc),
          divide_by_frequency: bool = True,
          type: str = 'fft',
          **kwargs) -> aeseries:
    """
    Computes the characteristic strain of the GW signal

    Parameters
    ----------
    comp : Literal['eq', 'h+eq', 'hxeq', 'pol' 'h+pol', 'hxpol'], optional
        line of sight. Can also be given the polarisation and it changes
        automatically to the los, by default 'eq'
    time_range : list | None, optional
        if given cuts considers  only that interval, by default None
    window_type : str, optional
        if not None applies a window to the signal before
        performing a FFT, by default 'hanning'
    distance : aerray, optional
        distance of the observer from the GW source, 
        by default (10 * u.kpc)
    divide_by_frequency : bool, optional
        if True returns the characteristic strain divided by the sqrt of
        the frequency, by default True
    type : str, optional
        type of characteristic strain to return 'fft' or 'energy,
        by default 'fft'

    Returns
    -------
    aeseries
        characteristic strain and frequency
    """

    if self._Simulation__gws is None:
        self.GW_Amplitudes()
    GWs = self._Simulation__gws
    GWs.set_distance(distance)
    kwargs['apply_window'] = True if window_type is not None else False
    kwargs['window_type'] = window_type
    
    GWs.set_fft_config(**kwargs)
    if 'eq' in comp:
        los = 'eq'
    else:
        los = 'pol'
    if time_range:
        istart = np.argmax(GWs.time >= time_range[0])
        istop = np.argmax(GWs.time >= time_range[1])
        return GWs[istart:istop].get_characteristic_strain(los=los,
                                                           mode=type,
                                                           divide_by_frequency=divide_by_frequency)
    return GWs.get_characteristic_strain(los=los,
                                         mode=type,
                                         divided_by_frequency=divide_by_frequency)

@sum_tob
def GW_spectrogram(self,
                   distance: aerray | None = None,
                   window_size: aerray = aerray(10, u.ms),
                   tob_corrected: bool = True,
                   scale_to:Literal['magnitude', 'psd'] = 'magnitude',
                    **kwargs) -> aeseries | list[aeseries]:
    """
    Computes the spectrograms of the GW strain(s).

    Parameters
    ----------
    distance : aerray | None, optional
        distance of the observer from the GW source, by default None
    window_size : aerray, optional
        size of the window for the STFT, by default aerray(10, u.ms)
    tob_corrected : bool, optional
        centres the time to the core bounce time, by default True
    scale_to : Literal['magnitude', 'psd'], optional
        how to scale the STFT, by default 'magnitude'

    Returns
    -------
    aeseries | list[aeseries]
        list of series containing the spectrograms of the four or single
        polarisations. Each aeseries contains:
            time: timeseries in s
            frequency: aray of the frequencies in Hz
            Zxx: magnitude
        In the 3D simulation the polarisations are ordered as:
            h_pl_e, h_pl_p, h_cr_e, h_cr_p
    """
    if self._Simulation__gws is None:
        self.GW_Amplitudes() 
    GW_strain = self._Simulation__gws
    if distance is not None:
        GW_strain.set_distance(distance)
        return GW_spectrogram(self.dim, GW_strain.get_polarisation(comp='all',
                                                                   dimensionless=True),
                              window_size, scale_to, **kwargs)
    return GWs_spectrogram(self.dim, GW_strain.get_polarisation(comp='all'),
                           window_size, scale_to, **kwargs)

def Deltah(self,
           peak: Literal['bounce', 'max'] ='bounce',
           comp: Literal['h+eq', 'h+pol','hxeq', 'hxpol']='h+eq',
           time_range: list[aerray|float] | None = None,
           tol: float = 0.05,
           distance: aerray | None = None,
           return_indices: Literal[True, False] = False,
           detrend: bool = True) -> aerray | tuple[aerray, list[int]]:
    """
    Returns the Delta h of the gravitational wave strain as defined
    in Richers et al. 2017 (https://arxiv.org/pdf/1701.02752.pdf).
    Basically the difference between the first maximum and minimum 
    postbounce, the first peak that appears is not considered.
    In case the highest peak is selected the amplitude returned is 
    the maxima between left and right.
    
    Parameters
        ----------
        peak : Literal['bounce', 'max'], optional
            peak to find, if bounce is selected the time interval is 
            neglected, otherwise the most prominent in the selected 
            interval will be considered by default 'bounce'
        comp : str, optional
            the strain to use, possibilities ['h+eq', 'h+pol',
            'hxeq', 'hxpol', 'all'], by default 'h+eq'
        time_range : list[aerray | float] | None, optional
            time range to considered to find the peak. If None the full
            strain is analysed, by default None
        tol : float, optional
            used by the maximum peak. If a secondary crest of valley if
            lower (in absolute value) of this percentage of the maximum
            then it is treated as no intersection with zero was found,
            by default 0.05
        distance : aerray | None, optional
            returns the value at that given distance,
            by default None
        return_indices : Literal[False], optional
            If the indices have to be returned, by default False
        detrend : bool, optional
            If the strain has to be detrended of any  memory effect before
            proceeding

        Returns
        -------
        aerray | tuple[aerray, list[int]]
            returns the value of the range strain width as an aerray. If
            the flag return_indices is set to True returns also the indices
            of the intersections
    """
    if self._Simulation__gws is None:
        self.GW_Amplitudes()
    
    GW_strain = self._Simulation__gws
    if detrend:
        GW_strain = GW_strain.copy()
        GW_strain.detrend()
    if distance is not None:
        GW_strain.set_distance(distance)
        dimensionless = True
    else:
        dimensionless = False
        
    return GW_strain.range_strain_width(peak = peak,
                                        comp = comp,
                                        time_range = time_range,
                                        tol = tol,
                                        dimensionless = dimensionless,
                                        return_indices = return_indices)

def GWs_peak_frequencies(self, 
                         peak: Literal['bounce', 'max'] = 'bounce',
                         comp: Literal['h+eq', 'h+pol', 'hxeq', 'hxpol'] = 'h+eq',
                         time_range: list[aerray|float] | None = None,
                         tol: float = 0.05,
                         normalised: bool = False,
                         return_max: Literal[False]=False,
                         detrend: bool = True) -> aeseries | aerray:
    """
    Computes the frequency of of the oscillation associated with the
    strain range width.
    
    Parameters
    ----------
    peak : Literal['bounce', 'max'], optional
        peak to find, if bounce is selected the time interval is 
        neglected, otherwise the most prominent in the selected 
        interval will be considered by default 'bounce'
    comp : str, optional
        the strain to use, possibilities ['h+eq', 'h+pol',
        'hxeq', 'hxpol', 'all'], by default 'h+eq'
    time_range : list[aerray | float] | None, optional
        time range to considered to find the peak. If None the full
        strain is analysed, by default None
    tol : float, optional
        used by the maximum peak. If a secondary crest of valley if
        lower (in absolute value) of this percentage of the maximum
        then it is treated as no intersection with zero was found,
        by default 0.05
    normalised : bool, optional
        the fourier transform is normalised by its maximum,
        by default False
    return_max : Literal[False], optional
        If True returns only the maximum value of the frequency. In
        the other case the full Fourier transform is returned.
        by default False
    detrend : bool, optional
            If the strain has to be detrended of any  memory effect before
            proceeding

    Returns
    -------
    aeseries | aerray
        aseseries of the Fourier transform containing the frequency
        and the absolute value of the spectrum. aerray of the maximum
        frequency of the oscillation
    """
    if self._Simulation__gws is None:
        self.GW_Amplitudes()
    
    GW_strain = self._Simulation__gws
    if detrend:
        GW_strain = GW_strain.copy()
        GW_strain.detrend()
    return GW_strain.frequency_of_peak(peak = peak,
                                       comp = comp,
                                       time_range = time_range,
                                       normalised = normalised,
                                       return_max = return_max)

@smooth
@derive
@sum_tob
def GWs_luminosity(self,
              tob_corrected: bool = True,
              **kwargs) -> aeseries:
    """
    Computes the GW luminosity in erg/s

    Parameters
    ----------
    tob_corrected : bool, optional
        centres the time to the core bounce time, by default True

    Returns
    -------
    aeseries
        Contains the GW luminosity evolution over time
    """
    if self._Simulation__gws is None:
        self.GW_Amplitudes() 
    GW_strain = self._Simulation__gws
    return GW_strain.get_luminosity()

@smooth
@sum_tob
def GWs_energy(self,
               tob_corrected: bool = True,
              **kwargs) -> aeseries:
    """
    Computes the energy carried away as GWs

    Parameters
    ----------
    tob_corrected : bool, optional
        centres the time to the core bounce time, by default True

    Returns
    -------
    aeseries
        contains the evolution time and the GW energy.
    """
    if self._Simulation__gws is None:
        self.GW_Amplitudes() 
    GW_strain = self._Simulation__gws
    return GW_strain.get_energy()
 
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