from __future__ import annotations
from AeViz.units import u, aerray, aeseries
from AeViz.units.constants import constants as c
import numpy as np
from typing import Literal, overload
import scipy.signal.windows
import scipy.signal
from scipy.optimize import curve_fit
import warnings
from numpy.fft import rfftfreq, rfft
from AeViz.utils.files.string_utils import merge_strings, apply_symbol
from AeViz.utils.math_utils import IDL_derivative
import copy
import inspect

class GWstrain:
    """
    Class for the manimulation of gravitational waves generated from
    core collapse supernova (CCSN) explosion. It may contain just the
    equatorial '+' polarisation in case of 2D simulations or two,
    independent polarisations '+' and 'x' for polar and equatorial
    lines of sight in the 3D case.

    Parameters
    ----------
    data : np.ndarray, optional
        gravitational waves data, rows representing each timestep and
        colums oredered as time, h+eq, h+pol, hxeq, hxpol. If None, an
        empty object is created by default None
    time_units : aerray, optional
        the units to use for the time array, by default u.s
    strain_units : aerray, optional
        the units to use for the straint, if it is the dimensionless
        strain please use u.dimensionless_unscaled, by default u.cm
    distance : aerray, optional
        the distance at which the GWs are evaluated, by default
        (10 * u.kpc)
    tob_corrected : bool, optional
        if the time is already corrected by the time of bounce, ie, the
        time 0 is when bounce happens, by default True
    tob : aerray, optional
        the time of bounce, if it is None and tob_corrected is True it
        is set as np.abs(time[0]), by default None
    **kwargs : nd.array, optional
        the second derivatives of the components of the quadrupole
        moment they are needed to compute the energy in 3D. They must be
        named as tij, with ij the cartesian components
    """

    def __init__(self,
                 data:np.ndarray = None,
                 time_units:aerray = u.s,
                 strain_units:aerray = u.cm,
                 distance:aerray = (10 * u.kpc),
                 tob_corrected:bool = True,
                 tob:aerray = None,
                 **kwargs              
                 ) -> None:
        
        assert (data.ndim == 2 and data.shape[1] in [2, 5]), \
            "The data must be ordered as time, h+eq, "\
                "(h+pol, hxeq, hxpol)."
        self.time = data[:, 0] * time_units
        if tob is None and tob_corrected:
            self.tob = np.abs(self.time[0])
            self.tob.set(name='Time of bounce', label=r'$t_\mathrm{b}$')
        elif tob is not None:
            self.tob = tob.to(self.time.unit)
            self.tob.set(name='Time of bounce', label=r'$t_\mathrm{b}$')
        else:
            self.tob = None
        if tob_corrected:
            self.time.set(name='time', label=r'$t-t_\mathrm{b}$',
                          limits=[-0.005, self.time[-1].value])
            self.tob_corrected = True
        elif not tob_corrected and self.tob is not None:
            self.time -= self.tob
            self.time.set(name='time', label=r'$t-t_\mathrm{b}$',
                          limits=[-0.005, self.time[-1].value])
            self.tob_corrected = True
        else:
            self.time.set(name='time', label=r'$t$',
                          limits=[-0.005, self.time[-1].value])
            self.tob_corrected = False
        
        if strain_units == u.dimensionless_unscaled and distance is None:
            raise TypeError("The strain cannot be dimensionless if no "\
                "distance is provided.")
        self.distance = distance
        if data.shape[1] == 2:
            self.hple = data[:, 1] * strain_units
            self.hplp = np.zeros(len(self.hple)) * strain_units
            self.hcre = np.zeros(len(self.hple)) * strain_units
            self.hcrp = np.zeros(len(self.hple)) * strain_units
            self.sim_dim = 2
        elif data.shape[1] == 5:
            self.hple = data[:, 1] * strain_units
            self.hplp = data[:, 2] * strain_units
            self.hcre = data[:, 3] * strain_units
            self.hcrp = data[:, 4] * strain_units
            self.sim_dim = 3

        if strain_units == u.dimensionless_unscaled:
            self.hple *= self.distance.to(u.cm)
            self.hplp *= self.distance.to(u.cm)
            self.hcre *= self.distance.to(u.cm)
            self.hcrp *= self.distance.to(u.cm)
        else:
            self.hple = self.hple.to(u.cm)
            self.hplp = self.hplp.to(u.cm)
            self.hcre = self.hcre.to(u.cm)
            self.hcrp = self.hcrp.to(u.cm)

        self.hple.set(name='GW strain h+eq',
                      label=r'$\mathcal{D}h_{+,\mathrm{eq}}$',
                      limits=[-150, 150])
        self.hplp.set(name='GW strain h+pol',
                      label=r'$\mathcal{D}h_{+,\mathrm{pol}}$',
                      limits=[-150, 150])
        self.hcre.set(name='GW strain hxeq',
                      label=r'$\mathcal{D}h_{\times,\mathrm{eq}}$',
                      limits=[-150, 150])
        self.hcrp.set(name='GW strain hxpol',
                      label=r'$\mathcal{D}h_{\times,\mathrm{pol}}$',
                      limits=[-150, 150])

        ## Load the components of the quadrupole moment if present
        self.quadrupole_set = {'txx', 'tyy', 'tzz',
                               'txy', 'txz', 'tyz'}
        self.has_tensor = False
        if kwargs.keys() == self.quadrupole_set:
            self.has_tensor = True
            if strain_units == u.dimensionless_unscaled:
                uu = self.distance.to(u.cm) * u.cm
            else:
                uu = strain_units.to(u.cm) * u.cm
            for name, value in kwargs.items():
                vv = (value * uu)
                vv.set(name=name, label=r'$\ddot{t}_{' + name[1:] + '}$')
                setattr(self, name, vv)
        else:
            self.quadrupole_set = {}

        self.set_fft_config()
        ## Utilities
        self.is_detrended = False ## has the EMD be performed and the residual removed
        self.is_regular = False
        self.is_padded = False
        self.asds = {}
        self.SNR_ene = {}
        self.SNR_fft = {}
    
    def __getitem__(self, indices) -> GWstrain:
        """
        Slice the object
        
        Returns
        -------
        GWstrain
            A brand new object but shorter
        """
        if self.sim_dim == 2:
            outdata = GWstrain(
                np.stack((self.time[indices].value, self.hple[indices].value),
                         axis=-1),
                self.time.unit,
                self.hple.unit,
                self.distance,
                self.tob_corrected,
                self.tob
            )
            outdata.is_detrended = self.is_detrended
            outdata.fft_config = self.fft_config
        else:
            outdata = GWstrain(
                np.stack((self.time[indices].value,
                          self.hple[indices].value,
                          self.hplp[indices].value,
                          self.hcre[indices].value,
                          self.hcrp[indices].value),
                          axis=-1),
                self.time.unit,
                self.hple.unit,
                self.distance,
                self.tob_corrected,
                self.tob,
                **{
                    name: getattr(self, name)[indices] for name in
                    self.quadrupole_set
                }
            )
            outdata.is_detrended = self.is_detrended
            outdata.fft_config = self.fft_config
        return outdata

    def __dimensionless(self, key: str, dimensionless: bool) -> aerray:
        """
        Returns the strain as aerray.

        Parameters
        ----------
        key : str
            key to the strain to return
        dimensionless : bool
            if the strain has to be returned as dimensionless

        Returns
        -------
        aerray
            the strain either dimensionless or not
        """
        strain = getattr(self, key).copy()
        if dimensionless:
            lb, nm, lm = strain.label, strain.name, strain.limits
            lb = lb.replace(r'\mathcal{D}', '')
            lm = [(lm[0] * strain.unit / self.distance).to(u.dimensionless_unscaled).value,
                  (lm[1] * strain.unit / self.distance).to(u.dimensionless_unscaled).value]
            outstrain = (strain / self.distance).to(u.dimensionless_unscaled)
            outstrain.set(label=lb, name=nm, limits=lm)
            return outstrain
        return strain   
    
    def __regularise(self) -> None:
        """
        Interpolate the time and the GW strain on to a regularly spaced
        grid. If both dt and n are None the total number of points will
        be used. New attributes called "attr_ref" will be created.
        """
        if self.fft_config['dt'] is not None:
            n = int(((self.time[-1] - self.time[0]) / 
                     self.fft_config['dt']).value)
        elif self.fft_config['n'] is not None:
            n = self.fft_config['n']
        else:
            n = len(self.time)
        new_time = np.linspace(self.time[0].value, self.time[-1].value, n,
                        endpoint=True)
        self.time_ref = aerray(
            new_time,
            self.time.unit, self.time.name, self.time.label,
            limits=self.time.limits
        )
        for hh in ['hple_ref', 'hplp_ref', 'hcre_ref', 'hcrp_ref']:
            h = getattr(self, hh)
            ## Save the stuff
            lb, nm, lm = h.label, h.name, h.limits
            h = np.interp(new_time, self.time.value, h.value) * h.unit
            h.set(name=nm, label=lb, limits=lm)
            setattr(self, hh, h)
        
        if self.has_tensor:
            for tt in self.quadrupole_set:
                t = getattr(self, tt).copy()
                ## Save the stuff
                lb, nm, lm = t.label, t.name, t.limits
                t = np.interp(new_time, self.time.value, t.value) * t.unit
                t.set(name=nm, label=lb, limits=lm)
                setattr(self, f'{tt}_ref', t)
        self.is_regular = True
    
    def __windowing(self) -> None:
        """
        Applies a window function to all the polarisations.
        """
        if self.fft_config['window_type'] is None:
            win = np.ones
        else:
            try:
                win = getattr(scipy.signal.windows, self.fft_config['window_type'])
            except:
                win = getattr(np, self.fft_config['window_type'])
            win_kwargs = inspect.signature(win).parameters
            window_kwargs = {
                k: v for k, v in self.fft_config['window_kwargs'].items()
                if k in win_kwargs
            }
            try:
                win = win(len(self.time_ref), **window_kwargs)
            except:
                win = win(len(self.time_ref))
        ## renormalise the window to account for the lost power
        wind_norm = np.sqrt(np.sum(win ** 2) / len(win))
        win = win / wind_norm
        
        for hh in ['hple_ref', 'hplp_ref', 'hcre_ref', 'hcrp_ref']:
            h = getattr(self, hh)
            ## Save the stuff
            lb, nm, lm = h.label, h.name, h.limits
            h = win * h
            h.set(name=nm, label=lb, limits=lm)
            setattr(self, hh, h)
        self.is_windowed = True

    def __pad(self) -> None:
        """
        Pad the waveform with the selected value. This will reset the 
        data so proceed with a copy if you do not want to mess everything
        up.
        """
        length = self.fft_config['pad_length'].to(self.time.unit)
        if hasattr(self, 'time_ref'):
            time = self.time_ref
        else:
            time = self.time.copy()
        dt = time[1] - time[0]
        n = int((length / (dt)).value)
        time_before = np.linspace((time[0] - length).value,
                                  time[0].value, n,
                                  endpoint=False) * time.unit
        time_after = np.linspace((time[-1] + dt).value, 
                                 (time[-1] + dt + length).value, n,
                                 endpoint=True) * time.unit
        ## store the time labels
        lb, nm = time.label, time.name
        self.time_ref = np.concatenate(
            (time_before, time, time_after)
        )
        self.time.set(name=nm, label=lb, limits=[-0.005, self.time[-1].value])
        
        if isinstance(self.fft_config['pad_value'], list):
            padding_left = np.ones(n) * self.fft_config['pad_value'][0] * \
                self.hple.unit
            padding_right = np.ones(n) * self.fft_config['pad_value'][1] * \
                self.hple.unit
        else:
            padding_left = np.ones(n) * self.fft_config['pad_value'] * \
                self.hple.unit
            padding_right = np.ones(n) * self.fft_config['pad_value'] * \
                self.hple.unit

        for hh in ['hple_ref', 'hplp_ref', 'hcre_ref', 'hcrp_ref']:
            h = getattr(self, hh)
            ## Save the stuff
            lb, nm, lm = h.label, h.name, h.limits
            h = np.concatenate((
                padding_left, h, padding_right
            ))
            h.set(name=nm, label=lb, limits=lm)
            setattr(self, hh, h)

        self.is_padded = True
    
    def __compute_spectra(self) -> None:
        """
        Computes the spectra of the amplitudes of the strain with the
        option selected by the fft configuration
        """
        self.time_ref = self.time.copy()
        for hh in ['hple', 'hplp', 'hcre', 'hcrp']:
            setattr(self, f'{hh}_ref', getattr(self, hh).copy())
        if self.fft_config['regularise']:
            self.__regularise()
        if self.fft_config['apply_window']:
            self.__windowing()
        if self.fft_config['pad']:
            self.__pad()
        ## Compute the frequency
        dt = self.time_ref[1] - self.time_ref[0]
        self.frequency = aerray(
            rfftfreq(len(self.time_ref), dt.value),
            u.Hz,
            'frequency',
            r'$f$',
            limits=[1e1, 1e4],
            log=True
        )
        ## Compute the real fourier transform of the signal
        for hh in ['hple_ref', 'hplp_ref', 'hcre_ref', 'hcrp_ref']:
            h = getattr(self, hh)
            fft = rfft(getattr(self, hh)) * dt
            fft.set(name=f'FFT {h.name}', label=merge_strings(
                r'$\mathcal{D}$', apply_symbol(h.label.replace(r'\mathcal{D}',
                                                               ''))
            ))
            setattr(self, f'{hh[:4]}_fft', fft)
        
    def __recompute_spectra(self) -> None:
        """
        Recomputes the spectra of the strains and all the data using the
        spectra
        """
        self.__compute_spectra()
        if hasattr(self, 'dEdf_eq'):
            self.compute_dEdf()
        if hasattr(self, 'hchar_eq_ene'):
            self.compute_characteristic_strain()
        if len(self.SNR_fft) != 0:
            self.compute_SNR()
    
    def copy(self) -> GWstrain:
        """
        Creates a copy of this object

        Returns
        -------
        GWstrain
            Return a GWStrain object
        """
        return copy.deepcopy(self)
    
    def set_fft_config(self,
                       regularise: bool = True,
                       dt: aerray = None,
                       n: int = None,                       
                       pad: bool = False,
                       pad_value: float | list = 0,
                       pad_length: aerray=(1*u.s),
                       apply_window: bool = True,
                       window_type: str = 'hanning',
                       **window_kwargs
                       ) -> None:
        """
        Sets the configuration for the fast fast fourier transform in a
        dictionary. If the configuration has changed also recomputes the
        spectra. And all quantities copmputed with those.
        
        Parameters
        ----------
        regularise : bool, optional
            resamples the data on to a regular grid with linear 
            interpolation, by default True
        dt : aerray, optional
            spacing between to consecutive time values, by default None
        n : int, optional
            number of intervals of the time array, by default None
        pad : bool, optional
            pads the signal with constant values
        pad_value : float | list, optional
            The value to place before and after, by default 0
        pad_length : aerray, optional
            how much time to pad, by default (1*u.s)
        window : bool, optional
            applies a window on to the signal, by default True
        window_type: str, optional
            the type of window to consider, by default 'hann',
        **window_kwargs:
            optional parameters for the window type
        """
        self.fft_config = {
            'regularise': regularise,
            'dt': dt,
            'n': n,
            'pad': pad,
            'pad_value': pad_value,
            'pad_length': pad_length,
            'apply_window': apply_window,
            'window_type': window_type,
            'window_kwargs': window_kwargs             
        }
        
        if hasattr(self, 'frequency'):
            self.__recompute_spectra()
    
    def load_asd(self, asd:np.ndarray|aerray|aeseries,
                 frequency:np.ndarray|aerray=None,
                 name:str=None) -> None:
        """
        save a psd in the dictionary

        Parameters
        ----------
        asd : np.ndarray | aerray | aeseries
            The main asd. If it is not an aeserise, it needs the frequency
        frequency : np.ndarray | aerray, optional
            frequency of at which the asd is computed, by default None
        name : str, optional
            name of the detector, by default None
        """
        if isinstance(asd, aeseries):
            self.asds[asd.data.name] = asd
        else:
            if frequency is None and name is None:
                raise KeyError("Frequency and name cannot be None if "\
                    "the asd is different than aeseries.")
            if isinstance(frequency, np.ndarray):
                frequency = aerray(frequency, u.Hz, name, name)
            if isinstance(asd, np.ndarray):
                asd = aerray(asd, (u.Hz**(-0.5)), name, name)
            self.asds[name] = aeseries(asd,
                                       frequency=frequency)
    
    def get_polarisation(self, comp:Literal['h+eq', 'h+pol', 'hxeq',
                                           'hxpol', 'all']='h+eq',
                         dimensionless: bool = False,
                         refined: bool = False) -> aeseries:
        """
        Returns the given polarisation as a aeseries, needed for
        compatibility with AeViz

        Parameters
        ----------
        comp : str, optional
            the strain to return, possibilities ['h+eq', 'h+pol',
            'hxeq', 'hxpol', 'all'], by default 'h+eq'
        dimensionless: bool, optional
            if True returns the selected polarisation at a given distance
        refined: bool, optional
            if the refined polarisation has to be returned

        Returns
        -------
        aeseries
            the chosesn polarisation or list of as a aeseries.
            It consists of time plus polarisation
        """
        ref = '_ref' if refined else ''
        if comp == 'all':
            if self.sim_dim == 2:
                return aeseries(
                    self.__dimensionless(f'hple{ref}', dimensionless),
                    time=getattr(self, f'time{ref}').copy()
                )
            else:
                return [
                    aeseries(
                        self.__dimensionless(f'hple{ref}', dimensionless),
                        time=getattr(self, f'time{ref}').copy()
                    ),
                    aeseries(
                        self.__dimensionless(f'hplp{ref}', dimensionless),
                        time=getattr(self, f'time{ref}').copy()
                    ),
                    aeseries(
                        self.__dimensionless(f'hcre{ref}', dimensionless),
                        time=getattr(self, f'time{ref}').copy()
                    ),
                    aeseries(
                        self.__dimensionless(f'hcrp{ref}', dimensionless),
                        time=getattr(self, f'time{ref}').copy()
                    )                    
                ]
        elif comp == 'h+eq':
            return aeseries(
                self.__dimensionless(f'hple{ref}', dimensionless),
                time=getattr(self, f'time{ref}').copy()
                )
        elif comp == 'h+pol':
            return aeseries(
                self.__dimensionless(f'hplp{ref}', dimensionless),
                time=getattr(self, f'time{ref}').copy()
                )
        elif comp == 'hxeq':
            return aeseries(
                self.__dimensionless(f'hcre{ref}', dimensionless),
                time=getattr(self, f'time{ref}').copy()
                )
        elif comp == 'hxpol':
            return aeseries(
                self.__dimensionless(f'hcrp{ref}', dimensionless),
                time=getattr(self, f'time{ref}').copy()
                )
    
    def set_distance(self, distance: aerray) -> None:
        """
        Set the distance at which to get the dimensionless strain

        Parameters
        ----------
        distance : aerray
        """
        self.distance = distance
    
    def detrend(self, mode:Literal['emd', 'analytical'] = 'emd') -> None:
        """
        Removes a possible trend in each polarisation. This will
        overwrite the polarisations.
        If you need the original please
        make a copy.
        
        Parameters
        ----------
        mode : Literal['emd', 'analytical'], optional
               detrend method to use. emd performs and empirical mode
               decomposition (emd) and removes the residual of the
               procedure from the GW signal. This works because the 
               memory signal is non monotonic so it is not oscillatory
               and remains as the residual.
               If analytical is selectred it employs the analytical
               solution found by Richardson et al 2024
               https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.133.231401
               Between the two, if possible emd seems more reliable, and
               less model dependent.
               
               .. math:: 
                   h_{\rm fit} = \frac{L}{1+\exp{-k(t-t_0)}}[1-H(t-t_s)] + 
                   \frac{L}{2}(1+\cos[2\pi f_t(t-t_s)])
        """
        try:
            from AeViz.utils.physics.EMD_utils import remove_residuals
            pyemd = True
        except:
            pyemd = False
            pass
        def butter_filter(signal, fs, cutoff, order=4):
            nyquist = 0.5 * fs
            normal_cutoff = np.array(cutoff) / nyquist
            b, a = scipy.signal.butter(order, normal_cutoff, btype='low',
                                       analog=False)
            return scipy.signal.filtfilt(b, a, signal.value)
        
        def logistic_function(t, L, k, t0):
            return L / (1.0 + np.exp(-k * (t - t0)))
        
        if not self.is_detrended:
            if pyemd and mode == 'emd':
                for hh in ['hple', 'hplp', 'hcre', 'hcrp']:
                    h = getattr(self, hh)
                    h = remove_residuals(h, self.time)
                    setattr(self, hh, h)
            else:
                dt = self.time[1] - self.time[0]
                fs = (1/dt).to(u.Hz).value
                ts = self.time[-1]
                window = int(((10*u.ms) / dt).to(u.dimensionless_unscaled).value)
                polyorder = 3 if window > 3 else 1
                bi = np.argmax(self.time >= 0)
                for hh in ['hple', 'hplp', 'hcre', 'hcrp']:
                    h = getattr(self, hh)
                    ## Remove high frequencies from the strain. 50 Hz 
                    ## is a pretty safe bet, since then the signal should
                    ## change that much over time
                    h_filt = butter_filter(h, fs, 50) * h.unit
                    ## smooth the signal
                    h_filt = scipy.signal.savgol_filter(h_filt.value,
                                                        window,
                                                        polyorder) * h_filt.unit
                    ## Compute the coefficients
                    L = h_filt[-1]
                    ## find t0 as where h = L/2
                    t0 = self.time[np.argmin(np.abs(h_filt-L/2))]
                    ## Find the first intersection with zeros beginning
                    ## from the end
                    try:
                        itr = np.where(h_filt.value[1:] * h_filt.value[:-1] < 0 )[0][-1] + 1
                    except:
                        continue
                    tr = ts - self.time[itr]
                    k = 1 / tr * np.log(81)
                    ft = 0.1 * u.Hz
                    ## Since we do are not interested in the extrapolation
                    ## part, we just consider the first term
                    exponent = k * (self.time - t0)
                    ## Find the heuristic fit
                    h_fit_heu = \
                        L / (1 + np.exp(-(exponent).to(u.dimensionless_unscaled).value))
                    ## Fit with scipy
                    try:
                        popt, pcov = curve_fit(logistic_function,
                                               self.time.value.to(u.s),
                                               h_filt.value,
                                               p0 = [L.value, k.value, t0.value])
                        L, k, t0 = popt
                        h_fit = logistic_function(self.time.value, L, k, t0) * h.unit
                        setattr(self, hh, h - h_fit)
                    except:
                        setattr(self, hh, h - h_fit_heu)
            self.is_detrended = True
    
    def compute_luminosity(self) -> None:
        """
        Computes the GW luminosity and energy for 2 and 3D simulations.
        In 2D simulations we have,
        from https://iopscience.iop.org/article/10.1086/379822:
        .. math:: 
            L_{GW} = \frac{c^3}{32\pi G} \frac{15}{64\pi} \frac{dh}{dt}^2
        
        in 3D from Scheidegger et al 2008 
        https://www.aanda.org/articles/aa/abs/2008/40/aa8577-07/aa8577-07.html
        .. math:: 
            L_{GW} = \frac{2c^3}/{15G} (\dot I_{xx}^2 + \dot I_{yy}^2
            \dot I_{zz}^2 - \dot I_{xx}\dot I_{yy} - \dot I_{xx}\dot I_{zz} -
            \dot I_{yy}\dot I_{zz} +
            3 (\dot I_{xy}^2+\dot I_{xz}^2+\dot I_{yz}^2))
        """
        dt = self.time[1] - self.time[0]
        if self.sim_dim == 2:
            const = c.c ** 3 / (32 * np.pi * c.G) * 15 / (64 * np.pi)
            self.LGW = (const * 
                        IDL_derivative(self.time, self.hple) ** 2).to(u.erg/u.s)
        else:
            const = 2 * c.c ** 3 / 15 / u.G
            self.LGW = const * (IDL_derivative(self.time, self.txx) ** 2 +
                                IDL_derivative(self.time, self.tyy) ** 2 +
                                IDL_derivative(self.time, self.tzz) ** 2 -
                                (IDL_derivative(self.time, self.txx) * IDL_derivative(self.time, self.tyy) +
                                 IDL_derivative(self.time, self.txx) * IDL_derivative(self.time, self.tzz) +
                                 IDL_derivative(self.time, self.tyy) * IDL_derivative(self.time, self.tzz)) +
                                3 * (IDL_derivative(self.time, self.txy) ** 2 + 
                                     IDL_derivative(self.time, self.txz) ** 2 + 
                                     IDL_derivative(self.time, self.tyz) ** 2)
                                )
        
        self.LGW.set(name='GW luminosity', label=r'$L_{\rm GW}$')
        self.EGW = np.nancumsum(self.LGW * dt).to(u.erg)
        self.EGW.set(name='GW energy', label=r'$E_{\rm GW}$')

    def get_luminosity(self) -> aeseries:
        """
        Extract the GW luminosity.
        
        Returns
        -------
        aeseries
            Contains the luminosity and the evolution time
        
        """
        if not hasattr(self, 'LGW'):
            self.compute_luminosity()
        return aeseries(
            self.LGW,
            time = self.time
        )
        
    def get_energy(self) -> aeseries:
        """
        Extract the evolution of the GW energy over time.
        
        Returns
        -------
        aeseries
            Contains the energy and the evolution time
        
        """
        if not hasattr(self, 'LGW'):
            self.compute_luminosity()
        return aeseries(
            self.EGW,
            time = self.time
        )
    
    def compute_dEdf(self) -> None:
        r"""
        Computes the energy spectrum derivatives with respect to the
        frequency. This will store the frequency and the spectrum in two
        attributes. So if you regularise the strain please rerun this.
        We follow eq (45) in Kuroda et al 2014, 
        https://doi.org/10.1103/PhysRevD.89.044011
        .. math::
            \frac{dE}{df} = \frac{\pi}{4}\frac{c^3}{G}f^2
            (|\tilde A_+|^2+|\tilde A_\times|^2)
        """
        if not  hasattr(self, 'frequency'):
            self.__compute_spectra()
        ## Compute h+ and hx
        A_eq = np.abs(self.hple_fft) ** 2 + np.abs(self.hcre_fft) ** 2
        A_pl = np.abs(self.hplp_fft) ** 2 + np.abs(self.hcrp_fft) ** 2
        const = np.pi / 4 * c.c ** 3 / c.G * self.frequency ** 2
        self.dEdf_eq = (A_eq * const).to(u.erg/u.Hz)
        self.dEdf_pol = (A_pl * const).to(u.erg/u.Hz)
        self.dEdf_eq.set(name='energy frequency eq',
                         label=r'$\frac{\mathrm{d}E_\mathrm{eq}}{\mathrm{d}f}$')
        self.dEdf_pol.set(name='energy frequency eq',
                         label=r'$\frac{\mathrm{d}E_\mathrm{pol}}{\mathrm{d}f}$')

    def get_dEdf(self, los:Literal['eq', 'pol']='eq') -> aeseries:
        """
        Extract the energy spectrum for each frequency for the selected
        line of sight (equatorial or polar)

        Parameters
        ----------
        los : Literal['eq', 'pol'], optional
            line of sight along which to compute the energy spectrum from,
            by default 'eq'

        Returns
        -------
        aeseries
            contains the frequency and the energy spectrum
        """
        try:
            dEdf = getattr(self, f'dEdf_{los}')
        except:
            self.compute_dEdf()
            dEdf = getattr(self, f'dEdf_{los}')
        return aeseries(
            dEdf,
            frequency=self.frequency
        )
        
    def compute_characteristic_strain(self) -> None:
        r"""
        Computes the characteristic strain for the two observers.
        It is computed in two ways:
        From the energy spectrum following Kuroda et al 2014, 
        https://doi.org/10.1103/PhysRevD.89.044011:
        ..math::
            h_{\rm char} = \sqrt{\frac{2}{\pi^2}\frac{G}{c^3}\frac{1}{D^2}
            \frac{\mathrm{d}E}{\mathrm{d}f}}
        
        From the fourier transform following Moore et al 2015,
        https://doi.org/10.1103/PhysRevD.89.044011:
        ..math:: 
            h_{\rm char} = 2 f^2 |\tilde h(f)|
        """
        if not hasattr(self, 'frequency'):
            self.compute_dEdf()
        const = 2 / np.pi ** 2 * c.G / c.c ** 3 / self.distance ** 2
        self.hchar_eq_ene = np.sqrt((const * self.dEdf_eq)).to(u.dimensionless_unscaled)
        self.hchar_eq_ene.set(name='hchar pol', label=r'$h_{\rm char, eq}$')
        self.hchar_pol_ene = np.sqrt((const * self.dEdf_pol)).to(u.dimensionless_unscaled)
        self.hchar_pol_ene.set(name='hchar pol', label=r'$h_{\rm char, pol}$')
        dt = (self.time[1]- self.time[0]).to(u.s)
        
        self.hchar_eq_fft = (2 * self.frequency * np.sqrt(
            (np.abs(self.hple_fft) / self.distance) ** 2 + 
            (np.abs(self.hcre_fft) / self.distance) ** 2)).to(u.dimensionless_unscaled)
        self.hchar_eq_fft.set(name='hchar pol', label=r'$h_{\rm char, pol}$')
        self.hchar_pol_fft = (2 * self.frequency * np.sqrt(
            np.abs(self.hplp_fft / self.distance) ** 2 + 
            np.abs(self.hcrp_fft / self.distance) ** 2)).to(u.dimensionless_unscaled)
        self.hchar_pol_fft.set(name='hchar pol', label=r'$h_{\rm char, pol}$')
    
    def get_characteristic_strain(self, los:Literal['eq', 'pol']='eq',
                                  mode:Literal['energy', 'fft']='fft',
                                  divided_by_frequency:bool=True) -> aeseries:
        """
        Extract the characteristic strain for the selected line of sight.

        Parameters
        ----------
        los : Literal['eq', 'pol'], optional
            line of sight, equatorial or polar, by default 'eq'
        mode : Literal['energy', 'fft'], optional
            computation of the characteristic strain with the
            energy spectrum or the fft, by default 'fft'
        divided_by_frequency : bool
            divides the characteristc strain by the square root of the 
            frequency, useful for comparing with the detector psd,
            by default True

        Returns
        -------
        aeseries
            contains the frequency and the characteristic strain
        """
        if not hasattr(self, 'hchar_eq_ene'):
            self.compute_characteristic_strain()
        # constrauct the name
        mode = 'ene' if mode == 'energy' else mode
        nm = f'hchar_{los}_{mode}'
        hchar = getattr(self, nm).copy()
        if divided_by_frequency:
            nm, lb = hchar.name, hchar.label
            lb = merge_strings(lb, r'$/\sqrt{f}$')
            hchar /= np.sqrt(self.frequency)
            hchar.set(name=nm, label=lb, log=True)
        return aeseries(
            hchar,
            frequency = self.frequency
        )
    
    def compute_SNR(self, detector:str=None) -> None:
        r"""
        Computes the SNR for the given detector asd.
        The SNR is computed as
        ..math:: 
            \mathcal{\rho} = \sqrt{\int\mathrm{d}\ln f \frac{h_{\rm char}^2}{f ASD^2}}
        Parameters
        ----------
        detector : str, optional
            name of the detector, if None the SNR will be computed for 
            all the detectors, by default None
        """
        def interpolate_ASD(ASD, freq):
            new_ASD =  np.interp(freq.value, ASD.frequency.value, ASD.data.value,
                                 left = np.nan, right = np.nan)
            new_ASD = new_ASD.data
            new_ASD = aerray(new_ASD, (u.Hz**-0.5), ASD.data.name,
                             r'$ASD_\mathrm{' + ASD.data.name + '}$')
            return new_ASD
        
        if not self.is_detrended:
            warnings.warn("The strain is not detrended, "\
                "this may cause the SNR to be larger.")
        if not hasattr(self, 'hchar_eq_ene'):
            self.compute_characteristic_strain()
        df = (self.frequency[1] - self.frequency[0]) / self.frequency ** 2
        if detector is not None:
            if not detector in self.asds.keys():
                raise ValueError("Detector not recognised, must be in"\
                    f"{self.asds.keys()}")
            newASD = interpolate_ASD(self.asds[detector],
                                     self.frequency)
            self.SNR_ene[detector] = {
                'eq': np.sqrt(np.nancumsum(df * self.hchar_eq_ene ** 2 / newASD ** 2)),
                'pol': np.sqrt(np.nancumsum(df * self.hchar_pol_ene ** 2 / newASD ** 2))
            }
            self.SNR_fft[detector] = {
                'eq': np.sqrt(np.nancumsum(df * self.hchar_eq_fft ** 2 / newASD ** 2)),
                'pol': np.sqrt(np.nancumsum(df * self.hchar_pol_fft ** 2 / newASD ** 2))
            }
        else:
            for det in self.asds.keys():
                newASD = interpolate_ASD(self.asds[det],
                                     self.frequency)
                self.SNR_ene[det] = {
                    'eq': np.sqrt(np.nancumsum(df * self.hchar_eq_ene ** 2 / newASD ** 2)),
                    'pol': np.sqrt(np.nancumsum(df * self.hchar_pol_ene ** 2 / newASD ** 2))
                }
                self.SNR_fft[det] = {
                    'eq': np.sqrt(np.nancumsum(df * self.hchar_eq_fft ** 2 / newASD ** 2)),
                    'pol': np.sqrt(np.nancumsum(df * self.hchar_pol_fft ** 2 / newASD ** 2))
                }
    
    def get_SNR(self, detector:str, los:Literal['eq', 'pol'],
                mode:Literal['energy', 'fft'],
                evolution: bool = False) -> float | aeseries:
        """
        Returrns the value of the SNR for the distance and detector selected.

        Parameters
        ----------
        detector : str
            name of the detector
        los : Literal['eq', 'pol']
            line of sight along which compute the SNR
        mode : Literal['energy', 'fft']
            the mode of computing the characteristic strain
        evolution : bool, optional
            returns the frequency evolution of the SNR
        Returns
        -------
        float | value
            value of the SNR or the aeseries  containing the frequency
            evolution, and the frequency
        """
        mode = 'ene' if mode == 'energy' else mode
        if not detector in getattr(self, f'SNR_{mode}'):
            self.compute_SNR(detector)
        SNR = getattr(self, f'SNR_{mode}')[detector][los]
        if evolution:
            return aeseries(SNR,
                            frequency = self.frequency)
        return SNR[-1].value
    
    def compute_modes(self) -> None:
        """
        Computes the right and left-handed modes for the GW signal for 
        both the equatorial and polar lines of sight.
        The right-handed mode is defined as
        .. math::
            h_R=(h_+ - ih_\times)\sqrt{2}
        While the left-handed mode is:
        .. math::
            h_L=(h_+ + ih_\times)\sqrt{2}
        """
        self.hL_eq = (self.hple.value + self.hcre.unit * 1j) / np.sqrt(2) * self.hple.unit
        self.hL_eq.set(name='left-handed mode equatorial', label=r'$\mathcal{D}h_{L,eq}$')
        self.hR_eq = (self.hple.value - self.hcre.unit * 1j) / np.sqrt(2) * self.hple.unit
        self.hL_eq.set(name='right-handed mode equatorial', label=r'$\mathcal{D}h_{R,eq}$')
        
        self.hL_pol = (self.hplp.value + self.hcrp.unit * 1j) / np.sqrt(2) * self.hplp.unit
        self.hL_pol.set(name='left-handed mode polar', label=r'$\mathcal{D}h_{L,pol}$')
        self.hR_pol = (self.hplp.value - self.hcrp.unit * 1j) / np.sqrt(2) * self.hplp.unit
        self.hL_pol.set(name='right-handed mode polar', label=r'$\mathcal{D}h_{R,pol}$')
    
    def get_modes(self, mode:Literal['left', 'right'],
                  los:Literal['eq', 'pol'],
                  dimensionless=False) -> aeseries:
        """
        Extract the right and left handed mode of the GW polarisation.

        Parameters
        ----------
        mode : Literal['left', 'right']
            which mode to extract, left or right handed
        los : Literal['eq', 'pol']
            line of sight ot consider
        dimensionless : bool, optional
            if to return it divided by the distance, by default False
            

        Returns
        -------
        aeseries
            series containing the time and mode
        """

        name = f'h{mode[0].capitalize()}_{los}'
        if not hasattr(self, name):
            self.compute_modes()
        hh = getattr(self, name).copy()
        if dimensionless:
            nm, lb = hh.name, hh.label
            hh = (hh / self.distance).to(u.dimensionless_unscaled)
            hh.set(name=nm, label=lb.replace(r'\mathcal{D}', ''))
        return aeseries(
            hh,
            time = self.time.copy()
            )
  
    def get_stokes_parameters(self,
                              parameter:Literal['V', 'I', 'Q', 'U'],
                              los: Literal['eq', 'pol'],
                              time_range: list[aerray|float] | None = None,
                              spectrogram: bool=False,
                              window_size: aerray = (10*u.ms),
                              window: str = None,
                              overlap: float = 0.5,
                               **kwargs) -> aeseries:
        """
        Computes and returns the stokes parameters V, I, Q, and U. These
        are not computed from the hR and hL modes, rather from the
        strain.

        Parameters
        ----------
        parameter : Literal['V', 'I', 'Q', 'U']
            The stokes parameter to compute
        los : Literal['eq', 'pol']
            line of sight to consider
        time_range : list[aerray | float] | None, optional
            time range to considered to find the peak. If None the full
            strain is analysed, by default None
        spectrogram : bool, optional
            if it is true it computes a short time fourier transform in 
            place of a FFT, by default False
        window_size : aerray, optional
            the size of the window used to compute the STFT,
            by default 10 * u.ms
        overlap : float, optional
            the overlap between two subsequent STFT, by default 0.5
        window : str
            The name of the window, any supported by
            scipy.signal.windows, by default None
        **kwargs :
            Anything supported by the specific window

        Returns
        -------
        aeseries
            contains the stokes parameter for the selected line of sight,
            the frequency and the time if one selects the spectrogram
        """
        hp = getattr(self, f'hpl{los[0]}').copy()
        hc = getattr(self, f'hcr{los[0]}').copy()
        if time_range is not None:
            ini = np.argmax(self.time>=time_range[0])
            isto = np.argmax(self.time>=time_range[1]) + 1
            tm = self.time[ini:isto]
            hp = hp[ini:isto]
            hc = hc[ini:isto]
            
        dt = (tm[1] - tm[0]).to(u.s)
        fs = (1 / dt).to(u.Hz)
        ## Get the window
        if window is None:
            win = np.ones
        else:
            try:
                win = getattr(scipy.signal.windows, window)
            except:
                win = getattr(np, window)
        if spectrogram:
            win_len = int((window_size / dt).to(u.dimensionless_unscaled).value)
            hop = overlap * win_len
            try:
                win = win(hop, **kwargs)
            except:
                win = win(hop)
            
            SFT = scipy.signal.ShortTimeFFT(win, hop, fs,
                                            scale_to='magnitude')
            frequency = aerray(SFT.f, u.Hz, 'frequency', r'$f$', 
                               None, [0, 2000], False)
            time = aerray(SFT.t(len(tm)) + tm[0].value,
                          u.s, tm.name, tm.label, None,
                          tm.limits)
            hptilde = SFT.stft(hp.value)
            hctilde = SFT.stft(hc.value)
            un = hp.unit
        else:
            try:
                win = win(len(tm), **kwargs)
            except:
                win = win(len(tm))
            frequency = aerray(
                rfftfreq(len(tm), dt.value), u.Hz, "frequency",
                r'$f$', None, [1, 2000], True
            )
            hptilde = rfft(hp.value * win) * dt.value
            hctilde = rfft(hc.value * win) * dt.value
            un = hc.unit * u.s
        ## compute the selected stokes parameter
        if parameter == 'V':
            par = 1j * (hptilde * np.conjugate(hctilde) - hctilde * np.conjugate(hptilde)) / 2
            nm, lb = 'V', r'$V$'
        elif parameter == 'I':
            par = (hptilde * np.conjugate(hptilde) + hctilde * np.conjugate(hctilde)) / 2
            nm, lb = 'I', r'$I$'
        elif parameter == 'Q':
            par = (hptilde * np.conjugate(hptilde) - hctilde * np.conjugate(hctilde)) / 2
            nm, lb = 'Q', r'$Q$'
        elif parameter == 'U':
            par = (hptilde * np.conjugate(hctilde) + hctilde * np.conjugate(hptilde)) / 2
            nm, lb = 'U', r'$U$'
        
        par = aerray(par, un, nm, lb, 'jet', [par.min(), par.max()])
        if spectrogram:
            par.set(limits=[par.limits[0] * 1.1, par.limits[1] * 0.9])
            return aeseries(
                par,
                time=time,
                frequency=frequency
            )
        else:
            return aeseries(
                par,
                frequency=frequency
            )

    @overload
    def range_strain_width(self, peak: Literal['bounce', 'max'],
                           comp: Literal['h+eq', 'h+pol', 'hxeq', 'hxpol'],
                           time_range: list[aerray|float] | None,
                           tol: float,
                           dimensionless: bool,
                           return_indices: Literal[True, False]=False
                           ) -> aerray: ...
    
    @overload
    def range_strain_width(self, peak: Literal['bounce', 'max'],
                           comp: Literal['h+eq', 'h+pol', 'hxeq', 'hxpol'],
                           time_range: list[aerray|float] | None,
                           tol: float, 
                           dimensionless: bool,
                           return_indices: Literal[True, False]
                           ) -> tuple[aerray, list[int]]: ...

    def range_strain_width(self, peak: Literal['bounce', 'max']='bounce',
                           comp: Literal['h+eq', 'h+pol', 'hxeq', 'hxpol']='h+eq',
                           time_range: list[aerray|float] | None = None,
                           tol: float = 0.05,
                           dimensionless: bool = False,
                           return_indices: Literal[True, False]=False
                           ):
        """
        Computes the strain range width of the selected GW polarisation.
        It measures the difference between the a cosecutive velley and
        crests operatively defined as:
        last point: intersection point with the x after the peak
        first point: third to last intersection with the x before the crest

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
        dimensionless : bool, optional
            If the value has to be returned as a dimensionless strain,
            by default False
        return_indices : Literal[False], optional
            If the indices have to be returned, by default False

        Returns
        -------
        aerray | tuple[aerray, list[int]]
            returns the value of the range strain width as an aerray. If
            the flag return_indices is set to True returns also the indices
            of the intersections
        """
        
        ## Get the strain
        h = self.get_polarisation(comp=comp)
        if peak == 'bounce':
            ## Find the intersections with the x axis
            zeros = np.where(h.data.value[:-1] * h.data.value[1:] < 0)[0] + 1
            ## Find bounce time
            bouncei = np.argmax(h.time >= (0 * u.s))
            ## Find the zero crossing after bounce
            zero_afterbi = np.argmax(zeros >= bouncei)
            ## Find the valley between two successive zero crossings
            mini = np.argmin(
                h.data[zeros[zero_afterbi]:zeros[zero_afterbi+1] + 1]) + \
                    zeros[zero_afterbi]
            ## Find the maximum between the minimum zero crossing and
            ## and the subsequent one. It needs to be above the absolute
            ## tollerance
            abs_tol = np.abs(h.data[mini]) * tol
            for i in range(zero_afterbi + 2, len(zeros)-2, 2):
                maxi = np.argmax(
                    h.data[zeros[zero_afterbi+1]:zeros[i] + 1]
                ) + zeros[zero_afterbi+1]
                if h.data[maxi] > abs_tol:
                    break
            #save the indices
            firsti = 0 if zero_afterbi == 0 else zeros[zero_afterbi-1]
            indices = [firsti, mini, maxi, zeros[np.argmax(zeros>=maxi)]]
        elif peak == 'max':
            ini = 0
            hh = h
            if time_range is not None:
                ini = np.argmax(h.time>=time_range[0])
                h = h[ini:np.argmax(h.time>=time_range[1])]
            ## Find the max index and the zeros
            zeros = np.where(h.data.value[:-1] * h.data.value[1:] < 0)[0] + 1
            maxi = np.argmax(np.abs(h.data))
            ## consider the highest exursion as positive, always
            if h.data[maxi] < 0:
                h.data *= -1
            ## Find the zero intersection after the maximum
            zero_aftermax = np.argmax(zeros >= maxi)
            zeros_premax = zero_aftermax - 1
            ## We check all the points after the maximum that are not
            ## included in the tolerance
            abs_tol = h.data[maxi] * tol
            for i in range(zero_aftermax+1, len(zeros)-1, 2):
                mm_gw = np.max(h.data[zeros[i]:zeros[i+1]])
                lasti = zeros[i]
                if mm_gw > abs_tol:
                    break
            ## Do the same thing for the points before the beginning of the
            ## peak
            for i in range(zeros_premax-1, 1, -2):
                mm_gw = np.max(h.data[zeros[i-1]:zeros[i]])
                firsti = zeros[i]
                if mm_gw > abs_tol:
                    break
            ## Compare the two deltah to see which is the highest and
            ## save the indices accordingly
            if np.abs(np.min(h.data[firsti:zeros[zeros_premax]])) > \
                np.abs(np.min(h.data[zeros[zero_aftermax]:lasti])):
                    indices = [
                        firsti + ini,
                        np.argmax(h.data[firsti:zeros[zeros_premax]]) + \
                            firsti + ini,
                        maxi + ini,
                        zeros[zero_aftermax] + ini
                    ]
            else:
                indices = [
                        zeros[zeros_premax] + ini,
                        maxi + ini,
                        np.argmin(h.data[zeros[zero_aftermax]:lasti]) + \
                            zeros[zero_aftermax] + ini,
                        lasti + ini
                    ]
            h = hh
        ## Compute the delta h
        deltah = np.abs(h.data[indices[1]] - h.data[indices[2]])
        hlab = h.data.label
        if dimensionless:
            deltah = deltah / self.distance
            hlab = hlab.replace(r'\mathcal{D}', '')
        
        deltah.set(name='strain range width', 
                   label=merge_strings(r'$\Delta$', hlab))
        if return_indices:
            return (deltah, indices)
        return deltah

    @overload
    def frequency_of_peak(self, peak: Literal['bounce', 'max'],
                           comp: Literal['h+eq', 'h+pol', 'hxeq', 'hxpol'],
                           time_range: list[aerray|float] | None,
                           tol: float,
                           normalised: bool,
                           return_max: Literal[False, True]=False
                          ) -> aeseries: ...
    
    @overload
    def frequency_of_peak(self, peak: Literal['bounce', 'max'],
                           comp: Literal['h+eq', 'h+pol', 'hxeq', 'hxpol'],
                           time_range: list[aerray|float] | None,
                           tol: float,
                           normalised: bool,
                           return_max: Literal[False, True]
                          ) -> aerray: ...
    
    def frequency_of_peak(self, peak: Literal['bounce', 'max'] = 'bounce',
                           comp: Literal['h+eq', 'h+pol', 'hxeq', 'hxpol'] = 'h+eq',
                           time_range: list[aerray|float] | None = None,
                           tol: float = 0.05,
                           normalised: bool = False,
                           return_max: Literal[False]=False
                          ):
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

        Returns
        -------
        aeseries | aerray
            aseseries of the Fourier transform containing the frequency
            and the absolute value of the spectrum. aerray of the maximum
            frequency of the oscillation
        """
        ## Compute the indices of the  oscillation
        _, ind = self.range_strain_width(peak, comp, time_range, tol,
                                         False, True)
        h = self.get_polarisation(comp, False)
        ## Cut the polarisation
        h = h[ind[0]:ind[-1]]
        ## Pad before and after with 1 second of zeros and resample
        ## everything
        dt = (h.time[1] - h.time[0]).to(u.s)
        n = int((1 * u.s) / dt)
        time = np.concatenate(
            (np.linspace((h.time[0] - (1*u.s)).value,
                         h.time[0].value, n,
                         endpoint=False) * h.time.unit,
             h.time,
             np.linspace((h.time[-1] + dt).value, 
                                 (h.time[-1] + dt + 1*u.s).value, n,
                                 endpoint=True) * h.time.unit)
        )
        h = np.concatenate((
            np.zeros(n) * h.data.unit, h.data, np.zeros(n) * h.data.unit
        ))
        new_time = (np.linspace(time[0].value, time[-1].value, len(time),
                        endpoint=True) * time.unit).to(u.s)
        new_h = np.interp(new_time.value, time.to(u.s).value, h.value) * h.unit
        dt = new_time[1] - new_time[0]
        frequency = rfftfreq(len(new_time), dt.value) * u.Hz
        fft = np.abs(rfft(new_h.value) * dt) * new_h.unit
        if return_max:
            imax = np.nanargmax(fft)
            max_freq = frequency[imax]
            max_freq.set(name='Peak frequency', label=r'$f_{\rm peak}$')
            return max_freq
        else:
            frequency.set(name='frequency', label=r'$f$')
            fft.set(name="fourier transform", label=r'$\tilde h_{peak}$')
            if normalised:
                fft /= np.max(fft)
                fft.set(name="fourier transform", label=r'Normalised amplitude')
            return aeseries(
                fft,
                frequency=frequency
            )        