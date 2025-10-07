from __future__ import annotations
from AeViz.units.aerray import aerray
import numpy as np
from AeViz.units import u
from AeViz.utils.files.string_utils import apply_symbol, merge_strings
import numpy as np
import scipy.signal
import warnings
from typing import Literal
from AeViz.utils.math_utils import IDL_derivative
from scipy.integrate import cumulative_simpson, cumulative_trapezoid

def get_selection_indices(a, b):
    """
    Determines how `b` was derived from `a`.
    Returns:
    - A tuple of slice objects if `b` is a slice.
    - A tuple of (array of indices) if `b` is non-contiguous.
    - False if `b` is not found in `a`.
    """
    a=a.value
    b=b.value
    if not np.shares_memory(a, b):  # `b` must be a view of `a`
        return False

    selection = []
    for axis in range(a.ndim):
        # Get all indices of `b` along this axis
        a_indices = np.arange(a.shape[axis])  # Full range
        b_indices = np.where(np.isin(a.take(indices=a_indices, axis=axis), b))[0]

        if len(b_indices) == 0:
            return False  # b is not in a

        # Check for contiguous slice
        start, stop, step = b_indices[0], b_indices[-1] + 1, None
        diff = np.diff(b_indices)  # Check the step

        if np.all(diff == diff[0]):  # Check if step is constant
            step = diff[0]
            selection.append(slice(start, stop, step))
        else:
            selection.append(b_indices)  # Fancy indexing or irregular selection

    return tuple(selection)


class aeseries:    
    def __init__(self, data, **kwargs):
        assert isinstance(data, aerray), "The main quantity has to " \
            "be an array"
        assert data.ndim >= len(kwargs), f"Too many axes ({len(kwargs)}) " \
            f"for a quantity with {data.ndim} dimensions"
        if len(data.shape) != len(set(data.shape)):
            raise NotImplementedError(f"Two equivalent dimensions in quantity")
        self.data = data
        self.__axis_indices = {}
        self.__indices_axis = {}
        self.__axis_names = []
        self.__axis_units = {}
        
        for name, axis in kwargs.items():
            assert isinstance(axis, aerray), f"{name} must be an aerray"
            #assert axis.ndim == 1, f"Axes must be 1D aerrays, {name} is {axis.ndim}"
            setattr(self, name, axis)
            self.__axis_indices[name] = data.shape.index(len(axis))
            self.__axis_names.append(name)
            self.__indices_axis[data.shape.index(len(axis))] = name
            self.__axis_units[name] = axis.unit

    def __repr__(self):
        outstr = "aeseries(\n"
        outstr += f"\tdata: {repr(self.data)}\n"
        for name in self.__axis_names:
            outstr += f"\t{name}: {repr(getattr(self, name))}\n"
        outstr += ")"
        return outstr

    def __setattr__(self, name, value):
        if name in self.__dict__ and name != 'data':
            indices = get_selection_indices(self.__dict__[name], value)
            if indices:
                self.__dict__[name] = value
                ## Slice also the data
                new_ind = []
                for i in range(self.__axis_indices[name]):
                    new_ind.append(slice(0, self.data.shape[i], 1))
                new_ind.append(indices[0])
                for i in range(self.__axis_indices[name]+1, self.data.ndim):
                    new_ind.append(slice(0, self.data.shape[i], 1))
                self.__dict__['data'] = self.data[tuple(new_ind)]
            else:
                warnings.warn(f"New value of {name} is not a slice of its previous value.")
                if value.shape != self.__dict__[name].shape:
                    raise IndexError(f"Mismatch on shapes between old and new value of {name}") 
                self.__dict__[name] = value
        else:
            self.__dict__[name] = value
    
    def __getitem__(self, indices):
        if isinstance(indices, slice) and len(self.__axis_names) == 1:
            return aeseries(self.data[indices],
                        **{name: getattr(self, name)[indices]
                           for name in self.__axis_names})
        return aeseries(self.data[indices],
                        **{name: getattr(self, name)[indices[self.__axis_indices[name]]]
                           for name in self.__axis_names})
    
    ## Operations handling
    def __add__(self, other):
        """
        Redefine the addition, works between aerrays and aeseries.
        If we sum two aeseries only the main "data" value is summed up
        """
        if isinstance(other, aerray):
            if other.unit.is_equivalent(self.data.unit):
                return aeseries(
                    self.data + other,
                    **{name: getattr(self, name) for name in self.__axis_names}
                )
            elif other.unit in self.__axis_units.values() and \
                len(self.__axis_units) == len(set(self.__axis_units.values())):
                qname = list(self.__axis_units.keys())[list(self.__axis_units.values()).index(other.unit)]
                return aeseries(
                    self.data,
                    **{name: getattr(self, name) + other if name == qname else
                    getattr(self, name) for name in self.__axis_names}
                )  
            elif sum([uu.is_equivalent(other.unit) for uu in self.__axis_units.values()]) == 1:
                qname = list(self.__axis_units.keys())[[uu.is_equivalent(other.unit) for uu in self.__axis_units.values()].index(True)]
                return aeseries(
                    self.data,
                    **{name: getattr(self, name) + other if name == qname else
                    getattr(self, name) for name in self.__axis_names}
                )
            else:
                raise TypeError("Units do not match between aeseries and aerray.")
        elif isinstance(other, aeseries):
            if not all([hasattr(other, name) for name in self.__axis_names]):
                raise TypeError("aeserieses do not have the same attributes")
            if not all([np.all(getattr(self, name) == getattr(other, name))
                        for name in self.__axis_names]):
                raise TypeError("Mismatch in the two aeseries axes.")
            return aeseries(self.data + other.data, **{
                name: getattr(self, name) for name in self.__axis_names
            })
        return NotImplemented
    
    def __radd__(self, other):
        return other.__add__(self)
    
    def __iadd__(self, other):
        """
        In place addition, if summing two aeseries, only the "data" value is sum up
        """
        if isinstance(other, aerray):
            if other.unit.is_equivalent(self.data.unit):
                self.__dict__['data'] =  self.data + other
                return self
            elif other.unit in self.__axis_units.values() and \
                len(self.__axis_units) == len(set(self.__axis_units.values())):
                qname = list(self.__axis_units.keys())[list(self.__axis_units.values()).index(other.unit)]
                self.__dict__[qname] += other
                return self
            elif sum([uu.is_equivalent(other.unit) for uu in self.__axis_units.values()]) == 1:
                qname = list(self.__axis_units.keys())[[uu.is_equivalent(other.unit) for uu in self.__axis_units.values()].index(True)]
                self.__dict__[qname] += other
                return self
        elif isinstance(other, aeseries):
            if not all([hasattr(other, name) for name in self.__axis_names]):
                raise TypeError("aeserieses do not have the same attributes")
            if not all([np.all(getattr(self, name) == getattr(other, name))
                        for name in self.__axis_names]):
                raise TypeError("Mismatch in the two aeseries axes.")
            self.data = self.data + other.data
            return self
            
        raise TypeError("In place addition works only between compatible aeseries or aeseries-aerray with compatible units.")
    
    def __sub__(self, other):
        """
        Same as for addition but with the - sign
        """
        if isinstance(other, aerray):
            if other.unit.is_equivalent(self.data.unit):
                return aeseries(
                    self.data - other,
                    **{name: getattr(self, name) for name in self.__axis_names}
                )
            elif other.unit in self.__axis_units.values() and \
                len(self.__axis_units) == len(set(self.__axis_units.values())):
                qname = list(self.__axis_units.keys())[list(self.__axis_units.values()).index(other.unit)]
                return aeseries(
                    self.data,
                    **{name: getattr(self, name) - other if name == qname else
                    getattr(self, name) for name in self.__axis_names}
                )  
            elif sum([uu.is_equivalent(other.unit) for uu in self.__axis_units.values()]) == 1:
                qname = list(self.__axis_units.keys())[[uu.is_equivalent(other.unit) for uu in self.__axis_units.values()].index(True)]
                return aeseries(
                    self.data,
                    **{name: getattr(self, name) - other if name == qname else
                    getattr(self, name) for name in self.__axis_names}
                )
            else:
                raise TypeError("Units do not match between aeseries and aerray.")
        elif isinstance(other, aeseries):
            if not all([hasattr(other, name) for name in self.__axis_names]):
                raise TypeError("aeserieses do not have the same attributes")
            if not all([np.all(getattr(self, name) == getattr(other, name))
                        for name in self.__axis_names]):
                raise TypeError("Mismatch in the two aeseries axes.")
            return aeseries(self.data - other.data, **{
                name: getattr(self, name) for name in self.__axis_names
            })
        return NotImplemented
    
    def __isub__(self, other):
        """
        Same as for addition but with the - sign
        """
        if isinstance(other, aerray):
            if other.unit.is_equivalent(self.data.unit):
                self.data = self.data - other
                return self    
            elif other.unit in self.__axis_units.values() and \
                len(self.__axis_units) == len(set(self.__axis_units.values())):
                qname = list(self.__axis_units.keys())[list(self.__axis_units.values()).index(other.unit)]
                self.__dict__[qname] = self.__dict__[qname] - other
                return self
            elif sum([uu.is_equivalent(other.unit) for uu in self.__axis_units.values()]) == 1:
                qname = list(self.__axis_units.keys())[[uu.is_equivalent(other.unit) for uu in self.__axis_units.values()].index(True)]
                self.__dict__[qname] = self.__dict__[qname] - other
                return self
            else:
                raise TypeError("Units do not match between aeseries and aerray.")
        elif isinstance(other, aeseries):
            if not all([hasattr(other, name) for name in self.__axis_names]):
                raise TypeError("aeserieses do not have the same attributes")
            if not all([np.all(getattr(self, name) == getattr(other, name))
                        for name in self.__axis_names]):
                raise TypeError("Mismatch in the two aeseries axes.")
            self.data = self.data - other.data
            return self
        return NotImplemented
    
    def __mul__(self, other):
        """
        Multiplication: operates only on the main "data"
        """
        if isinstance(other, aeseries):
            if not all([hasattr(other, name) for name in self.__axis_names]):
                raise TypeError("aeserieses do not have the same attributes")
            if not all([np.all(getattr(self, name) == getattr(other, name))
                        for name in self.__axis_names]):
                raise TypeError("Mismatch in the two aeseries axes.")
            return aeseries(self.data * other.data, **{
                name: getattr(self, name) for name in self.__axis_names
            })
        else:
            return aeseries(self.data * other, **{
                name: getattr(self, name) for name in self.__axis_names
            })
    
    def __rmul__(self, other):
        return self.__mul__(other)
    
    def __imul__(self, other):
        """
        Multiplication: operates only on the main "data"
        """
        if isinstance(other, aeseries):
            if not all([hasattr(other, name) for name in self.__axis_names]):
                raise TypeError("aeserieses do not have the same attributes")
            if not all([np.all(getattr(self, name) == getattr(other, name))
                        for name in self.__axis_names]):
                raise TypeError("Mismatch in the two aeseries axes.")
            self.data = self.data * other.data
            return self
        else:
            self.data = self.data * other
            return self
            
    def __truediv__(self, other):
        """
        Multiplication: operates only on the main "data"
        """
        if isinstance(other, aeseries):
            if not all([hasattr(other, name) for name in self.__axis_names]):
                raise TypeError("aeserieses do not have the same attributes")
            if not all([np.all(getattr(self, name) == getattr(other, name))
                        for name in self.__axis_names]):
                raise TypeError("Mismatch in the two aeseries axes.")
            return aeseries(self.data / other.data, **{
                name: getattr(self, name) for name in self.__axis_names
            })
        else:
            return aeseries(self.data / other, **{
                name: getattr(self, name) for name in self.__axis_names
            })
    
    def __itruediv__(self, other):
        """
        Multiplication: operates only on the main "data"
        """
        if isinstance(other, aeseries):
            if not all([hasattr(other, name) for name in self.__axis_names]):
                raise TypeError("aeserieses do not have the same attributes")
            if not all([np.all(getattr(self, name) == getattr(other, name))
                        for name in self.__axis_names]):
                raise TypeError("Mismatch in the two aeseries axes.")
            self.data = self.data / other.data
            return self
        else:
            self.data = self.data / other
            return self
    
    def __pow__(self, exponent):
        """ power gets applied to the "data" component"""
        if not np.isscalar(exponent):
            raise TypeError("Exponent must be a scalar value.")
        return aeseries(self.data ** exponent, **{
                name: getattr(self, name) for name in self.__axis_names
            })
    
    ## Operator, oh could you help me place this call?
    
    ## We differentiate between series and other stupp
    def __eq__(self, other):
        if isinstance(other, aeseries):
            if not all([hasattr(other, name) for name in self.__axis_names]):
                raise TypeError("aeserieses do not have the same attributes")
            comparison = [np.all(self.data == other.data)]
            comparison2 = [np.all(getattr(self, name) == getattr(other, name))
                        for name in self.__axis_names]
            return all(comparison + comparison2)
        else:
            return self.data == other
    
    def __ne__(self, other):
        if isinstance(other, aeseries):
            if not all([hasattr(other, name) for name in self.__axis_names]):
                raise TypeError("aeserieses do not have the same attributes")
            comparison = [np.all(self.data == other.data)]
            comparison2 = [np.all(getattr(self, name) == getattr(other, name))
                        for name in self.__axis_names]
            return not all(comparison + comparison2)
        else:
            return self.data != other
    
    def __ge__(self, other):
        if isinstance(other, aeseries):
            if not all([hasattr(other, name) for name in self.__axis_names]):
                raise TypeError("aeserieses do not have the same attributes")
            if not all([np.all(getattr(self, name) == getattr(other, name))
                        for name in self.__axis_names]):
                raise TypeError("Mismatch in the two aeseries axes.")
            return self.data >= other.data
        else:
            return self.data >= other
    
    def __le__(self, other):
        if isinstance(other, aeseries):
            if not all([hasattr(other, name) for name in self.__axis_names]):
                raise TypeError("aeserieses do not have the same attributes")
            if not all([np.all(getattr(self, name) == getattr(other, name))
                        for name in self.__axis_names]):
                raise TypeError("Mismatch in the two aeseries axes.")
            return self.data <= other.data
        else:
            return self.data <= other
    
    def __lt__(self, other):
        if isinstance(other, aeseries):
            if not all([hasattr(other, name) for name in self.__axis_names]):
                raise TypeError("aeserieses do not have the same attributes")
            if not all([np.all(getattr(self, name) == getattr(other, name))
                        for name in self.__axis_names]):
                raise TypeError("Mismatch in the two aeseries axes.")
            return self.data < other.data
        else:
            return self.data < other
    
    def __gt__(self, other):
        if isinstance(other, aeseries):
            if not all([hasattr(other, name) for name in self.__axis_names]):
                raise TypeError("aeserieses do not have the same attributes")
            if not all([np.all(getattr(self, name) == getattr(other, name))
                        for name in self.__axis_names]):
                raise TypeError("Mismatch in the two aeseries axes.")
            return self.data > other.data
        else:
            return self.data > other
            
    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """
        Handle NumPy functions like np.sin, np.exp.
        """
        input_values = []
        for i in inputs:
            if isinstance(i, aeseries):
                input_values.append(i.data)
            else:
                input_values.append(i)
        result = getattr(ufunc, method)(*input_values, **kwargs)
        return aeseries(result, **{
                name: getattr(self, name) for name in self.__axis_names
            })
    
    def __array_function__(self, func, types, args, kwargs):
        if func == np.fft.fft:
            return self.fft(**kwargs)
        elif func == np.fft.rfft:
            return self.rfft(**kwargs)
        elif func in [scipy.signal.stft, scipy.signal.ShortTimeFFT]:
            return self.stft(**kwargs)
        else:
            return aeseries(func, **{
                name: getattr(self, name) for name in self.__axis_names
            })
    
    def __array__(self, dtype=None):
        """Allows NumPy to recognize aeseries as an array-like object."""
        return np.asarray(self.data, dtype=dtype)
    
    def return_axis_names(self):
        return self.__axis_names

    def to(self, unit):
        """
        converts "data" to the selected unit
        """
        return aeseries(self.data.to(unit), **{
                name: getattr(self, name) for name in self.__axis_names
            })
    
    def fft(self, norm:Literal['backward', 'ortho', 'forward']='backward',
            norm_by_max=False, time_range=None, axis=-1,
            check_spacing=False, windowing:Literal['bartlett', 'blackman',
                                                   'hamming', 'hanning',
                                                   'kaiser']=None):
        """
        Returns an aeseries with the fourier transform of the signal.
        Time is necessary to perform this.add()
        **kwargs:
        norm{“backward”, “ortho”, “forward”}:  Indicates which direction
                    of the forward/backward pair of transforms is scaled
                    and with what normalization factor
        norm_by_max: bool if True the FFT is normalized by the maximum
                     value
        time_range: list of float or aerrays, crop the signal between
                    the two times
        axis: Axis over which to compute the FFT
        check_spacing: if true check if the time array is equally spaced,
                       can cause issue for small timesteps
        windowing{'bartlett', 'blackman', 'hamming', 'hanning', 'kaiser'}:
                Window to apply to the signal, default is no window
        returns the absolute value of the fft
        """
        if 'time' not in self.__axis_names:
            raise AttributeError("aeseries does not have the time attribute.")
        if self.time.shape != self.data.shape:
            raise TypeError("Mismacth in data and time shapes.")
        if check_spacing:
            if not np.all(np.diff(self.time) == self.time[1] - self.time[0]):
                raise ValueError("dt is not constant")
        indices = slice(0, len(self.time), 1)
        if time_range:
            if not len(time_range) == 2:
                warnings.warn("Ignoring time range, non matching input values")
            else:
                time_range.sort()
                indices = slice(np.argmax(self.time >= time_range[0]),
                                np.argmax(self.time >= time_range[1]), 1)
        if windowing is None:
            window = 1.
        else:
            window = getattr(np, windowing)(indices[1] - indices[0])
        freq = np.fft.fftfreq(len(self.time[indices]),
                              np.mean(np.diff(self.time[indices].to(u.s).value)))
        tilde_signal = np.abs(np.fft.fft(self.data[indices].value * window,
                                         norm=norm, axis=axis))
        if norm_by_max:
            tilde_signal /= np.abs(tilde_signal).max()
            units = u.dimensionless_unscaled
        else:
            units = self.data.unit * u.s
        ffreq = aerray(freq, u.Hz, 'frequency', r'$f$')
        ttilde_signal = aerray(tilde_signal, units,
                               merge_strings('fft_', self.data.name),
                               apply_symbol(self.data.label, r'\tilde'),
                               self.data.cmap,
                               [tilde_signal.min(), tilde_signal.max()],
                               self.data.log)
        return aeseries(ttilde_signal, frequency=ffreq)
    
    def rfft(self, norm:Literal['backward', 'ortho', 'forward']='backward',
             norm_by_max=False, time_range=None, axis=-1,
             check_spacing=False, windowing:Literal['bartlett', 'blackman',
                                                   'hamming', 'hanning',
                                                   'kaiser']=None):
        """
        Returns an aeseries with the fourier transform of the signal.
        Time is necessary to perform this.add()
        **kwargs:
        norm{“backward”, “ortho”, “forward”}:  Indicates which direction
                    of the forward/backward pair of transforms is scaled
                    and with what normalization factor
        norm_by_max: bool if True the FFT is normalized by the maximum
                     value
        time_range: list of float or aerrays, crop the signal between
                    the two times
        axis: Axis over which to compute the FFT
        check_spacing: if true check if the time array is equally spaced,
                       can cause issue for small timesteps
        windowing{'bartlett', 'blackman', 'hamming', 'hanning', 'kaiser'}:
                Window to apply to the signal, default is no window
        returns the absolute value of the fft
        """
        if 'time' not in self.__axis_names:
            raise AttributeError("aeseries does not have the time attribute.")
        if self.time.shape != self.data.shape:
            raise TypeError("Mismacth in data and time shapes.")
        if check_spacing:
            if not np.all(np.diff(self.time) == self.time[1] - self.time[0]):
                raise ValueError("dt is not constant")
        indices = slice(0, len(self.time), 1)
        if time_range:
            if not len(time_range) == 2:
                warnings.warn("Ignoring time range, non matching input values")
            else:
                time_range.sort()
                indices = slice(np.argmax(self.time >= time_range[0]),
                                np.argmax(self.time >= time_range[1]), 1)
        if windowing is None:
            window = 1.
        else:
            window = getattr(np, windowing)(len(self.time[indices]))
        freq = np.fft.rfftfreq(len(self.time[indices]),
                              np.mean(np.diff(self.time[indices].to(u.s).value)))
        tilde_signal = np.abs(np.fft.rfft(self.data[indices].value * window,
                                          norm=norm, axis=axis))
        if norm_by_max:
            tilde_signal /= np.abs(tilde_signal).max()
            units = u.dimensionless_unscaled
        else:
            units = self.data.unit * u.s
        ffreq = aerray(freq, u.Hz, 'frequency', r'$f$', None, [0, 2000])
        ttilde_signal = aerray(tilde_signal, units,
                               merge_strings('fft_', self.data.name),
                               apply_symbol(self.data.label, r'\tilde'),
                               self.data.cmap,
                               [tilde_signal.min(), tilde_signal.max()],
                               True)
        return aeseries(ttilde_signal, frequency=ffreq)

    def stft(self, window_size=aerray(10, u.ms), check_spacing=False,
             time_range=None, scale_to:Literal['magnitude', 'psd']='magnitude',
             windowing:Literal['bartlett', 'blackman', 'hamming', 'hanning',
                                'kaiser']='hanning', overlap=0.5,
             highpass=None, lowpass=None, bandpass=None):
        """
        Returns an aeseries with the short time fourier transform of the signal.
        Time is necessary to perform this.
        **kwargs:
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
        filters are avaiable: lowpass, highpass and bandpass
        returns the absolute value of the fft
        """
        def butter_filter(signal, fs, cutoff, btype, order=4):
            nyquist = 0.5 * fs
            normal_cutoff = np.array(cutoff) / nyquist
            b, a = scipy.signal.butter(order, normal_cutoff, btype=btype,
                                       analog=False)
            signal[:] = scipy.signal.filtfilt(b, a, signal.value)
            return signal
        
        if 'time' not in self.__axis_names:
            raise AttributeError("aeseries does not have the time attribute.")
        if self.time.shape != self.data.shape:
            raise TypeError("Mismacth in data and time shapes.")
        if check_spacing:
            if not np.all(np.diff(self.time) == self.time[1] - self.time[0]):
                raise ValueError("dt is not constant")
        istart = 0
        istop = len(self.time)
        if time_range:
            if not len(time_range) == 2:
                warnings.warn("Ignoring time range, non matching input values")
            else:
                time_range.sort()
                istart = np.argmax(self.time >= time_range[0])
                istop = np.argmax(self.time >= time_range[1])
        indices = slice(istart, istop, 1)
        win_len = 0
        i = istart+1
        while self.time[i] - self.time[istart] < window_size:
            win_len += 1
            i += 1
        hop = int(overlap * win_len)
        window = getattr(np, windowing)(win_len)
        fs = 1 / np.mean(np.diff(self.time[indices].to(u.s).value))
        signal = self.data[indices].copy()
        if highpass is not None:
            signal = butter_filter(signal, cutoff=highpass, fs=fs, btype='high')
        if lowpass is not None:
            signal = butter_filter(signal, cutoff=lowpass, fs=fs, btype='low')
        if isinstance(bandpass, list):
            if len(bandpass) == 2:
                signal = butter_filter(signal, cutoff=bandpass, fs=fs,
                                       btype='band')
        SFT = scipy.signal.ShortTimeFFT(window, hop, fs, scale_to=scale_to)
        freq = SFT.f
        Zxx = np.abs(SFT.stft(signal.value))
        tm = SFT.t(len(self.time)) + self.time[0].value
        ffreq = aerray(freq, u.Hz, 'frequency', r'$f$', None, [0, 2000])
        ttm = aerray(tm, u.s, self.time.name, self.time.label,
                     self.time.cmap, self.time.limits, self.time.log)
        if scale_to == 'magnitude':
            uu = self.data.unit
            lab = merge_strings(r'Amplitude $|$', self.data.label, r'$|(f)$')
        else:
            uu = self.data.unit ** 2 / u.Hz
            lab = merge_strings(r'PSD$_{$', self.data.label, r'$}(f)$')
        ZZxx = aerray(Zxx, uu, self.data.name+f'_{scale_to}', lab,
                      self.data.cmap, [np.min(Zxx[np.nonzero(Zxx)]) * 1.1, Zxx.max() * 0.9],
                      True)
        return aeseries(ZZxx, time=ttm, frequency=ffreq)
    
    def parameters(self):
        """
        Returns the parameters of the aeseries
        """
        return self.__axis_names
    
    def derive(self, axis):
        """
        Derive the data along the selected axis
        """
        if axis not in self.__axis_names:
            raise AttributeError("Axis not found in the aeseries")
        return aeseries(IDL_derivative(getattr(self, axis), self.data,
                                           axis=self.__axis_indices[axis]),
                        **{name: getattr(self, name) for name in self.__axis_names})
        
    def integrate(self, axis):
        """
        linearly integrate along the selected axis
        """
        if axis not in self.__axis_names:
            raise AttributeError("Axis not found in the aeseries")
        xvar = getattr(self, axis)
        dvar_label = r'$\mathrm{d}' + \
            xvar.label.replace('$', '').split('-')[0] + r'$'
        integral = cumulative_simpson(self.data.value, x=xvar.value,
                                      axis=self.__axis_indices[axis])
        integral = aerray(
            np.insert(integral, 0, integral[0], self.__axis_indices[axis]),
            unit = (self.data.unit * xvar.unit),
            name = self.data.name + '_' + xvar.name + '_integrated',
            label = merge_strings(r'$\int$', dvar_label, self.data.label),
            cmap = self.data.cmap,
            limits = [integral.min() * 1.1, integral.max() * 0.9],
            log = self.data.log
        )
        return aeseries(integral,
                        **{name: getattr(self, name) for name in self.__axis_names})