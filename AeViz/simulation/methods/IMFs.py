from AeViz.utils.EMD_utils import (get_IMFs, HHT_spectra,
                                    instantaneous_amplitude,
                                    instantaneous_frequency)
import numpy as np
from PyEMD import EMD
from typing import Literal
import warnings

"""
Functions to process gravitational from a simulation in
spherical coordinates waves with the EEMD and EMD
decomposition.
These functions are not meant to be used standalone, but rather to be
imported into the Simulation class.
"""

## -----------------------------------------------------------------
## IMFs
## -----------------------------------------------------------------

def IMFs(self, strain:Literal['h+eq', 'hxeq', 'h+pol', 'hxpol']='h+eq',
            mode:Literal['EMD', 'EEMD']='EMD', min_imfs=0, max_imfs=10,
            start_time=-0.05, end_time=None, tob_corrected=True, **kwargs):
    """
    Returns the Intrinsic Mode Functions of the GWs strain.
    If mode is 'EEMD' the IMFs are calculated with the Ensemble
    Empirical Mode Decomposition, otherwise with the Hilbert-Huang
    Transform.
    If tob_corrected is False the IMFs are not corrected for the
    time of bounce.
    Returns: time, IMFs, residue
    """
    if self.dim == 1:
        warnings.warn("No GWs in spherical symmetry.")
        return None, None, None
    elif self.dim == 2:
        strain = 'h+eq'
    if mode == 'EEMD':
        time, IMFs, residue = get_IMFs(self.storage_path, strain)
        if not tob_corrected:
            time += self.tob
    else:
        GWs = self.GW_Amplitudes(tob_corrected)
        time = GWs[:, 0]
        if strain == 'h+eq':
            h = GWs[:, 1]
        elif strain == 'h+pol':
            h = GWs[:, 2]
        elif strain == 'hxeq':
            h = GWs[:, 3]
        elif strain == 'hxpol':
            h = GWs[:, 4]
        if end_time is not None:
            index_end = np.argmax(time >= end_time)
        else:
            index_end = -20
        if start_time is not None:
            index_start = np.argmax(time >= start_time)
        else:
            index_start = 0
        time = time[index_start:index_end]
        h = h[index_start:index_end]
        emd = EMD()
        emd.emd(S=h, T=time, max_imf=max_imfs)
        IMFs, residue = emd.get_imfs_and_residue()
    if len(IMFs) < max_imfs:
        max_imfs = len(IMFs)
    return time, IMFs[min_imfs:max_imfs, :], residue

def instantaneous_frequency(self, time=None, IMFs=None, 
                            strain:Literal['h+eq', 'hxeq', 'h+pol',
                                            'hxpol']='h+eq',
                            mode:Literal['EMD', 'EEMD']='EMD', min_imfs=0,
                            max_imfs=10, tob_corrected=True, **kwargs):
    """
    Returns the instantaneous frequency of the GWs strain.
    If mode is 'EEMD' the IMFs are calculated with the Ensemble
    Empirical Mode Decomposition, otherwise with the Hilbert-Huang
    Transform.
    If tob_corrected is False the IMFs are not corrected for the
    time of bounce.
    Returns: time, instantaneous frequency
    """
    if time is None or IMFs is None:
        time, IMFs, *_ = self.IMFs(strain=strain, mode=mode,
                                    min_imfs=min_imfs, max_imfs=max_imfs,
                                    tob_corrected=tob_corrected)
    if IMFs is None:
        return None
    if len(IMFs) < max_imfs:
        max_imfs = len(IMFs)
    return time, \
        instantaneous_frequency(IMFs, time, **kwargs)

def instantaneous_amplitude(self, strain:Literal['h+eq', 'hxeq',
                                                    'h+pol', 'hxpol']='h+eq',
                            mode:Literal['EMD', 'EEMD']='EMD', min_imfs=0,
                            max_imfs=10, tob_corrected=True, **kwargs):
    """
    Returns the instantaneous amplitude of the GWs strain.
    If mode is 'EEMD' the IMFs are calculated with the Ensemble
    Empirical Mode Decomposition, otherwise with the Hilbert-Huang
    Transform.
    If tob_corrected is False the IMFs are not corrected for the
    time of bounce.
    Returns: time, instantaneous amplitude
    """
    time, IMFs, *_ = self.IMFs(strain=strain, mode=mode,
                                min_imfs=min_imfs, max_imfs=max_imfs,
                                tob_corrected=tob_corrected)
    if IMFs is None:
        return None
    if len(IMFs) < max_imfs:
        max_imfs = len(IMFs)
    return time, \
        instantaneous_amplitude(IMFs, **kwargs)

def HH_spectrum(self, strain:Literal['h+eq', 'hxeq', 'h+pol',
                                        'hxpol']='h+eq',
                mode:Literal['EMD', 'EEMD']='EMD', min_imfs=0, max_imfs=10,
                time_bins=None, freq_bins=100, tob_corrected=True,
                **kwargs):
    """
    Returns the Hilbert-Huang spectrum of the GWs strain.
    If mode is 'EEMD' the IMFs are calculated with the Ensemble
    Empirical Mode Decomposition, otherwise with the Hilbert-Huang
    Transform.
    If tob_corrected is False the IMFs are not corrected for the
    time of bounce.
    Returns: spectrogram, frequencies, time
    """
    time, IMFs, *_ = self.IMFs(strain=strain, mode=mode, max_imfs=max_imfs,
                                min_imfs=min_imfs,
                                tob_corrected=tob_corrected)
    if IMFs is None:
        return None
    return HHT_spectra(IMFs, time, time_bins, freq_bins, **kwargs)
