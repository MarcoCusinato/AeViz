from AeViz.utils.physics.EMD_utils import (get_IMFs, HHT_spectra,
                                    instant_ampl,
                                    instant_freq,
                                    change_label,
                                    change_name)
from AeViz.simulation.methods import *
from AeViz.utils.files.file_utils import create_series
from AeViz.units import aerray, u
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
@smooth
@derive
@sum_tob
def IMFs(self, comp:Literal['h+eq', 'hxeq', 'h+pol', 'hxpol']='h+eq',
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
    kwargs.setdefault('res', True)
    if self.dim == 1:
        warnings.warn("No GWs in spherical symmetry.")
        return None, None, None
    elif self.dim == 2:
        comp = 'h+eq'
    if mode == 'EEMD':
        time, IMFs, residue = get_IMFs(self.storage_path, comp)
    else:
        GWs = self.GW_Amplitudes(comp=comp)
        if end_time is not None:
            index_end = np.argmax(GWs.time >= end_time)
        else:
            index_end = -20
        if start_time is not None:
            index_start = np.argmax(GWs.time >= start_time)
        else:
            index_start = 0
        GWs = GWs[index_start:index_end]
        emd = EMD()
        emd.emd(S=GWs.data.value, T=GWs.time.value, max_imf=max_imfs)
        time = GWs.time.value
        IMFs, residue = emd.get_imfs_and_residue()
    if len(IMFs) < max_imfs:
        max_imfs = len(IMFs)
    IMFs = IMFs[min_imfs:max_imfs, :]
   
    time = aerray(time, u.s, 'time', r'$t-t_\mathrm{b}$', None, 
                  [-0.005, time[-1]])
    out_IMFs = []
    i = 1
    for IMF in IMFs:
        out_IMFs.append(
            aerray(IMF, u.cm, f'IMF_{i}', r'IMF$_{' + str(i) +r'}$', None, [-150, 150],
                   False)
        )
        i += 1
    residue = aerray(residue, u.cm, 'residue', r'Res', None, [-150, 150], False)
    IMFs = create_series(time, out_IMFs)[0]
    if kwargs['res']:
        IMFs.append(create_series(time, residue)[0])
    return IMFs

def instantaneous_frequency(self, IMFs=None, strain:Literal['h+eq', 'hxeq',
                                                        'h+pol', 'hxpol']='h+eq',
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
    if IMFs is None:
        IMFs = self.IMFs(strain=strain, mode=mode,
                         min_imfs=min_imfs, max_imfs=max_imfs,
                         tob_corrected=tob_corrected)
    if IMFs is None:
        return None
    if not isinstance(IMFs, list):
        IMFs = [IMFs]
    ifreq = instant_freq(IMFs, **kwargs)
    freqs = []
    if len(IMFs) > 1:
        for i in range(0, len(IMFs)):
            freqs.append(
                aerray(ifreq[:, i], u.Hz, change_name(IMFs[i].data.name, 'IF'),
                       change_label(IMFs[i].data.label, 'IF'))
            )
    else:
        freqs = [aerray(ifreq[:, i], u.Hz, change_name(IMFs[0].data.name, 'IF'),
                       change_label(IMFs[0].data.label, 'IF'))]
    return create_series(IMFs[0].time, freqs)[0]

def instantaneous_amplitude(self, IMFs=None, strain:Literal['h+eq', 'hxeq',
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
    if IMFs is None:
        IMFs = self.IMFs(strain=strain, mode=mode,
                         min_imfs=min_imfs, max_imfs=max_imfs,
                         tob_corrected=tob_corrected)
    if IMFs is None:
        return None
    if not isinstance(IMFs, list):
        IMFs = [IMFs]
    iampl = instant_ampl(IMFs, **kwargs)
    ampl = []
    if len(IMFs) > 1:
        for i in range(0, len(IMFs)):
            ampl.append(
                aerray(iampl[:, i], u.Hz, change_name(IMFs[i].data.name, 'IA'),
                       change_label(IMFs[i].data.label, 'IA'))
            )
    else:
        ampl = [aerray(iampl[:, i], u.Hz, change_name(IMFs[0].data.name, 'IA'),
                       change_label(IMFs[0].data.label, 'IA'))]
    return create_series(IMFs[0].time, ampl)[0]

def HH_spectrum(self, IMFs=None, strain:Literal['h+eq', 'hxeq', 'h+pol',
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
    if IMFs is None:
        IMFs = self.IMFs(strain=strain, mode=mode, max_imfs=max_imfs,
                     min_imfs=min_imfs, tob_corrected=tob_corrected)[:-1]
    if IMFs is None:
        return None
    return HHT_spectra(IMFs, time_bins, freq_bins, **kwargs)
