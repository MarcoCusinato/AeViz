import numpy as np
import h5py, os
from collections import defaultdict
from numpy.fft import fft, ifft, fftfreq
from PyEMD import EMD
from scipy.signal import hilbert
from scipy.signal import savgol_filter
from sparse import COO
from AeViz.utils.decorators import EMD_smooth


def polish_signal(GWs, frequency_cut):
    """
    Remove high frequency noise from the signal.
    """
    dt = GWs[1, 0] - GWs[0, 0]
    frequency_spectrum = fftfreq(GWs.shape[0], dt)
    if frequency_cut is None:
        frequency_cut = np.max(frequency_spectrum)
        return GWs, frequency_cut
    
    
    mask = (np.abs(frequency_spectrum) > frequency_cut)
    for i in range(1, GWs.shape[1]):
        fourier_transform = fft(GWs[:, i])
        fourier_transform[mask] = 0
        GWs[:, i] = np.real(ifft(fourier_transform))
    return GWs, frequency_cut

def remove_residuals(GWs):
    """
    Remove residuals from the signal.
    """
    emd = EMD()
    for i in range(1, GWs.shape[1]):
        emd.emd(GWs[:, i], GWs[:, 0])
        _, residual = emd.get_imfs_and_residue()
        GWs[:, i] -= residual
    return GWs

def save_IMFs(res, path, time, IMFs, residue, args):
    """
    Save the IMFs and the residual.
    """
    file_name = os.path.join(path, 'IMFs_' + args.strain + '_' + str(res) + '.h5')
    with h5py.File(file_name, 'w') as f:
        f.create_dataset('time', data=time)
        f.create_dataset('f_cut', data=args.cut_freq)
        f.create_dataset('Noise_seed', data=args.nseed)
        f.create_dataset('sigma', data=args.sigma)
        f.create_dataset('Realizations', data=args.nres)
        f.create_dataset('IMFs', data=IMFs)
        f.create_dataset('Residue', data=residue)
        f.create_dataset('res_removed', data=args.remove_residual)

def compact_IMFs(path, strain, files):
    all_imfs = defaultdict(list)
    residue = []
    seeds = []
    resolutions = []
    sigma = []
    fcut = []
    res_removed = []
    folder_path = os.path.join(path, 'EEMD')
    IMF_name = 'IMFs_' + strain + '_'
    for split in range(len(files)):
        with h5py.File(os.path.join(folder_path, IMF_name + str(split) + '.h5'),
                       'r') as f:
            IMFs = f['IMFs'][...]
            time = f['time'][...]
            seeds.append(f['Noise_seed'][...])
            resolutions.append(f['Realizations'][...])
            sigma.append(f['sigma'][...])
            residue.append(f['Residue'][...])
            fcut.append(f['f_cut'][...])
            res_removed.append(f['res_removed'][...])
        for (i_imf, imf) in enumerate(IMFs):
            all_imfs[i_imf].append(imf)
    all_imfs = dict(all_imfs)
    residue = np.array(residue).mean(axis=0)
    IMFs = np.array([np.array(imfs).mean(axis=0) for imfs in all_imfs.values()])
    
    out_file = os.path.join(path, 'IMFs_' + strain + '.h5')
    if os.path.exists(out_file):
        choice = input('File already exists, do you want to overwrite it? (y/out_file_name): ')
        if choice == 'y':
            pass
        else:
            out_file = os.path.join(path, choice)
    with h5py.File(out_file, 'w') as f:
        f.create_dataset('time', data=time)
        f.create_dataset('IMFs', data=IMFs)
        f.create_dataset('Residue', data=residue)
        f.create_dataset('Noise_seed', data=seeds)
        f.create_dataset('Realizations', data=resolutions)
        f.create_dataset('sigma', data=sigma)
        f.create_dataset('f_cut', data=fcut)
        f.create_dataset('res_removed', data=res_removed)

## METHODS for the HHT spectrum
def HHT_hist_bins(time, scale='linear', mode='time', bins=None):
    """
    Compute the bins for the HHT histogram.
    """
    if mode == 'time':
        min_value = time.min()
        max_value = time.max()
    elif mode == 'frequency':
        min_value = 0
        max_value = fftfreq(time.size, time[1] - time[0]).max()
    else:
        raise ValueError('Mode must be either time or frequency.')
    if bins is None:
        bins = np.sqrt(time.size) + 1
    if mode == 'time' and bins > 1001:
        bins = 1001
    elif mode == 'frequency' and bins > 100:
        bins = 100
    bins = int(bins)
    if scale == 'linear':
        bin_edge = np.linspace(min_value, max_value, bins)
    elif scale == 'log':
        bin_edge =  np.logspace(np.log10(min_value),
                                np.log10(max_value), bins)
    else:
        raise ValueError('Scale must be either linear or log.')

    bin_centre = (bin_edge[1:] + bin_edge[:-1]) / 2
    return bin_edge, bin_centre

@EMD_smooth
def instantaneous_amplitude(IMF, **kwargs):
    """
    Compute the istantaneous amplitude of the signal.
    """
    if IMF.ndim == 2:
        analytic_signal = np.zeros_like(IMF)
        for i in range(IMF.shape[0]):
            analytic_signal[i, :] = np.abs(hilbert(IMF[i, :]))
    else:
        analytic_signal = np.abs(hilbert(IMF))
    return analytic_signal

@EMD_smooth
def instantaneous_frequency(IMF, time, **kwargs):
    """
    Compute the istantaneous frequency of the signal.
    """
    dt = time[1] - time[0]
    sample_frequency = 1 / dt
    if IMF.ndim == 2:
        phase = np.zeros_like(IMF)
        for i in range(IMF.shape[0]):
            analytic_signal = hilbert(IMF[i, :])
            iphase = np.unwrap(np.angle(analytic_signal), axis=0)
            phase[i, :] = savgol_filter(iphase, 3, 1, deriv=1, axis=0)
    else:
        analytic_signal = hilbert(IMF)
        phase = np.unwrap(np.angle(analytic_signal), axis=0)
        phase = savgol_filter(phase, 3, 1, deriv=1, axis=0)
    return phase / (2.0 * np.pi) * sample_frequency

def HHT_spectra_single_if(IF, IA, frequency_edges):
    mask = ((IF >= frequency_edges[0]) & (IF < frequency_edges[-1]))
    f_inds = np.digitize(IF, frequency_edges) - 1
    t_ind = np.arange(f_inds.shape[0])
    coordinates = np.c_[f_inds.flatten(), t_ind.flatten()].T
    IA = IA.flatten()[mask]
    coordinates = coordinates[:, mask]
    spectra = COO(coordinates, IA, shape=(f_inds.shape[0], IF.size))
    spectra = spectra.todense()
    return spectra

def HHT_spectra(IMF, time, tbins=None, fbins=None, **kwargs):
    
    IF = instantaneous_frequency(IMF, time, **kwargs)
    IA = instantaneous_amplitude(IMF, **kwargs)
    _, time_centre = HHT_hist_bins(time, mode='time',
                                            bins=tbins)
    freq_edges, freq_centre = HHT_hist_bins(time, mode='frequency',
                                            bins=fbins)
    if IMF.ndim == 2:
        for i in range(IMF.shape[0]):
            sp = HHT_spectra_single_if(IF[i, :], IA[i, :], freq_edges)
            if i == 0:
                spectra = sp
            else:
                spectra += sp
    else:
        spectra = HHT_spectra_single_if(IF, IA, freq_edges)
    tinx = np.ones(time_centre.size, dtype=bool)
    finx = np.ones(freq_centre.size, dtype=bool)
    return spectra[np.ix_(finx, tinx)], freq_centre[finx], \
        time_centre[tinx]

def get_IMFs(storage_path, strain):
    """
    Loads the IMFs and the residue from the storage path.
    """
    save_path = os.path.join(storage_path, 'IMFs_' + strain + '.h5')
    if not os.path.exists(save_path):
        raise FileNotFoundError('File not found. Please compute the IMFs first.')
    with h5py.File(save_path, 'r') as f:
        IMFs = f['IMFs'][...]
        time = f['time'][...]
        res = f['Residue'][...]
    return time, IMFs, res