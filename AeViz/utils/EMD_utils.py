import numpy as np
import h5py, os
from collections import defaultdict
from numpy.fft import fft, ifft, fftfreq
from PyEMD import EMD

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
    for i in range(GWs.shape[1]):
        fourier_transform = fft(GWs[:, i])
        fourier_transform[mask] = 0
        GWs[:, i] = np.real(ifft(fourier_transform))
    return GWs, frequency_cut

def remove_residuals(GWs):
    """
    Remove residuals from the signal.
    """
    emd = EMD()
    for i in range(GWs.shape[1]):
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