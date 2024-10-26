from AeViz.spherical_harmonics.spherical_harmonics import SphericalHarmonics
import numpy as np
from AeViz.utils.utils import (check_existence, progressBar, checkpoints)
from AeViz.utils.file_utils import save_hdf
from AeViz.grid.grid import grid
import os, h5py

def Harmonics_decomposition_rho(simulation, file_name, theta, phi, dOmega, SpH,
                                lmax = 4):
    rho = simulation.rho(file_name)
    out_array = np.zeros((np.math.factorial(lmax) + 1, rho.shape[-1]))
    harm_index = 0
    for l in range( lmax + 1 ):
        for m in range( -l, l + 1 ):
            Ylm = SpH.Ylm_norm(m, l, theta, phi)
            out_array[harm_index, :] = np.sum( rho * Ylm[..., None] * dOmega[..., None],
                                        axis=tuple(range(simulation.dim-1)) )
            harm_index += 1
    return out_array

def Harmonics_decomposition_rho_msum(simulation, file_name, theta, phi, dOmega,
                                     SpH, lmax = 40):
    rho = simulation.rho(file_name)
    out_array = np.zeros((lmax+1, rho.shape[-1]))
    harm_index = 0
    for l in range( lmax + 1 ):
        for m in range( -l, l + 1 ):
            Ylm = SpH.Ylm_norm(m, l, theta, phi)
            out_array[harm_index, :] += np.sum( rho * Ylm[..., None] * dOmega[..., None],
                                        axis=tuple(range(simulation.dim-1)) )
        harm_index += 1
    return out_array
   
def calculate_rho_decomposition(simulation, save_checkpoints=True, msum=False):
    if msum:
        lmax = 40
        fname = 'rho_decomposition_SpH_msum.h5'
    else:
        lmax = 4
        fname = 'rho_decomposition_SpH.h5'
    if check_existence(simulation, fname):
        time, decomposition, processed_hdf = read_rho_decomposition(simulation, 
                                                                    lmax, msum)
        if processed_hdf[-1].decode("utf-8") == simulation.hdf_file_list[-1]:
            return True
        else:
            start_point = len(processed_hdf)
            processed_hdf = [ff.decode("utf-8") for ff in processed_hdf]
            print('Checkpoint found for the harmonics decomposition file, ' \
                  'starting from the beginning.\nPlease wait...')
    else:
        start_point = 0
        processed_hdf = []
        print('No checkpoint found for the harmonics decomposition file, ' \
              'starting from the beginning.\nPlease wait...')
    if (checkpoints[simulation.dim] == False) or (not save_checkpoints):
        checkpoint = len(simulation.hdf_file_list)
    else:
        checkpoint = checkpoints[simulation.dim]
    
    ## Set up the spherical harmonics
    SpH = SphericalHarmonics()
    ## Get the angular component
    dOmega = simulation.cell.dOmega(simulation.ghost)
    theta = simulation.cell.theta(simulation.ghost)
    phi = simulation.cell.phi(simulation.ghost)
    findex = start_point
    check_index = 0
    progress_index = 0
    total_points = len(simulation.hdf_file_list) - start_point
    for file in simulation.hdf_file_list[start_point:]:
        progressBar(progress_index, total_points,
                    suffix='Computing spherical harmonics...')
        if msum:
            in_data = Harmonics_decomposition_rho_msum(simulation, file, theta,
                                                       phi, dOmega, SpH)
        else:
            in_data = Harmonics_decomposition_rho(simulation, file, theta, phi,
                                                  dOmega, SpH)
        try:
            time = np.concatenate((time, simulation.time(file)))
            decomposition = np.concatenate((decomposition, in_data[..., None]),
                                           axis=-1)
        except:
            time = simulation.time(file)
            decomposition = in_data[..., None]
        processed_hdf.append(file)
        if (check_index >= checkpoint) and save_checkpoints:
            print('Checkpoint reached, saving...\n')
            save_decomposition(simulation, decomposition, time, processed_hdf,
                               lmax, msum)
            
            check_index = 0
        check_index += 1
        progress_index += 1
        findex += 1
    print('Computation completed, saving...')
    save_decomposition(simulation, decomposition, time, processed_hdf, lmax,
                       msum)
    return True


def save_decomposition(simulation, decomposition, time, processed_hdf, lmax, msum):
    keys = ['time']
    quantity = [time]
    if msum:
        file_name = 'rho_decomposition_SpH_msum.h5'
    else:
        file_name = 'rho_decomposition_SpH.h5'
    dec_index = 0
    for l in range(lmax + 1):
        if msum:
            keys.append('rho_l' + str(l))
            quantity.append(decomposition[dec_index, ...])
            dec_index += 1
        else:
            for m in range(-l, l + 1):
                keys.append('rho_l' + str(l) + 'm' + str(m))
                quantity.append(decomposition[dec_index, ...])
                dec_index += 1
    keys.append('processed')
    quantity.append(processed_hdf)
    save_hdf(os.path.join(simulation.storage_path, file_name),
                keys, quantity)
    
def read_rho_decomposition(simulation, lmax, msum):
    if msum:
        fname = 'rho_decomposition_SpH_msum.h5'
        data_dim = lmax + 1
    else:
        fname = 'rho_decomposition_SpH.h5'
        data_dim = np.math.factorial(lmax) + 1
    decomposition_data = h5py.File(os.path.join(simulation.storage_path, 
                                            fname), 'r')
    data = [
        decomposition_data['time'][...]
    ]

    dec_data = np.zeros((data_dim,
                         len(simulation.cell.radius(simulation.ghost)),
                         len(decomposition_data['time'][...])))
    dec_index = 0
    for l in range(lmax + 1):
        if msum:
            key = 'rho_l' + str(l)
            dec_data[dec_index, ...] = decomposition_data[key][...]
            dec_index += 1
        else:
            for m in range(-l, l + 1):
                key = 'rho_l' + str(l) + 'm' + str(m)
                dec_data[dec_index, ...] = decomposition_data[key][...]
                dec_index += 1
    data.append(dec_data)
    data.append(decomposition_data['processed'][...])
    decomposition_data.close()
    return data
    
def get_sph_profile(simulation, l, m=None):
    if m is None:
        fname = 'rho_decomposition_SpH_msum.h5'
        key = 'rho_l' + str(l)
    else:
        fname = 'rho_decomposition_SpH.h5'
        key = 'rho_l' + str(l) + 'm' + str(m)
    decomposition_data = h5py.File(os.path.join(simulation.storage_path, 
                                            fname), 'r')
    data = decomposition_data[key][...]
    time = decomposition_data['time'][...]
    decomposition_data.close()
    return time, data


def get_sph_profiles_r(simulation, l, m=None, zero_norm=True,
                       rhomin=None, rhomax=None, r=None):
    rr = [rhomin, rhomax, r]
    assert rr.count(None) < 3, "Please provide at least one of the three " \
        "arguments: rmin, rmax, r"
    if m is None:
        time, r00 = get_sph_profile(simulation, 0)
    else:
        time, r00 = get_sph_profile(simulation, 0, 0)
    _, rlm = get_sph_profile(simulation, l, m)
    if zero_norm:
        rlm /= r00
    if r is not None:
        radius = simulation.cell.radius(simulation.ghost)
        rindex = np.argmax(radius >= r)
        return time, rlm[rindex, ...]
    else:
        if rhomin is None:
            rhomin = 0
        if rhomax is None:
            rhomax = r00.max()
        mask = (r00 >= rhomin) & (r00 <= rhomax)
        rlm[~mask] = 0
        return time, rlm[mask].sum(axis=0)
