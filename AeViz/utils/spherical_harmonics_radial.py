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
            Ylm = SpH.Ylm_norm( m, l, theta, phi )
            try:
                out_array[harm_index, :] = np.sum( rho * Ylm * dOmega[..., None],
                                        axis=tuple(range(simulation.dim-1)) ) 
            except:
                out_array[harm_index, :] = np.sum( rho * Ylm[..., None] * dOmega[..., None],
                                        axis=tuple(range(simulation.dim)) ) 
            harm_index += 1
    return out_array
   
def calculate_rho_decomposition(simulation, save_checkpoints=True):
    if check_existence(simulation, 'rho_decomposition_SpH.h5'):
        time, decomposition, processed_hdf = read_rho_decomposition(simulation, 
                                                                    4)
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
                               4)
            
            check_index = 0
        check_index += 1
        progress_index += 1
        findex += 1
    print('Computation completed, saving...')
    save_decomposition(simulation, decomposition, time, processed_hdf, 4)
    return True


def save_decomposition(simulation, decomposition, time, processed_hdf, lmax):
    keys = ['time']
    quantity = [time]
    dec_index = 0
    for l in range(lmax + 1):
        for m in range(-l, l + 1):
            keys.append('rho_l' + str(l) + 'm' + str(m))
            quantity.append(decomposition[dec_index, ...])
    keys.append('processed')
    quantity.append(processed_hdf)
    save_hdf(os.path.join(simulation.storage_path, 'rho_decomposition_SpH.h5'),
                keys, quantity)
    
def read_rho_decomposition(simulation, lmax):
    decomposition_data = h5py.File(os.path.join(simulation.storage_path, 
                                            'rho_decomposition_SpH.h5'), 'r')
    data = [
        decomposition_data['time'][...]
    ]
    dec_data = np.zeros((np.math.factorial(lmax) + 1,
                         len(simulation.cell.radius(simulation.ghost)),
                         len(decomposition_data['time'][...])))
    dec_index = 0
    for l in range(lmax + 1):
        for m in range(-l, l + 1):
            key = 'rho_l' + str(l) + 'm' + str(m)
            dec_data[dec_index, ...] = decomposition_data[key][...]
            dec_index += 1
    data.append(dec_data)
    data.append(decomposition_data['processed'][...])
    decomposition_data.close()
    return data
    
def get_sph_profile(simulation, l, m):
    decomposition_data = h5py.File(os.path.join(simulation.storage_path, 
                                            'rho_decomposition_SpH.h5'), 'r')
    key = 'rho_l' + str(l) + 'm' + str(m)
    data = decomposition_data[key][...]
    decomposition_data.close()
    return data
