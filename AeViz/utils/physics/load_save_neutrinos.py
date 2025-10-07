from typing import Literal
import numpy as np
from AeViz.utils.utils import progressBar, check_existence, checkpoints
from AeViz.units.aeseries import aerray, aeseries
from AeViz.units import u
import h5py, os
from AeViz.utils.files.file_utils import save_hdf, create_series

FILE_H5 = 'neutrino_luminosity.h5'

def calculate_luminosity(simulation, save_checkpoints=True, rmax=5e7, \
                            **kwargs):
    if check_existence(simulation, FILE_H5):
        time, L, processed_hdf = read_luminosity(simulation)
        ## Retrocompatibility option
        if processed_hdf is None:
            if len(simulation.hdf_file_list) == len(time):
                save_hdf(os.path.join(simulation.storage_path, FILE_H5),
                        ['time', 'luminosity', 'processed'],
                        [time, L, simulation.hdf_file_list])
                return create_series(time, L['nue'], L['nua'], L['nux'])
            else:
                start_point = 0
                time = 0
                L = {'nue' : 0.0, 'nua' : 0.0, 'nux' : 0.0}
                processed_hdf = []
        elif processed_hdf[-1].decode("utf-8") == simulation.hdf_file_list[-1]:
            return create_series(time, L['nue'], L['nua'], L['nux'])
        else:
            start_point = len(processed_hdf)
            processed_hdf = [ff.decode("utf-8") for ff in processed_hdf]
            print('Checkpoint found for neutrino luminosity, starting' \
                ' from checkpoint.\nPlease wait...')
    else:
        start_point = 0
        processed_hdf = []
        print('No checkpoint found for neutrino luminosity, starting from' \
            ' the beginning.\nPlease wait...')
    if (checkpoints[simulation.dim] == False) or (not save_checkpoints):
        checkpoint = len(simulation.hdf_file_list)
    else:
        checkpoint = checkpoints[simulation.dim]

    idx = np.argmin(simulation.cell.radius(simulation.ghost) < (rmax * u.cm))
    check_index = 0
    progress_index = 0
    total_points = len(simulation.hdf_file_list) - start_point
    for file in simulation.hdf_file_list[start_point:]:
      progressBar(progress_index, total_points, suffix='Computing...')

      Lum = simulation.neutrino_luminosity(file, comp='all', **kwargs)

      Lnuetot = np.nanmean(Lum[0][..., idx:idx+1]) * 4.0 * np.pi * rmax * rmax
      Lnuatot = np.nanmean(Lum[1][..., idx:idx+1]) * 4.0 * np.pi * rmax * rmax
      Lnuxtot = np.nanmean(Lum[2][..., idx:idx+1]) * 4.0 * np.pi * rmax * rmax

      try:
        time = np.concatenate((time, simulation.time(file)))
        L = {key: np.concatenate((L[key], Ltot), axis=-1)
                  for (key, Ltot) in zip(['nue', 'nua', 'nux'],
                                         [Lnuetot, Lnuatot, Lnuxtot])}

      except:
        time = simulation.time(file)
        L = {key: Ltot for (key, Ltot) in zip(['nue', 'nua', 'nux'],
                                         [Lnuetot, Lnuatot, Lnuxtot])}
      processed_hdf.append(file)
      if (check_index >= checkpoint and save_checkpoints):
          print('Checkpoint reached, saving...\n')
          save_hdf(os.path.join(simulation.storage_path, FILE_H5),
                    ['time', 'luminosity', 'processed'],
                    [time, L, processed_hdf])
          check_index = 0
      check_index += 1
      progress_index += 1
    
    print('Computation completed, saving...')
    save_hdf(os.path.join(simulation.storage_path, FILE_H5),
                ['time', 'luminosity', 'processed'],
                [time, L, processed_hdf])
    print('Done!')

    time, L, _ = read_luminosity(simulation)
    return create_series(time, L['nue'], L['nua'], L['nux'])
    
def read_luminosity(simulation):
    data_h5 = h5py.File(os.path.join(simulation.storage_path, \
                                  FILE_H5), 'r')
    data = [
            aerray(data_h5['time'][...], u.s, 'time',
                    r'$t-t_\mathrm{b}$', None,
                    [-0.005, data_h5['time'][-1]]),
            {
            'nue' : aerray(data_h5['luminosity/nue'][...], u.erg / u.s, 'Lnue', \
                               r'$L_{\nu_e}$', None, [1e48, 1e53]),
            'nua' : aerray(data_h5['luminosity/nua'][...], u.erg / u.s, 'Lnua', \
                               r'$L_{\nu_a}$', None, [1e48, 1e53]),
            'nux' : aerray(data_h5['luminosity/nux'][...], u.erg / u.s, 'Lnux', \
                               r'$L_{\nu_x}$', None, [1e48, 1e53])
            }]
    if 'processed' in data_h5:
        data.append(data_h5['processed'][...])
    else:
        data.append(None)
    data_h5.close()
    return data
