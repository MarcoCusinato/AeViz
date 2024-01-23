import os, h5py
import numpy as np
import pandas as pd

dict_flavour = {'0': 'nue', '1': 'nua', '2': 'nux'}

def load_file(path_folder, file_name):
    """
    Load data files, three attempts are made
    - Attempt 1: use of loadtxt method of NumPy
    - Attempt 2: use of loadtxt method of NumPy, reading only the first n columns,
                 with n specified on the 2nd row of the .txt header file
    - Attempt 3: use of Pandas to generate missing data
    - Attempt 4: read file line by line filling a matrix with dimension NxM, where N
                 is the number of rows in the file and M the maximum number of colums.
                 It assumes the last line as the one with the more colums.
    """
    path = os.path.join(path_folder, file_name)
    assert os.path.exists(path), "Selected file does not exists"
    try:
        data = np.loadtxt(path)
    except:
        try:
            col_number = int(np.loadtxt(path.replace('.dat', '.txt'), skiprows=1) + 2)
            data = np.loadtxt(path, usecols=tuple(range(col_number)))
        except:
            try:
                head = list(np.genfromtxt(path.replace('.dat', '.txt'), dtype=str,
                                        delimiter=',', skip_footer=1))
                data_str = pd.read_table(path, dtype=str, sep='\s+', names=head,
                                        usecols=range(col_number))
                data_str = data_str.fillna('0')
                data_str = data_str.to_numpy()
                index_list = []
                for i in range(data_str.shape[0]):
                    try:
                        data_str[i,:].astype('float')
                    except:
                        index_list.append(i)
                data_str = np.delete(data_str, index_list,0)
                data = data_str.astype('float')
            except:
                with open(path, 'r') as f:
                    lines = f.readlines()
                    data = np.zeros( ( len( lines ), len( lines[-1].split() ) ) )
                    for (i, line) in zip( range( len( lines ) ), lines ):
                        line_data = np.array( line.split() ).astype(float)
                        data[i, :len( line_data )] = line_data
    return data

def find_column_changing_line(path_folder, file_name):
    """
    Loads a data file and returns the number of the line at which
    there is a numbero of columns change.
    """
    path = os.path.join(path_folder, file_name)
    assert os.path.exists(path), "Selected file does not exists"
    number_of_colums = None
    line_number = 1
    with open(path, 'r') as f:
        for line in f:
            columns = len(line.strip().split())
            if number_of_colums is None:
                number_of_colums = columns
            if number_of_colums != columns:
                line_number += 1
                break
            line_number += 1
    if line_number < 3:
        line_number = None
    return line_number


def save_h5(save_path, data_radius, data_average,
            indices, ghost_cells, corrected_for_tob):
        """
        Save data in h5 format
        """
        file_out = h5py.File(save_path, 'w')
        if 'neutrino' in save_path:
                neutrino_save(file_out, data_radius, data_average, indices,
                              ghost_cells, corrected_for_tob)
        else:
                standard_save(file_out, data_radius, data_average, indices,
                              ghost_cells, corrected_for_tob)
        file_out.close()

def standard_save(file_out, data_radius, data_average, indices, ghost_cells,
                  corrected_for_tob):
        """
        Save radii data:
                Contains the following dataset
                - radii: NumPy array with dimensions (phi, theta, time),
                         theta and phi with ghost cells (default 1),
                         contains the raidius value for every cell
                - time: NumPy array with dimension (time), contains time valuse for each radius value
                - min: NumPy array with dimension (time), contains minimum radius value
                - max: NumPy array with dimension (time), contains maximum radius value
                - average: NumPy array with dimension (time), contains average radius value
                - indices: NumPy array with dimension (time), contains file indices
                - tob_correction: bool, if time value is tob corrected
        """
        file_out.create_dataset('radii', data = data_radius)
        file_out.create_dataset('time', data = data_average[:, 0])
        file_out.create_dataset('min', data = data_average[:, 1])
        file_out.create_dataset('max', data = data_average[:, 2])
        file_out.create_dataset('average', data = data_average[:, 3])
        file_out.create_dataset('indices', data = indices)
        file_out.create_dataset('tob_correction', data = corrected_for_tob)
        g_group = file_out.create_group('ghost_cells')
        for (key, value) in ghost_cells.items():
                g_group.create_dataset(key, data = value)

def neutrino_save(file_out, data_radius, data_average, indices, ghost_cells,
                  corrected_for_tob):
        """
        Save nutrino sphere radii data:
                Contains three groups, one for each neutrino flavour, with
                the following dataset
                - radii: NumPy array with dimensions (phi, theta, time),
                         theta and phi with ghost cells (default 1),
                         contains the raidius value for every cell
                - time: NumPy array with dimension (time), contains time valuse for each radius value
                - min: NumPy array with dimension (time), contains minimum radius value
                - max: NumPy array with dimension (time), contains maximum radius value
                - average: NumPy array with dimension (time), contains average radius value
                - indices: NumPy array with dimension (time), contains file indices
                - tob_correction: bool, if time value is tob corrected
        """
        file_out.create_dataset('time', data = data_average[:, 0])
        nu_group = file_out.create_group('radii')
        max_group = file_out.create_group('max')
        min_group = file_out.create_group('min')
        average_group = file_out.create_group('average')
        for flavour in range(3):
                key = dict_flavour[str(flavour)]
                nu_group.create_dataset(key, data = data_radius[..., flavour, :])
                max_group.create_dataset(key, data = data_average[:, 3 * flavour + 1])
                min_group.create_dataset(key, data = data_average[:, 3 * flavour + 2])
                average_group.create_dataset(key, data = data_average[:, 3 * flavour + 3])
        file_out.create_dataset('indices', data = indices)
        file_out.create_dataset('tob_correction', data = corrected_for_tob)
        g_group = file_out.create_group('ghost_cells')
        for (key, value) in ghost_cells.items():
                g_group.create_dataset(key, data = value)
