from __future__ import annotations
import os, h5py
import numpy as np
import pandas as pd
import inspect
from AeViz.units import aeseries, aerray, u
import requests
from AeViz.utils.utils import units_from_string, check_existence
from AeViz.simulation.simulation import Simulation

def list_module_functions(module):
    """
    Lists all the functions decleared in a specific module
    Arguments:
        - imported module
    Returns:
        list of all the functions decleared in the module
    """
    return [
        (name, obj) for name, obj in inspect.getmembers(module, inspect.isfunction)
        if obj.__module__ == module.__name__ 
    ]

def load_file(path_folder, file_name):
    """
    Load data files, three attempts are made
    - Attempt 1: use of loadtxt method of NumPy
    - Attempt 2: use of loadtxt method of NumPy, reading only the first 
                 n columns, with n specified on the 2nd row of the .txt
                 header file
    - Attempt 3: use of Pandas to generate missing data
    - Attempt 4: read file line by line filling a matrix with dimension 
                 NxM, where N is the number of rows in the file and M 
                 the maximum number of colums. It assumes the last line 
                 as the one with the more colums.
    - Attempt 5: read file line by line filling a matrix with dimension
                 NxM, where N is the number of rows in the file and M
                 the maximum number of colums. This time we process all
                 the lines to get the maximum number of colums.
                 Then we fill up the matrix with the data.
    - Attempt 6: pray
    """
    path = os.path.join(path_folder, file_name)
    assert os.path.exists(path), "Selected file does not exists"
    try:
        data = np.loadtxt(path)
    except:
        try:
            col_number = int(np.loadtxt(path.replace('.dat', '.txt'), 
                                        skiprows=1) + 2)


            data = np.loadtxt(path, usecols=tuple(range(col_number)))
        except:
            try:
                head = list(np.genfromtxt(path.replace('.dat', '.txt'),
                                          dtype=str, delimiter=',',
                                          skip_footer=1))
                data_str = pd.read_table(path, dtype=str, sep=r'\s+',
                                         names=head, usecols=range(col_number))
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
                try:
                    with open(path, 'r') as f:
                        lines = f.readlines()
                        data = np.zeros((len(lines), len(lines[-1].split())))
                        for (i, line) in zip(range(len(lines)), lines):
                            line_data = np.array( line.split() ).astype(float)
                            data[i, :len( line_data )] = line_data
                except:
                    with open(path, 'r') as f:
                        n_cols = -1
                        for line in f:
                            n_cols = max(n_cols, len(line.split()))
                        data = np.zeros((len(lines), n_cols))
                        for (i, line) in zip(range(len(lines)), lines):
                            line_data = np.array( line.split() ).astype(float)
                            data[i, :len( line_data )] = line_data
                        
    return data

def find_column_changing_line(path_folder: str,
                              file_name: str) -> list:
    """
    Loads a data file and returns the list of line numbers at which
    number of the line at which the number of columns changes.

    Parameters
    ----------
    path_folder : str
        path of the folder in which the file is located
    file_name : str
        name of the file to load

    Returns
    -------
    list
        list of row numbers
    """
    #default_column = 0
    #if column is None:
    #    column = default_column
    path = os.path.join(path_folder, file_name)
    assert os.path.exists(path), "Selected file does not exists"
    number_of_colums = None
    line_number = 1
    line_change = []
    with open(path, 'r') as f:
        for line in f:
            columns = len(line.strip().split())
            if number_of_colums is None:
                number_of_colums = columns
            if number_of_colums != columns:
                line_number += 1
                number_of_colums = columns
                line_change.append(line_number)
                #break
            line_number += 1
    #if len(line_change) > 1:
    #    line_number = line_change[column]
    #if len(line_change) == 1:
    #    line_number = line_change[0]
    #if line_number < 3:
    #    line_number = None
    return line_change

def save_hdf(save_path, dataset_keywords, dataset_values):
    """
    Save data in hdf format
    dataset_keywords: list of strings, keywords for the datasets
    dataset_values: list of whatever you want, these are the values for
                    the datasets
    """
    def save_with_attr(arr, key, group):
        d = group.create_dataset(key if isinstance(key, str) else str(key),
                                 data = arr.value)
        d.attrs['name'] = arr.name
        d.attrs['label'] = arr.label
        d.attrs['cmap'] = arr.cmap
        d.attrs['lim0'] = arr.limits[0]
        d.attrs['lim1'] = arr.limits[1]
        d.attrs['log'] = arr.log
        d.attrs['unit'] = str(arr.unit)
    assert len(dataset_keywords) == len(dataset_values), \
        "Number of keywords and values do not match"
    with h5py.File(save_path, 'w') as file_out:
        for (key, value) in zip(dataset_keywords, dataset_values):
            if type(value) == dict:
                group = file_out.create_group(key)
                for (k, v) in value.items():
                    if type(v) == dict:
                        subgroup = group.create_group(k)
                        for (kk, vv) in v.items():
                            if isinstance(vv, aerray):
                                save_with_attr(vv, kk, subgroup)
                            else:
                                subgroup.create_dataset(
                                    kk if isinstance(kk, str) else str(kk),
                                    data = vv)
                    else:
                        if isinstance(v, aerray):
                            save_with_attr(v, k, group)
                        else:
                            group.create_dataset(
                                k if isinstance(k, str) else str(k),
                                data = v)
            elif isinstance(value, aerray):
                save_with_attr(value, key, file_out)
            else:
                file_out.create_dataset(key, data = value)

def save_merge_dictionary_hdf(simulation: Simulation,
                              value_dictionary: dict,
                              save_path: str,
                              file_name: str,
                              **kwargs
                              ) -> None:
    """
    Merges the postprocessing located in a file with the newly run one.
    In case no postprocessing is found, the file is created.
    The output file needs to be organised as follows: two lists
    containing time and processed_hdf files, a dictionary containing
    local and global dictionaries of lists.
    

    Parameters
    ----------
    simulation : Simulation
        simulation object, needed to save in the correct path
    value_dictionary : dict
        time, dictionary of quantities, file processed list
    save_path : str
        path to the save folder
    file_name : str
        name of the file to load and save
    """
    time, out_dictionary, processed_hdf = value_dictionary
    time = np.concatenate(time)
    out_dictionary['global'] = {kk: np.concatenate(vv) for (kk, vv) 
                                in out_dictionary['global'].items()}
    out_dictionary['local'] = {kk: np.stack(vv, axis=-1) for (kk, vv) 
                                    in out_dictionary['local'].items()}
    if check_existence(simulation, os.path.join(save_path, file_name)):
        old_t, old_dict, old_proc = load_hdf_to_dictionary(save_path,
                                                           file_name)
        processed_hdf = old_proc.extend(processed_hdf)
        nm, lb, lg, cm = old_t.name, old_t.label, old_t.log, old_t.cmap
        time = np.concatenate(old_t, time)
        time.set(limits = [-0.005, time.value[-1]])
        time.set(name=nm, label=lb, limits=lm, cmap=cm, log=lg)
        for key in out_dictionary['global'].keys():
            nm = old_dict['global'][key].name
            lb = old_dict['global'][key].label
            lm = old_dict['global'][key].limits
            lg = old_dict['global'][key].log
            cm = old_dict['global'][key].cmap
            out_dictionary['global'][key] = \
                np.concatenate(old_dict['global'][key],
                               out_dictionary['global'][key])
        for key in out_dictionary['local'].keys():
            nm = old_dict['local'][key].name
            lb = old_dict['local'][key].label
            lm = old_dict['local'][key].limits
            lg = old_dict['local'][key].log
            cm = old_dict['local'][key].cmap
            out_dictionary['local'][key] = \
                np.concatenate(old_dict['local'][key],
                                out_dictionary['local'][key])
    else:
        time.set(name='time', label=r'$t-t_{\rm b}$',
                 limits=[-0.005, time.value[-1]], log=False)
        for key in out_dictionary['global'].keys():
            out_dictionary['global'][key].set(**kwargs['global'][key])
        for key in out_dictionary['local'].keys():
            out_dictionary['local'][key].set(**kwargs['local'][key])
            
    save_hdf(os.path.join(save_path, file_name),
             ['time', 'local', 'global', 'processed_hdf'],
             [time, out_dictionary['local'], out_dictionary['global'],
              processed_hdf])

def load_hdf_to_dictionary(path: str,
                           file_name:str
                           ) -> tuple[aerray, dict[aerray], list[str]]:
    """
    Loads a hdf file split into 'global' and 'local' datasets into a
    dictionary of aerrays.

    Parameters
    ----------
    path : str
        path to the file to load
    file_name : str
        name of the file to load

    Returns
    -------
    tuple[aerray, dict[aerray], list[str]]
        aerray containing the time and a dictionary containing two
        dictionaries of aerrays.
    """
    with h5py.File(os.path.join(path, file_name), 'r') as f:
        time = aerray(f['time'][:], units_from_string(f['time'].attrs['unit']),
                        f['time'].attrs['name'], f['time'].attrs['label'],
                        f['time'].attrs['cmap'], [f['time'].attrs['lim0'],
                                                f['time'].attrs['lim1']],
                        f['time'].attrs['log'])
        out_dictionary = {'global': {},
                            'local': {}}
        for key in f['global'].keys():
            att = f[f'global/{key}'].attrs
            out_dictionary['global'][key] = aerray(f[f'global/{key}'][:],
                                                    units_from_string(att['unit']),
                                                    att['name'],
                                                    att['label'],
                                                    att['cmap'],
                                                    [att['lim0'], att['lim1']],
                                                    att['log'])
        for key in f['local'].keys():
            att = f[f'local/{key}'].attrs
            out_dictionary['local'][key] = aerray(f[f'local/{key}'][...],
                                                    units_from_string(att['unit']),
                                                    att['name'],
                                                    att['label'],
                                                    att['cmap'],
                                                    [att['lim0'], att['lim1']],
                                                    att['log'])
        processed_hdf = [ff.decode("utf-8") for ff in f['processed_hdf']]
    return time, out_dictionary, processed_hdf

def load_dataset(path: str,
              file_name: str,
              dset_name: str,
              return_series: bool = True) -> aerray | aeseries:
    with h5py.File(os.path.join(path, file_name), 'r') as f:
        dset = f[dset_name]
        if len(dset.attrs) > 0:
            att = dset.attrs
            arr = aerray(dset[...],
                         units_from_string(att['unit']),
                         att['name'],
                         att['label'],
                         att['cmap'],
                         [att['lim0'], att['lim1']],
                         att['log'])
        else:
            arr = dset[...]
            if all((isinstance(a, bytes) for a in arr)):
                arr = [a.decode("utf-8") for a in arr]
        if isinstance(arr, aerray) and return_series:
            dset = f['time']
            att = dset.attrs
            time = aerray(dset[...],
                            units_from_string(att['unit']),
                            att['name'],
                            att['label'],
                            att['cmap'],
                            [att['lim0'], att['lim1']],
                            att['log'])
            return create_series(time, arr)
        return arr

def create_series(time, *args):
    """
    Creates as many aeseries as argument.
    """
    ghost_cells = False
    if type(args[-1]) == dict:
        if 'r_l' in args[-1]:
            ghost_cells = args[-1]
            args = args[:-1]
    series = []
    for arg in args:
        if type(arg) == dict:
            ddict = {}
            for key in arg.keys():
                if isinstance(arg[key], aerray):
                    ddict[key] = aeseries(arg[key], time=time.copy())
                elif type(arg[key]) == dict:
                    dddict = {}
                    for kk in arg[key].keys():
                        dddict[kk] = aeseries(arg[key][kk], time=time.copy())
                    ddict[key] = dddict                   
            series.append(ddict)
        elif type(arg) == list:
            llist = []
            for a in arg:
                llist.append(aeseries(a, time=time.copy()))
            series.append(llist)
        else:
            series.append(aeseries(arg, time=time.copy()))
    if ghost_cells:
        series.append(ghost_cells)
    return series

def load_asd(path, detector):
    url = {
        'LIGOO3H': "https://dcc.ligo.org/public/0165/T2000012/002/aligo_O3actual_H1.txt",
        'LIGOO3L': "https://dcc.ligo.org/public/0165/T2000012/002/aligo_O3actual_L1.txt",
        'LIGOO4': "https://dcc.ligo.org/public/0165/T2000012/002/aligo_O4low.txt",
        'LIGOO4High': "https://dcc.ligo.org/public/0165/T2000012/002/aligo_O4high.txt",
        'LIGO': "https://dcc.ligo.org/public/0165/T2000012/002/AplusDesign.txt",
        'VirgoO3': "https://dcc.ligo.org/public/0165/T2000012/002/avirgo_O3actual.txt",
        'VirgoO4': "https://dcc.ligo.org/public/0165/T2000012/002/avirgo_O4high_NEW.txt",
        'Virgo': "https://dcc.ligo.org/public/0165/T2000012/002/avirgo_O5low_NEW.txt",
        'VirgoO5High': "https://dcc.ligo.org/public/0165/T2000012/002/avirgo_O5low_NEW.txt",
        'KAGRA': "https://dcc.ligo.org/public/0165/T2000012/002/kagra_128Mpc.txt",
        'KAGRA80': "https://dcc.ligo.org/public/0165/T2000012/002/kagra_80Mpc.txt",
        'ET': "https://apps.et-gw.eu/tds/?call_file=ET-0000A-18_ETDSensitivityCurveTxtFile.txt",
        'CE': ""
    }
    assert detector in url.keys(), f"The detector should be one of {list(url.keys())}"
    file_path = os.path.join(path, 'psds', detector + '.txt')
    if not os.path.exists(file_path):
        if not os.path.exists(os.path.join(path, 'psds')):
            os.mkdir(os.path.join(path, 'psds'))
        resp = requests.get(url[detector])
        with open(file_path, 'wb') as f:
            f.write(resp.content)
    psd = np.loadtxt(file_path)
    frequency = aerray(psd[:, 0], u.Hz, 'frequency', r'$f$', None, [10, 4e3],
                       True)
    i10 = np.argmax(frequency >= 10)
    i4000 = np.argmax(frequency >= 4000)
    if detector == 'ET':
        asd = aerray(psd[:, 3], (u.Hz**(-0.5)), detector, detector, None,
                     [psd[i10:i4000, 3].min(), psd[i10:i4000, 3].max()], True)
    else:
        asd = aerray(psd[:, 1], (u.Hz**(-0.5)), detector, detector, None,
                     [psd[i10:i4000, 1].min(), psd[i10:i4000, 1].max()], True)
    
    return aeseries(asd,
                    frequency=frequency)