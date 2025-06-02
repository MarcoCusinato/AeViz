import os, h5py
import numpy as np
import pandas as pd
import inspect
from AeViz.units import aeseries, aerray, u
import requests

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

def find_column_changing_line(path_folder, file_name, column=None):
    """
    Loads a data file and returns the number of the line at which
    there is a numbero of columns change.
    """
    default_column = 0
    if column is None:
        column = default_column
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
    if len(line_change) > 1:
        line_number = line_change[column]
    if line_number < 3:
        line_number = None
    return line_number

def save_hdf(save_path, dataset_keywords, dataset_values):
    """
    Save data in hdf format
    dataset_keywords: list of strings, keywords for the datasets
    dataset_values: list of whatever you want, these are the values for
                    the datasets
    """
    assert len(dataset_keywords) == len(dataset_values), \
        "Number of keywords and values do not match"
    file_out = h5py.File(save_path, 'w')
    for (key, value) in zip(dataset_keywords, dataset_values):
        if type(value) == dict:
            group = file_out.create_group(key)
            for (k, v) in value.items():
                if type(v) == dict:
                    subgroup = group.create_group(k)
                    for (kk, vv) in v.items():
                        if isinstance(vv, aerray):
                            subgroup.create_dataset(kk, data = vv.value)
                        else:
                            subgroup.create_dataset(kk, data = vv)
                else:
                    if isinstance(v, aerray):
                        group.create_dataset(k, data = v.value)
                    else:
                        group.create_dataset(k, data = v)
        else:
            file_out.create_dataset(key, data = value)
    file_out.close()
    
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
        'ET': "https://apps.et-gw.eu/tds/?call_file=ET-0000A-18_ETDSensitivityCurveTxtFile.txt"
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