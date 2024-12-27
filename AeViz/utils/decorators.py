import h5py, os
import matplotlib.pyplot as plt
import numpy as np
from AeViz.utils.math_utils import IDL_derivative
from scipy.ndimage import convolve
from functools import wraps

def hdf_isopen(func):
    """
    Takes as input the Simulation object and either the file name, or
    file index or time. If the file is not open, it opens it.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        if type(args[1]) is int:
            file = args[0].hdf_file_list[args[1]]
        elif type(args[1]) is float:
            file = args[0].find_file_from_time(args[1])
        else:
            file = args[1]
        if file != args[0]._Simulation__opened_hdf_file:
            if args[0]._Simulation__data_h5 is not None:
                args[0]._Simulation__data_h5.close()
            if file not in args[0].hdf_file_list:
                raise ValueError("Selected file does not exist.")
            args[0]._Simulation__opened_hdf_file = file
            args[0]._Simulation__data_h5 = h5py.File(
                os.path.join(args[0]._Simulation__hdf_path, file), 'r')
        return func(*args, **kwargs)
    return wrapper

def fig_window_open(func):
    """
    Check if the window figure has been closed. If it has it destroys it.
    """
    @wraps(func)
    def check_opens(*args, **kwargs):
        if args[0].fig is not None:
            if not plt.fignum_exists(args[0].fig.number):
                args[0].Close()
        return func(*args, **kwargs)
    return check_opens

def derive(func):
    """
    Decorator to calculate the derivative of a function.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        data = func(*args, **kwargs)
        if 'der' not in kwargs:
            return data
        if type(data) == tuple:
            data = list(data)
        if kwargs['der'] in ['r', 'radius']:
            r = args[0].cell.radius(args[0].ghost)
            if type(data) == list:
                return data
            elif len(r) in data.shape:
                return IDL_derivative(r, data, 'radius')
            else:
                return data
        elif kwargs['der'] in ['theta', 'th']:
            if args[0].dim < 1:
                return data
            theta = args[0].cell.theta(args[0].ghost)
            if type(data) == list:
                return data
            elif len(theta) in data.shape:
                return IDL_derivative(theta, data, 'theta')
            else:
                return data
        elif kwargs['der'] in ['phi', 'ph']:
            if args[0].dim < 2:
                return data
            phi = args[0].cell.phi(args[0].ghost)
            if type(data) == list:
                return data
            elif len(phi) in data.shape:
                return IDL_derivative(phi, data, 'phi')
            else:
                return data
        elif kwargs['der'] in ['t', 'time']:
            if type(data) == list:
                for i in range(1, len(data)):
                    if type(data[i]) != dict:
                        data[i] = IDL_derivative(data[0], data[i])
                return data
            else:
                data[:, 1:] = IDL_derivative(data[:, 0], data[:, 1:], axis=0)
                return data
        else:
            return data

    return wrapper

def smooth(func):
    """
    Decorator that add a smothing of the data.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        data = func(*args, **kwargs)
        if 'smooth' not in kwargs:
            return data
        if 'smooth_window' in kwargs:
            window_points = kwargs['smooth_window']
        else:
            window_points = 5
        if kwargs['smooth'] == 'gauss':
            x = np.linspace(-window_points // 2, window_points // 2,
                            window_points)
            window = np.exp(-(x**2) / (2 * 2 ** 2))
        else:
            window = np.ones(window_points) / window_points   
        if type(data) == tuple:
            data = list(data)
        if type(data) == list:
            for i in range(1, len(data)):
                if type(data[i]) != dict:
                    data[i] = np.convolve(data[i], window, mode='same')
            return data
        elif len(data.shape) == 1:
            return np.convolve(data, window, mode='same')
        elif data.shape[0] not in [len(args[0].cell.radius(args[0].ghost)),
                                   len(args[0].cell.theta(args[0].ghost)),
                                   len(args[0].cell.phi(args[0].ghost)),
                                   len(args[0].hdf_file_list)]:
            for i in range(data.shape[1]):
                data[:, i] = np.convolve(data[:, i], window, mode='same')
            return data
        elif len(args[0].hdf_file_list) == data.shape[0]:
            for i in range(data.shape[1]):
                data[:, i] = np.convolve(data[:, i], window, mode='same')
            return data
        else:
            if kwargs['smooth'] == 'gauss':
                grid = np.meshgrid(*[np.linspace(-window_points // 2,
                                                 window_points // 2,
                                                 window_points)
                                     for _ in range(data.ndim)],
                                   indexing='ij')
                coords = np.stack(grid, axis=-1)
                exponent = -np.sum((coords) ** 2, axis=-1) / \
                            (2 * 2 ** 2)
                window = np.exp(exponent)
                window /= np.sum(window)
            else:
                window = np.ones([window_points] * data.ndim)
                window /= np.sum(window)
            return convolve(data, window, mode='wrap')

    return wrapper

def EMD_smooth(func):
    """
    Decorator that add a smothing of the data for EMDs.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        data = func(*args, **kwargs)
        if 'smooth' not in kwargs:
            return data
        if 'smooth_window' in kwargs:
            window_points = kwargs['smooth_window']
        else:
            window_points = 5
        if kwargs['smooth'] == 'gauss':
            x = np.linspace(-window_points // 2, window_points // 2,
                            window_points)
            window = np.exp(-(x**2) / (2 * 2 ** 2))
        else:
            window = np.ones(window_points) / window_points   
        if data.ndim == 2:
            for i in range(0, data.shape[0]):
                data[i, :] = np.convolve(data[i, :], window, mode='same')
        else:
            data = np.convolve(data, window, mode='same')
        
        return data

    return wrapper