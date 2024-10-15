import h5py, os
import matplotlib.pyplot as plt
from AeViz.utils.math_utils import IDL_derivative

def hdf_isopen(func):
    """
    Takes as input the Simulation object and either the file name, or
    file index or time. If the file is not open, it opens it.
    """
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