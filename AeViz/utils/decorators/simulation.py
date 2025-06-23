import h5py, os
from AeViz.utils.math_utils import IDL_derivative
from . import wraps, np, aerray, aeseries
import inspect
from AeViz.utils.files.string_utils import merge_strings
import warnings
from AeViz.units.aerray import apply_monkey_patch, remove_monkey_patch
try:
    from astropy.convolution import (convolve, Gaussian1DKernel,
                                     Gaussian2DKernel, Box1DKernel,
                                     Box2DKernel)
except:
    remove_monkey_patch()
    from astropy.convolution import (convolve, Gaussian1DKernel,
                                     Gaussian2DKernel, Box1DKernel,
                                     Box2DKernel)
    apply_monkey_patch()

def hdf_isopen(func):
    """
    Takes as input the Simulation object and either the file name, or
    file index or time. If the file is not open, it opens it.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        if type(args[1]) is int:
            file = args[0].hdf_file_list[args[1]]
        elif type(args[1]) is float or isinstance(args[1], aerray):
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

def derive(func):
    """
    Decorator to calculate the derivative of a function.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        data = func(*args, **kwargs)
        if 'der' not in kwargs and 'integrate' not in kwargs:
            return data
        if type(data) == tuple:
            data = list(data)
        if not type(data) == list:
            data = [data]
        for i, dd in enumerate(data):
            if isinstance(dd, aerray):
                if kwargs['der'] in ['time', 't']:
                    raise AttributeError('No time provided in aerray')
                elif kwargs['der'] in ['r', 'radius']:
                    r = args[0].cell.radius(args[0].ghost)
                    if not len(r) in dd.shape:
                        raise AttributeError(f'No dimension to derive. '\
                                             f'Radius has {len(r)} points '\
                                                f'and data has {dd.shape}')
                    data[i] = IDL_derivative(r, dd, axis=dd.shape.index(len(r)))
                elif kwargs['der'] in ['theta', 'th']:
                    if args[0].dim < 2:
                        raise AttributeError('Sorry, you do not have a '\
                                             'second dimension')
                    theta = args[0].cell.theta(args[0].ghost)
                    if not len(theta) in dd.shape:
                        raise AttributeError(f'No dimension to derive. '\
                                             f'Theta has {len(theta)} points '\
                                                f'and data has {dd.shape}')
                    data[i] = IDL_derivative(theta, dd,
                                        axis=dd.shape.index(len(theta)))
                elif kwargs['der'] in ['phi', 'ph']:
                    if args[0].dim < 3:
                        raise AttributeError('Sorry, you do not have a '\
                                             'third dimension')
                    phi = args[0].cell.phi(args[0].ghost)
                    if not len(phi) in dd.shape:
                        raise AttributeError(f'No dimension to derive. '\
                                             f'Phi has {len(phi)} points '\
                                                f'and data has {dd.shape}')
                    data[i] = IDL_derivative(phi, dd, axis=dd.shape.index(len(phi)))
            elif isinstance(dd, aeseries):
                if 'der' in kwargs:
                    data[i] = dd.derive(kwargs['der'])
                else:
                    data[i] = dd.integrate(kwargs['integrate'])
        if len(data) == 1:
            return data[0]
        return data
    return wrapper

def smooth(func):
    """
    Decorator that applies a convolution filter to the data, does not work with
    3D arrays.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        data = func(*args, **kwargs)
        if 'smooth' not in kwargs:
            return data
        if 'smooth_window' in kwargs:
            window_points = kwargs['smooth_window']
            if window_points % 2 == 0:
                window_points += 1
        else:
            window_points = 5
        if kwargs['smooth'] == 'gauss':
            smooth_type = 'gauss'
        else:
            smooth_type = 'box'

        if type(data) == tuple:
            data = list(data)
        if type(data) != list:
            data = [data]
        
                
        for i, dd in enumerate(data):
            try:
                dim = dd.ndim
            except:
                dim = dd.data.ndim
            if smooth_type == 'gauss':
                if dim == 1:
                    kernel = Gaussian1DKernel(2, x_size = window_points)
                elif dim == 2:
                    kernel = Gaussian2DKernel(2, 2, x_size=window_points,
                                              y_size=window_points)
                else:
                    kernel = np.ones(window_points ** dim).reshape(tuple([window_points] * dim))
                    kernel /= kernel.sum()
            elif smooth_type == 'box':
                if dim == 1:
                    kernel = Box1DKernel(window_points)
                elif dim == 2:
                    kernel = Box2DKernel(window_points)
                else:
                    kernel = np.ones(window_points ** dim).reshape(tuple([window_points] * dim))
                    kernel /= kernel.sum()
            if isinstance(dd, aerray):
                data[i] = aerray(convolve(dd.value, kernel, boundary='extend'),
                                 dd.unit, dd.name, dd.label, dd.cmap, dd.limits,
                                 dd.log)
            elif isinstance(dd, aeseries):
                data[i].data = aerray(convolve(dd.data.value, kernel, boundary='extend'),
                                      dd.data.unit, dd.data.name, dd.data.label,
                                      dd.data.cmap, dd.data.limits, dd.data.log)
            else:
                raise TypeError("Type not supported.")
        if len(data) == 1:
            data = data[0]
        return data
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
            if window_points % 2 == 0:
                window_points += 1
        else:
            window_points = 5
        if kwargs['smooth'] == 'gauss':
            smooth_type = 'gauss'
        else:
            smooth_type = 'box'
        if smooth_type == 'gauss':
            kernel = Gaussian1DKernel(2, x_size = window_points)
        elif smooth_type == 'box':
            kernel = Box1DKernel(window_points)
        if data.ndim == 2:
            for i in range(data.shape[1]):
                data[:, i] = convolve(data[:, i], kernel, boundary='extend')
        else:
            data = convolve(data, kernel, boundary='extend')
        return data

    return wrapper

def subtract_tob(func):
    """
    subtract the tob from an aeseries
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        sig = inspect.signature(func)
        bound_args = sig.bind_partial(*args, **kwargs)
        bound_args.apply_defaults()
        data = func(*args, **kwargs)
        if not bound_args.arguments['tob_corrected']:
            return data
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if type(data) == tuple:
                data = list(data)
            if type(data) == list:
                for dd in data:
                    nm, lb = dd.time.name, dd.time.label
                    dd.time -= args[0].tob
                    dd.time.set(name=nm,
                                label=merge_strings(lb,r'$-$', r'$t_\mathrm{b}$'),
                                limits=[-0.005, dd.time.value[-1]])
            else:
                nm, lb = data.time.name, data.time.label
                data.time -= args[0].tob
                data.time.set(name=nm,
                            label=merge_strings(lb,r'$-$', r'$t_\mathrm{b}$'),
                            limits=[-0.005, data.time.value[-1]])
        return data
    return wrapper

def sum_tob(func):
    """
    sum the tob from an aeseries
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        def handle_dict(a):
            for k, v in a.items():
                if isinstance(v, aeseries):
                    nm, lb = v.time.name, v.time.label
                    v.time += args[0].tob
                    v.time.set(name=nm,
                                label=merge_strings(lb,r'$+$', r'$t_\mathrm{b}$'),
                                limits=[0, v.time.value[-1]])
                else:
                    continue
            return a
        sig = inspect.signature(func)
        bound_args = sig.bind_partial(*args, **kwargs)
        bound_args.apply_defaults()
        data = func(*args, **kwargs)
        if bound_args.arguments['tob_corrected']:
            return data
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if type(data) == tuple:
                data = list(data)
            if type(data) == dict:
                data = handle_dict(data)
            elif type(data) == list:
                for dd in data:
                    if isinstance(dd, aeseries):
                        nm, lb = dd.time.name, dd.time.label
                        dd.time += args[0].tob
                        dd.time.set(name=nm,
                                    label=merge_strings(lb,r'$+$', r'$t_\mathrm{b}$'),
                                    limits=[0, dd.time.value[-1]])
                    elif isinstance(dd, dict):
                        dd = handle_dict(dd)
            else:
                nm, lb = data.time.name, data.time.label
                data.time += args[0].tob
                data.time.set(name=nm,
                            label=merge_strings(lb,r'$+$', r'$t_\mathrm{b}$'),
                            limits=[0, data.time.value[-1]])
        return data
    return wrapper

def mask_points(func):
    """
    Decorator to mask the quantity with the selected value.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        data = func(*args, **kwargs)
        if not 'mask' in kwargs:
            return data
        if 'unbound' in kwargs['mask']:
            mask = (args[0].MHD_energy(args[1]) + \
                2 * args[0].gravitational_energy(args[1]) > 0)
        else:
            raise NotImplementedError(f"{kwargs['mask']} not available as input.")
        
        if 'not' not in kwargs['mask']:
            mask = mask
        data[mask] = np.nan
        return data
    return wrapper