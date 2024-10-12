import h5py, os
import matplotlib.pyplot as plt

def hdf_isopen(func):
    """
    Takes as input the Simulation object and either the file name, or
    file index or time. If the file is not open, it opens it.
    """
    def wrapper(*args):
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
        return func(*args)
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