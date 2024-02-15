import h5py, os

def hdf_isopen(func):
    """
    Takse as input the Simulation object and either the file name, or
    file index or time. If the file is not open, it opens it.
    """
    def wrapper(*args):
        if type(args[1]) is int:
            args[1] = args[0].hdf_file_list[args[1]]
        elif type(args[1]) is float:
            args[1] = args[0].find_file_from_time(args[1])
        if args[1] != args[0]._Simulation__opened_hdf_file:
            if args[0]._Simulation__data_h5 is not None:
                args[0]._Simulation__data_h5.close()
            if args[1] not in args[0].hdf_file_list:
                raise ValueError("Selected file does not exist.")
            args[0]._Simulation__opened_hdf_file = args[1]
            args[0]._Simulation__data_h5 = h5py.File(
                os.path.join(args[0]._Simulation__hdf_path, args[1]), 'r')
        return func(*args)
    return wrapper