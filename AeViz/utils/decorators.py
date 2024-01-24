import h5py, os

def hdf_isopen(func):
    def wrapper(*args):
        if args[1] != args[0]._Simulation__opened_hdf_file:
            if args[0]._Simulation__data_h5 is not None:
                args[0]._Simulation__data_h5.close()
            if args[1] not in args[0].hdf_file_list:
                raise ValueError("Selected file does not exist.")
            args[0]._Simulation__opened_hdf_file = args[1]
            args[0]._Simulation__data_h5 = h5py.File(os.path.join(args[0]._Simulation__hdf_path, args[1]), 'r')
        return func(*args)
    return wrapper