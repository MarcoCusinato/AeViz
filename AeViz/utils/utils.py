import sys, os, h5py
import numpy as np
from AeViz.utils.file_utils import save_hdf

def progressBar(count_value, total, suffix=''):
    bar_length = 100
    filled_up_Length = int(round(bar_length * count_value / float(total)))
    percentage = round(100.0 * count_value/float(total),1)
    bar = '#' * filled_up_Length + '-' * (bar_length - filled_up_Length)
    sys.stdout.write('[%s] %s%s ...%s\r' %(bar, percentage, '%', suffix))
    sys.stdout.flush()
    
def check_existence(simulation, file_name):
    """
    Check if the radius calculation has already been performed or at
    least partially performed.
    """
    if os.path.exists(os.path.join(simulation.storage_path, 
                                       file_name)):
        return True
    else:
        return False
    
def time_array(simulation):
    """
    Get the time array of the local simulation output.
    """
    if check_existence(simulation, 'time.h5'):
        data = h5py.File('time.h5', 'r')
        time_array = data['time'][:]
        data.close()
        if len(time_array) == len(simulation.hdf_file_list):
            return time_array
        else:
            start_time = len(time_array)
    else:
        start_time = 0
    
    for file_name in simulation.hdf_file_list[start_time:]:
        try:
            time_array = np.concatenate((time_array,
                                         simulation.time(file_name)))
        except:
            time_array = simulation.time(file_name)
    save_hdf('time.h5', ['time'], ['time'])
    return time_array
        
        