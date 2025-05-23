import os

## CHECKPOINTS FOR COMPUTING LOCAL QUANTITIES
checkpoints = {
    1: False,
    2: 400,
    3: 20
}

def progressBar(count_value, total, suffix=''):
    """
    Display a progress bar in the terminal.
    """
    import sys
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
    from AeViz.units.aerray import aerray
    from AeViz.utils.files.file_utils import save_hdf
    from AeViz.units import u
    import numpy as np
    import h5py
    if check_existence(simulation, 'time.h5'):
        data = h5py.File(os.path.join(simulation.storage_path, 'time.h5'), 'r')
        time_array = aerray(data['time'][...], u.s, 'time', r'$t$', None,
                            [None, None])
        if 'processed' in data.keys():
            processed_hdf = data['processed'][...]
            processed_hdf = [ff.decode("utf-8") for ff in processed_hdf]
            data.close()
        else:
            processed_hdf = simulation.hdf_file_list[:len(time_array)]
            data.close()
            save_hdf(os.path.join(simulation.storage_path, 'time.h5'),
                    ['time', 'processed'], [time_array.value, processed_hdf])
        
        if processed_hdf[-1] == simulation.hdf_file_list[-1]:
            return time_array
        else:
            start_time = len(processed_hdf)
    else:
        start_time = 0
        processed_hdf = []
    progress_index = 0
    total_index = len(simulation.hdf_file_list[start_time:])
    for file_name in simulation.hdf_file_list[start_time:]:
        try:
            time_array = np.concatenate((time_array,
                                         simulation.time(file_name)))
        except:
            time_array = simulation.time(file_name)
        progressBar(progress_index, total_index, 'Storing timeseries')
        progress_index += 1
        processed_hdf.append(file_name)
    save_hdf(os.path.join(simulation.storage_path, 'time.h5'),
             ['time', 'processed'], [time_array.value, processed_hdf])
    time_array.set('time', r'$t$', None, [None, None])
    return time_array
