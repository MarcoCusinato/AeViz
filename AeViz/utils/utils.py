import sys, os, h5py
import numpy as np
from AeViz.utils.file_utils import save_hdf
from AeViz.units import u
from AeViz.units.aerray import aerray
import re

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
        data = h5py.File(os.path.join(simulation.storage_path, 'time.h5'), 'r')
        time_array = aerray(data['time'][...], u.s, 'time', r'$t$', None,
                            [None, None])
        data.close()
        if len(time_array) == len(simulation.hdf_file_list):
            return time_array
        else:
            start_time = len(time_array)
    else:
        start_time = 0
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
    save_hdf(os.path.join(simulation.storage_path, 'time.h5'),
             ['time'], [time_array.value])
    time_array.set('time', r'$t$', None, [None, None])
    return time_array
        
def merge_strings(*args):
    """
    Marges the strings keeping in mind that they can have latex sintax
    in them
    """
    if len(args) == 0:
        return None
    if len(args) == 1:
        return args
    assert all([type(ar) == str for ar in args]), "Can only be concatenating strings"
    out_string = r''
    for ar in args:
        if out_string.endswith('$') and ar.startswith('$'):
            out_string = out_string[:-1] + ar[1:]
        else:
            out_string += ar
    return out_string

def apply_symbol(latex_str: str, symbol: str = "\\tilde"):
    if latex_str is None:
        return None
    # Check if the string starts and ends with $
    in_math_mode = latex_str.startswith('$') and latex_str.endswith('$')
    
    # Remove surrounding $ if present
    core_str = latex_str[1:-1] if in_math_mode else latex_str
    
    # Find the first letter or word before an underscore
    match = re.match(r"([^_]+)(.*)", core_str)
    
    if not match:
        return latex_str  # Return unchanged if no valid match
    
    first_part, rest = match.groups()
    
    # Apply the symbol
    modified = f"{symbol}{{{first_part}}}{rest}"
    
    # Restore math mode if needed
    return f"${modified}$" if in_math_mode else modified
        
