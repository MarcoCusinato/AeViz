from scidata.par.read_parfile import get_indices_from_parfile
from scidata.quantities.quantities import SimulationAnalysis
import h5py, os

def check_file_to_load(path):
    if path.endswith('.hdf5') or path.endswith('.h5') or path.endswith('.hdf'):
        return 'hdf5'
    elif not os.path.exists(path):
        return 'sim'
    
class Data:
    def __init__(self):
        self.data = None
        self.is_open = self.is_open()

    def Load(self, path, simulation_path=None, dim=None):
        self.data_type = check_file_to_load(path)
        if self.data_type == 'hdf5':
            self.data = h5py.File(path, 'r')
            dir_path = os.path.dirname(os.path.abspath(path))
            try:
                self.hydroTHD_index, _ = get_indices_from_parfile('start.pars', os.path.join(dir_path, '../pars'))
            except:
                self.hydroTHD_index = None
                Warning.warn('No start.pars file found in the parent directory of the hdf5 file. No auto-detection of the indices. Please set them manually.')
        elif self.data_type == 'sim':
            self.data = SimulationAnalysis(path, dim, simulation_path)
    
    def is_open(self):
        if self.data is None:
            return False
        else:
            return True
    
    def Close(self):
        if self.is_open:
            if type(self.data) is h5py.File:
                self.data.close()
            elif type(self.data) is SimulationAnalysis:
                del self.data
            self.data = None
    

        
