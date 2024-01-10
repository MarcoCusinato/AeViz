from scidata.par.read_parfile import get_indices_from_parfile
from scidata.quantities.quantities import SimulationAnalysis
from scidata.units.units import units
import h5py, os
import numpy as np

u = units()

def check_file_to_load(path):
    if path.endswith('.hdf5') or path.endswith('.h5') or path.endswith('.hdf'):
        return 'hdf5'
    elif not os.path.exists(path):
        return 'sim'

def return_index(hydro_dict, name):
    convert_dict = {
        'RHO': {'type': 'hydro', 'index': hydro_dict['hydro']['I_RH']},
        'ENE': {'type': 'hydro', 'index': hydro_dict['hydro']['I_EN']},
        'VX': {'type': 'hydro', 'index': hydro_dict['hydro']['I_VX']},
        'VY': {'type': 'hydro', 'index': hydro_dict['hydro']['I_VY']},
        'VZ': {'type': 'hydro', 'index': hydro_dict['hydro']['I_VZ']},
        'YE': {'type': 'hydro', 'index': hydro_dict['hydro']['I_YE']},
        #'YZ': {'type': 'hydro', 'index': hydro_dict['hydro']['I_YZ']},
        'ERR': {'type': 'thd', 'index': hydro_dict['thd']['I_EOSERR']},
        'LRTZ': {'type': 'thd', 'index': hydro_dict['thd']['I_LRTZ']},
        'DENS': {'type': 'thd', 'index': hydro_dict['thd']['I_DENS']},
        'EINT': {'type': 'thd', 'index': hydro_dict['thd']['I_EINT']},
        'ENTH': {'type': 'thd', 'index': hydro_dict['thd']['I_ENTH']},
        'PELE': {'type': 'thd', 'index': hydro_dict['thd']['I_PELE']},
        'TELE': {'type': 'thd', 'index': hydro_dict['thd']['I_TELE']},
        'NELE': {'type': 'thd', 'index': hydro_dict['thd']['I_NELE']},
        'PION': {'type': 'thd', 'index': hydro_dict['thd']['I_PION']},
        'TION': {'type': 'thd', 'index': hydro_dict['thd']['I_TION']},
        'NION': {'type': 'thd', 'index': hydro_dict['thd']['I_NION']},
        'VELX': {'type': 'thd', 'index': hydro_dict['thd']['I_VELX']},
        'VELY': {'type': 'thd', 'index': hydro_dict['thd']['I_VELY']},
        'VELZ': {'type': 'thd', 'index': hydro_dict['thd']['I_VELZ']},
        'T': {'type': 'thd', 'index': hydro_dict['thd']['I_TMPR']},
        'ENTR': {'type': 'thd', 'index': hydro_dict['thd']['I_ENTR']},
        'GAMM': {'type': 'thd', 'index': hydro_dict['thd']['I_GAMM']},
        'HEAT': {'type': 'thd', 'index': hydro_dict['thd']['I_HEAT']},
        'DELP': {'type': 'thd', 'index': hydro_dict['thd']['I_DELP']},
        'SMOMX': {'type': 'thd', 'index': hydro_dict['thd']['I_SMOMX']},
        'SMOMY': {'type': 'thd', 'index': hydro_dict['thd']['I_SMOMY']},
        'SMOMZ': {'type': 'thd', 'index': hydro_dict['thd']['I_SMOMZ']},
        'PGAS': {'type': 'thd', 'index': hydro_dict['thd']['I_PGAS']},
        'VSOUND': {'type': 'thd', 'index': hydro_dict['thd']['I_CSND']},
        'X_n': {'type': 'thd', 'index': hydro_dict['thd']['I_COMP'][0]},
        'X_p': {'type': 'thd', 'index': hydro_dict['thd']['I_COMP'][1]},
        'X_alpha': {'type': 'thd', 'index': hydro_dict['thd']['I_COMP'][2]},
        'X_h': {'type': 'thd', 'index': hydro_dict['thd']['I_COMP'][3]},
        'Abar': {'type': 'thd', 'index': hydro_dict['thd']['I_COMP'][4]},
        'Zbar': {'type': 'thd', 'index': hydro_dict['thd']['I_COMP'][5]},
        'CPOT_e': {'type': 'thd', 'index': hydro_dict['thd']['I_CPOT'][0]},
        'CPOT_n': {'type': 'thd', 'index': hydro_dict['thd']['I_CPOT'][1]},
        'CPOT_p': {'type': 'thd', 'index': hydro_dict['thd']['I_CPOT'][2]},
        'CPOT_nu': {'type': 'thd', 'index': hydro_dict['thd']['I_CPOT'][3]},
        'BHEX': {'type': 'thd', 'index': hydro_dict['thd']['I_BHEX']}
        }
    return convert_dict[name]
    
class Data:
    def __init__(self):
        self.loaded_data = None
        self.is_open = self.is_open()
        self.hydroTHD_index, self.gh_cells = None, None
        self.radius, self.theta, self.phi = None, None, None
        self.sim_dim = None

    def Load(self, path, simulation_path=None, dim=None):
        self.data_type = check_file_to_load(path)
        if self.data_type == 'hdf5':
            self.loaded_data = h5py.File(path, 'r')
            dir_path = os.path.dirname(os.path.abspath(path))
            self.radius = u.convert_to_km(self.loaded_data['X']['znc'][...])
            self.theta = self.loaded_data['Y']['znc'][...]
            self.phi = self.loaded_data['Z']['znc'][...]
            if self.theta.size == 1:
                self.sim_dim = 1
            elif self.phi.size == 1:
                self.sim_dim = 2
            else:
                self.sim_dim = 3
            try:
                self.hydroTHD_index, self.gh_cells = get_indices_from_parfile('start.pars', os.path.join(dir_path, '../pars'))
            except:
                self.hydroTHD_index = None
                Warning.warn('No start.pars file found in the parent directory of the hdf5 file. No auto-detection of the indices. Please set them manually.')
        elif self.data_type == 'sim':
            self.loaded_data = SimulationAnalysis(path, dim, simulation_path)
    
    def is_open(self):
        if self.loaded_data is None:
            return False
        else:
            return True
    
    def Close(self):
        if self.is_open:
            if type(self.loaded_data) is h5py.File:
                self.loaded_data.close()
            elif type(self.loaded_data) is SimulationAnalysis:
                del self.loaded_data
            self.loaded_data = None
    
    def __get_data_from_name(self, name):
        if self.data_type == 'hdf5':
            data = np.squeeze(np.array(self.loaded_data[return_index(self.hydroTHD_index, name)['type']]['data'])[..., return_index(self.hydroTHD_index, name)['index']])
            if name == 'YE':
                data /= np.squeeze(np.array(self.loaded_data['hydro']['data'])[:,:,:,self.hydroTHD_index['hydro']['I_RH']])
            return data
        elif self.data_type == 'sim':
            raise NotImplementedError('Not implemented yet.')
    

        
