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
        'RHO': {'type': 'hydro/data', 'index': hydro_dict['hydro']['I_RH']},
        'ENE': {'type': 'hydro/data', 'index': hydro_dict['hydro']['I_EN']},
        'VX': {'type': 'hydro/data', 'index': hydro_dict['hydro']['I_VX']},
        'VY': {'type': 'hydro/data', 'index': hydro_dict['hydro']['I_VY']},
        'VZ': {'type': 'hydro/data', 'index': hydro_dict['hydro']['I_VZ']},
        'YE': {'type': 'hydro/data', 'index': hydro_dict['hydro']['I_YE']},
        'YN': {'type': 'hydro/data', 'index': hydro_dict['hydro']['I_YN']},
        'ERR': {'type': 'thd/data', 'index': hydro_dict['thd']['I_EOSERR']},
        'LRTZ': {'type': 'thd/data', 'index': hydro_dict['thd']['I_LRTZ']},
        'DENS': {'type': 'thd/data', 'index': hydro_dict['thd']['I_DENS']},
        'EINT': {'type': 'thd/data', 'index': hydro_dict['thd']['I_EINT']},
        'ENTH': {'type': 'thd/data', 'index': hydro_dict['thd']['I_ENTH']},
        'PELE': {'type': 'thd/data', 'index': hydro_dict['thd']['I_PELE']},
        'TELE': {'type': 'thd/data', 'index': hydro_dict['thd']['I_TELE']},
        'NELE': {'type': 'thd/data', 'index': hydro_dict['thd']['I_NELE']},
        'PION': {'type': 'thd/data', 'index': hydro_dict['thd']['I_PION']},
        'TION': {'type': 'thd/data', 'index': hydro_dict['thd']['I_TION']},
        'NION': {'type': 'thd/data', 'index': hydro_dict['thd']['I_NION']},
        'VELX': {'type': 'thd/data', 'index': hydro_dict['thd']['I_VELX']},
        'VELY': {'type': 'thd/data', 'index': hydro_dict['thd']['I_VELY']},
        'VELZ': {'type': 'thd/data', 'index': hydro_dict['thd']['I_VELZ']},
        'T': {'type': 'thd/data', 'index': hydro_dict['thd']['I_TMPR']},
        'ENTR': {'type': 'thd/data', 'index': hydro_dict['thd']['I_ENTR']},
        'GAMM': {'type': 'thd/data', 'index': hydro_dict['thd']['I_GAMM']},
        'HEAT': {'type': 'thd/data', 'index': hydro_dict['thd']['I_HEAT']},
        'DELP': {'type': 'thd/data', 'index': hydro_dict['thd']['I_DELP']},
        'JX': {'type': 'thd/data', 'index': hydro_dict['thd']['I_SMOMX']},
        'JY': {'type': 'thd/data', 'index': hydro_dict['thd']['I_SMOMY']},
        'JZ': {'type': 'thd/data', 'index': hydro_dict['thd']['I_SMOMZ']},
        'PGAS': {'type': 'thd/data', 'index': hydro_dict['thd']['I_PGAS']},
        'VSOUND': {'type': 'thd/data', 'index': hydro_dict['thd']['I_CSND']},
        'X_n': {'type': 'thd/data', 'index': hydro_dict['thd']['I_COMP'][0]},
        'X_p': {'type': 'thd/data', 'index': hydro_dict['thd']['I_COMP'][1]},
        'X_alpha': {'type': 'thd/data', 'index': hydro_dict['thd']['I_COMP'][2]},
        'X_h': {'type': 'thd/data', 'index': hydro_dict['thd']['I_COMP'][3]},
        'Abar': {'type': 'thd/data', 'index': hydro_dict['thd']['I_COMP'][4]},
        'Zbar': {'type': 'thd/data', 'index': hydro_dict['thd']['I_COMP'][5]},
        'CPOT_e': {'type': 'thd/data', 'index': hydro_dict['thd']['I_CPOT'][0]},
        'CPOT_n': {'type': 'thd/data', 'index': hydro_dict['thd']['I_CPOT'][1]},
        'CPOT_p': {'type': 'thd/data', 'index': hydro_dict['thd']['I_CPOT'][2]},
        'CPOT_nu': {'type': 'thd/data', 'index': hydro_dict['thd']['I_CPOT'][3]},
        'BHEX': {'type': 'thd/data', 'index': hydro_dict['thd']['I_BHEX']},
        'NUEE': {'type': 'neutrinogrey/egrey', 'index': [0, 0] },
        'NUEX': {'type': 'neutrinogrey/egrey', 'index': [0, 1] },
        'NUEY': {'type': 'neutrinogrey/egrey', 'index': [0, 2] },
        'NUEZ': {'type': 'neutrinogrey/egrey', 'index': [0, 3] },
        'NUAE': {'type': 'neutrinogrey/egrey', 'index': [1, 0] },
        'NUAX': {'type': 'neutrinogrey/egrey', 'index': [1, 1] },
        'NUAY': {'type': 'neutrinogrey/egrey', 'index': [1, 2] },
        'NUAZ': {'type': 'neutrinogrey/egrey', 'index': [1, 3] },
        'NUXE': {'type': 'neutrinogrey/egrey', 'index': [2, 0] },
        'NUXX': {'type': 'neutrinogrey/egrey', 'index': [2, 1] },
        'NUXY': {'type': 'neutrinogrey/egrey', 'index': [2, 2] },
        'NUXZ': {'type': 'neutrinogrey/egrey', 'index': [2, 3] },
        'BX': {'type': 'mag_vol/data', 'index': 0},
        'BY': {'type': 'mag_vol/data','index': 1},
        'BZ': {'type': 'mag_vol/data', 'index': 2},
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
            if return_index(self.hydroTHD_index, name)['index'] is None:
                raise ValueError('The selected quantity is not present in the hdf5 file.')
            index = return_index(self.hydroTHD_index, name)['index']
            if type(index) == list:
                data = np.squeeze(np.array(self.loaded_data[return_index(self.hydroTHD_index, name)['type']])[..., index[0], index[1]])
            else:
                data = np.squeeze(np.array(self.loaded_data[return_index(self.hydroTHD_index, name)['type']])[..., index])
            if name in ['YE', 'YN', 'VX', 'VY', 'VZ']:
                data /= np.squeeze(np.array(self.loaded_data['hydro']['data'])[..., self.hydroTHD_index['hydro']['I_RH']])
            return data
        elif self.data_type == 'sim':
            raise NotImplementedError('Not implemented yet.')
    

        
