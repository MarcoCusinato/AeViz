from AeViz.utils.parfiles import get_indices_from_parfile
from AeViz.simulation.simulation import Simulation
from AeViz.units.units import units
import h5py, os
import numpy as np
from AeViz.cell.cell import cell as cl
from AeViz.cell.ghost import ghost
from AeViz.load_utils.utils import check_file_to_load, return_index

u = units()
    
class Data(object):
    """
    Class to load either an hdf5 file or an entire simulation.
    """
    def __init__(self):
        self.loaded_data = None
        self.hydroTHD_index, self.gh_cells = None, None
        self.sim_dim = None

    def Load(self, path, simulation_path=None, dim=None):
        """
        Load into the memory either an hdf5 file or an entire simulation.
        """
        self.data_type = check_file_to_load(path)
        if self.data_type == 'hdf5':
            self.__load_hdf(path)
        elif self.data_type == 'sim':
            self.loaded_data = Simulation(path, simulation_path, dim)
            self.cell = self.loaded_data.cell
            self.ghost = self.loaded_data.ghost
            self.sim_dim =self.loaded_data.dim
    
    def __load_hdf(self, path):
        """
        Load an hdf5 file into the memory using its path. 
        """
        self.loaded_data = h5py.File(path, 'r')
        dir_path = os.path.dirname(os.path.abspath(path))
        self.cell = cl(radius=np.stack((self.loaded_data['X']['znl'][...],
                                        self.loaded_data['X']['znc'][...],
                                        self.loaded_data['X']['znr'][...]),
                                        axis=-1),
                  theta=np.stack((self.loaded_data['Y']['znl'][...],
                                  self.loaded_data['Y']['znc'][...],
                                  self.loaded_data['Y']['znr'][...]), axis=-1),
                  phi=np.stack((self.loaded_data['Z']['znl'][...],
                                self.loaded_data['Z']['znc'][...],
                                self.loaded_data['Z']['znr'][...]), axis=-1))
        self.sim_dim = self.cell.dim
        try:
            self.hydroTHD_index, self.gh_cells = get_indices_from_parfile(
                'start.pars', os.path.join(dir_path, '../pars'))
            self.ghost = ghost(self.gh_cells)
        except:
            self.hydroTHD_index = None
            Warning.warn('No start.pars file found in the parent directory '\
                        ' of the hdf5 file. No auto-detection of the indices.'\
                            ' Please set them manually.')
    

    def is_open(self):
        if self.loaded_data is None:
            return False
        else:
            return True
    
    def Close(self):
        if self.is_open():
            if type(self.loaded_data) is h5py.File:
                self.loaded_data.close()
            elif type(self.loaded_data) is Simulation:
                del self.loaded_data
            self.loaded_data = None
    
    def __get_profile(self, name):
        if self.data_type == 'hdf5':
            raise ValueError('The profile method is not implemented for hdf5' \
                             ' files.')
        elif self.data_type == 'sim':
            return self.loaded_data.radial_profile(name)

    def __get_data_from_name(self, name, file=None):
        if self.data_type == 'hdf5':
            if return_index(self.hydroTHD_index, name)['index'] is None:
                raise ValueError('The selected quantity is not present in the'\
                                 ' hdf5 file.')
            index = return_index(self.hydroTHD_index, name)['index']
            if type(index) == list:
                data = np.squeeze(self.loaded_data[return_index(
                    self.hydroTHD_index, name)['type']][..., index[0],
                                                         index[1]])
            else:
                data = np.squeeze(self.loaded_data[return_index(
                    self.hydroTHD_index, name)['type']][..., index])
            if name in ['YE', 'YN', 'VX', 'VY', 'VZ']:
                data /= np.squeeze(self.loaded_data['hydro']['data']\
                                   [..., self.hydroTHD_index['hydro']['I_RH']])
            return self.ghost.remove_ghost_cells(data, self.sim_dim)
        elif self.data_type == 'sim':
            if file is not None:
                return getattr(self.loaded_data, name)(file)
            else:
                return getattr(self.loaded_data, name)()
    

        
