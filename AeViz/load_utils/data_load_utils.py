from AeViz.utils.files.parfiles import (get_indices_from_parfile,
                                  get_stencils,
                                  get_simulation_info)
from AeViz.simulation.simulation import Simulation
from AeViz.units.units import units
import h5py, os
import numpy as np
from AeViz.cell.cell import cell as cl
from AeViz.cell.ghost import ghost
from AeViz.load_utils.utils import (check_file_to_load, return_index,
                                    return_neutrino_flux,
                                    return_neutrino_mean_ene,
                                    return_integrated_neutrinos,
                                    return_angular_momentum,
                                    return_PNS_kick)
from AeViz.utils.files.path_utils import local_storage_folder
from AeViz.grid.grid import grid

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
            self.sim_dim = self.loaded_data.dim
        self.save_path = self.__save_path()
    
    def __save_path(self):
        if self.data_type == 'hdf5':
            return os.path.join(local_storage_folder(), 'hdf_plots')
        elif self.data_type == 'sim':
            return self.loaded_data.storage_path

    def __load_hdf(self, path):
        """
        Load an hdf5 file into the memory using its path. 
        """
        self.loaded_data = h5py.File(path, 'r')
        dir_path = os.path.dirname(os.path.abspath(path))
        try:
            startfile = os.path.join(dir_path, '../pars/start.pars')
            self.hydroTHD_index = get_indices_from_parfile(startfile)
            self.gh_cells = get_stencils(startfile)
            geom, dim, _, qts = get_simulation_info(startfile)
            self.ghost = ghost(self.gh_cells)
        except:
            raise ValueError('No start.pars file detected in the parent '\
                             'directory of the hdf5 file.')

        self.cell = cl(radius=np.stack((self.loaded_data['X']['znl'][...],
                                        self.loaded_data['X']['znc'][...],
                                        self.loaded_data['X']['znr'][...]),
                                        axis=-1),
                  theta=np.stack((self.loaded_data['Y']['znl'][...],
                                  self.loaded_data['Y']['znc'][...],
                                  self.loaded_data['Y']['znr'][...]), axis=-1),
                  phi=np.stack((self.loaded_data['Z']['znl'][...],
                                self.loaded_data['Z']['znc'][...],
                                self.loaded_data['Z']['znr'][...]), axis=-1),
                dim=dim, geom=geom, neu=qts['neudim'])
        self.sim_dim = self.cell.dim

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
    
    def __get_profile(self, name, **kwargs):
        if self.data_type == 'hdf5':
            raise ValueError('The profile method is not implemented for hdf5' \
                             ' files.')
        elif self.data_type == 'sim':
            if name == 'rho_spherical_harmonics':
                return self.loaded_data.rho_spherical_harmonics(**kwargs)
            elif name == 'fourier_amplitude':
                return self.loaded_data.fourier_amplitude(**kwargs)
            else:
                return self.loaded_data.radial_profile(name, **kwargs)

    def __get_data_from_name(self, name, file=None, **kwargs):
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
                return getattr(self.loaded_data, name)(file, **kwargs)
            return getattr(self.loaded_data, name)(**kwargs)