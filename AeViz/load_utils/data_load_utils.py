from AeViz.utils.parfiles import get_indices_from_parfile
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
from AeViz.utils.path_utils import local_storage_folder

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
                if name == 'BX':
                    out = getattr(self.loaded_data,
                                  'magnetic_fields')(file, **kwargs)[..., 0]
                elif name == 'BY':
                    out = getattr(self.loaded_data, 'magnetic_fields')(file)\
                        [..., 1]
                elif name == 'BZ':
                    out = getattr(self.loaded_data,
                                  'magnetic_fields')(file, **kwargs)[..., 2]
                elif name == 'total_magnetic_energy':
                    out = getattr(self.loaded_data,
                                  'magnetic_energy')(file, **kwargs)[0]
                elif name == 'poloidal_magnetic_energy':
                    out = getattr(self.loaded_data,
                                  'magnetic_energy')(file, **kwargs)[1]
                elif name == 'toroidal_magnetic_energy':
                    out = getattr(self.loaded_data,
                                  'magnetic_energy')(file, **kwargs)[2]
                elif 'nu' in name and 'moment' in name:
                    out = return_neutrino_flux(self.loaded_data, name, file,
                                               **kwargs)
                elif 'nue_mean_ene' == name or 'nua_mean_ene' == name or \
                    'nux_mean_ene' == name:
                    out = return_neutrino_mean_ene(self.loaded_data,
                                                   name, file, **kwargs)
                else:
                    out = getattr(self.loaded_data, name)(file, **kwargs)
                return out
            else:
                if ('explosion' in name) or ('gain' in name) or \
                    ('innercore' in name) or \
                    ('PNS' in name and 'radius' not in name):
                    return self.__get_energy_data(name, **kwargs)
                elif 'nu_integrated_' in name:
                    return  return_integrated_neutrinos(self.loaded_data, name,
                                                        **kwargs)
                elif 'kick_velocity_' in name:
                    return return_PNS_kick(self.loaded_data, name)
                    
                return getattr(self.loaded_data, name)(**kwargs)
    
    def __plane_cut(self, data, indextheta = None, indexphi = None):
        if (indexphi == None and indextheta == None) or \
            self.sim_dim == 2:
            return data
        if indexphi is not None:
            return np.concatenate([np.flip(data[indexphi, :, :], axis=0),
                                   data[(indexphi + data.shape[0] // 2) % 
                                        data.shape[0], :, :]], axis=0)
        elif indextheta is not None:
            
            return data[:, indextheta, :]
            
        return data
    
    def __get_GW_decomposition_data(self, name):
        if self.sim_dim == 2:
            _, t, AE220, f_h, nuc_h, conv_h, out_h = self.loaded_data.AE220()
            return t, AE220, f_h, nuc_h, conv_h, out_h
        else:
            raise TypeError('Not implemented for 3D simulations.')
    
    def __get_energy_data(self, name, **kwargs):
        qt = name.split('_')[-1]
        if 'explosion' in name:
            data = getattr(self.loaded_data, 'explosion_mass_ene')(**kwargs)
            if qt == 'mass':
                return data[0], data[1]
            elif qt == 'ene':
                return data[0], data[2]
            elif qt == 'kin':
                return data[0], data[3]
            elif qt == 'mag':
                return data[0], data[4]
        elif 'gain' in name:
            data = getattr(self.loaded_data, 'gain_mass_nu_heat')(**kwargs)
            if qt == 'mass':
                return data[0], data[1]
            elif qt == 'ene':
                return data[0], data[2]
        elif 'innercore' in name:
            data = getattr(self.loaded_data, 'innercore_mass_ene')(**kwargs)
            if qt == 'mass':
                return data[0], data[1]
            elif qt == 'ene':
                return data[0], data[6]
            elif qt == 'kin':
                return data[0], data[2]
            elif qt == 'mag':
                return data[0], data[3]
            elif qt == 'rot':
                return data[0], data[4]
            elif qt == 'grav':
                return data[0], data[5]
            elif qt == 'T/W':
                return data[0], data[7]
        elif 'PNS' in name:
            if 'PNS_angular_mom' in name:
                return return_angular_momentum(self.loaded_data, name,
                                               **kwargs)
            data = getattr(self.loaded_data, 'PNS_mass_ene')(**kwargs)
            if qt == 'mass':
                return data[0], data[1]
            elif qt == 'ene':
                return data[0], data[6]
            elif qt == 'kin':
                return data[0], data[2]
            elif qt == 'mag':
                return data[0], data[3]
            elif qt == 'rot':
                return data[0], data[4]
            elif qt == 'grav':
                return data[0], data[5]
            elif qt == 'conv':
                return data[0], data[7]

    def __get_1D_radii_data(self, name, **kwargs):
        name = name.split('_')
        rad_type = name[-1]
        if not 'neutrino' in name:
            name = '_'.join(name[:-1])
            flavour = None
        else:
            flavour = name[-2]
            name = '_'.join(name[:-2])
        time, _, max_r, min_r, avg_r, _ = getattr(self.loaded_data,
                                                  name)(**kwargs)
        if flavour:
            max_r = max_r[flavour]
            min_r = min_r[flavour]
            avg_r = avg_r[flavour]
        if rad_type == 'max':
            return time, max_r
        elif rad_type == 'min':
            return time, min_r
        elif rad_type == 'avg':
            return time, avg_r
        else:
            return time, max_r, min_r, avg_r
        
