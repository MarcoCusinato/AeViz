from AeViz.cell.cell import cell
from AeViz.cell.ghost import ghost
from AeViz.units import u
from AeViz.utils.parfiles import (load_parfile,
                                  get_stencils,
                                  get_indices_from_parfile,
                                  get_simulation_info,
                                  get_initial_parameters)
from AeViz.utils.path_utils import (pltf, simulation_local_storage_folder, 
                                    find_simulation)
from AeViz.utils.decorators import hdf_isopen
from AeViz.utils.file_utils import list_module_functions
from AeViz.utils.utils import time_array
import numpy as np
import types, os


class Simulation:
    def __init__(self, simulation_name, simulation_folder_path=None,
                 dim = None):
        """
        Simulation class initialization. This is the class that should
        be called when loading up a simulation.
        Methods are loaded at runtime based on the simulation type.
        """
        # Let's set up the main parameters that define a simulation,
        # dimension, geometry, path and name
        # Also let's load the indices of the specific quantities
        self.simulation_name = simulation_name
        self.path = find_simulation(self.simulation_name, pltf(),
                                    simulation_folder_path)
        parfile = load_parfile('start.pars', os.path.join(self.path, 'pars'))
        self.GEOM, self.dim, self.relativistic, \
            self.evolved_qts = get_simulation_info(parfile)
        self.ghost_cells = get_stencils(parfile)
        self.hydroTHD_index = get_indices_from_parfile(parfile)
        self.cell = cell(self.path, self.dim)
        self.ghost = ghost(self.ghost_cells)
        self.storage_path = simulation_local_storage_folder(pltf(), 
                                                self.simulation_name, self.dim)
        del parfile
        if self.GEOM == 2:
            # if in spherical symmetry we suppose we are simulating a
            # SN, so try to load the information about the progenitor
            try:
                self.initial_parameters = get_initial_parameters(
                                                        os.path.join(self.path,
                                                                       'pars'))
            except:
                self.initial_parameters = None
        
        ## Private attributes
        ## Paths to the main simulation folders
        self.__log_path = os.path.join(self.path, 'log')
        self.__hdf_path = os.path.join(self.path, 'outp-hdf')
        self.__grid_path = os.path.join(self.path, 'grid')
        ## Name of the files containing global simulation data
        self.__integrated_nu_path = 'neu.dat'
        self.__rho_max_path = 'rho.dat'
        self.__grw_path = 'grw.dat'
        self.__mag_data = 'mag.dat'
        ## Opened file name
        self.__data_h5 = None
        self.__opened_hdf_file = ''
        self.hdf_file_list = self.__get_hdf_file_list()
        ## Load the methods based on the simulation type
        self.__load_hydro_methods()
        self.__load_THD_methods()
        self.__load_composition_methods()
        self.__load_global_methods()
        if self.evolved_qts['magdim'] > 0:
            self.__load_magnetic_fields_methods()
        if self.evolved_qts['neudim'] > 0:
            self.__load_neutrinos_methods()
        if self.GEOM == 2:
            self.__load_supernova_methods()
            self.__load_instabilities_methods()
            self.__load_spherical_harmonics_methods()
            self.__load_spherical_symmetry_methods()
            if self.dim > 1:
                self.__load_GWs_methods()
            self.tob = self.time_of_bounce()

    ## -----------------------------------------------------------------
    ## UTILITIES
    ## -----------------------------------------------------------------

    ## FILE LIST AND FILE SEARCH
    def __get_hdf_file_list(self):
        """
        List of all the 'timestep' files in the outp-hdf folder.
        """
        file_list = os.listdir(self.__hdf_path)
        #remove x00 files
        file_list = [x for x in file_list if x.startswith('h')]
        #return the sorted list of files
        file_list.sort()
        return file_list

    def find_file_from_time(self, time_to_find, time_in_ms=True, 
                            return_index=False, tob_corrected=True):
        """
        Returns the name of the file corresponding to the given time. If
        return_index is True, returns the index of the file in the 
        hdf_file_list. If time_in_ms is True, the time to giveis in ms,
        otherwise is in s.
        """
        if time_in_ms:
            time_to_find = u.convert_to_s(time_to_find)
        
        file_list = self.hdf_file_list
        time = time_array(self)
        if not tob_corrected and self.GEOM == 2:
            time += self.tob
        index = np.argmax(time>=time_to_find)
        if return_index:
            return file_list[index], index
        return file_list[index]
    
    ## ERROR
    @hdf_isopen
    def error(self, file_name, **kwargs):
        return self.ghost.remove_ghost_cells(np.squeeze(
            self.__data_h5['thd/data'][..., self.hydroTHD_index['thd']
                                          ['I_EOSERR']]), self.dim)
    
    ## TIME
    @hdf_isopen
    def time(self, file_name, tob_corrected=True):
        """
        Time of the simulation. If tob_corrected is True, the time is 
        corrected for the time of bounce.
        """
        if tob_corrected:
            return np.array(self.__data_h5['Parameters/t']) - self.tob
        return np.array(self.__data_h5['Parameters/t'])
    
    def __load_hydro_methods(self):
        import AeViz.simulation.methods.hydro
        funcs = list_module_functions(AeViz.simulation.methods.hydro)
        for name, obj in funcs:
            method_name = name
            if name == 'omega' and self.GEOM != 2 and self.dim < 2:
                continue
            if 'velocity' in name:
                if self.evolved_qts['veldim'] < 1:
                    continue
                if 'theta' in name and self.evolved_qts['veldim'] < 2:
                    continue
                if 'phi' in name and self.evolved_qts['veldim'] < 3:
                    continue
                if self.GEOM == 1:
                    method_name = method_name.replace('radial', 'x').\
                        replace('theta', 'y').replace('phi', 'z')
            setattr(self, method_name, types.MethodType(obj, self))
            
    def __load_THD_methods(self):
        import AeViz.simulation.methods.thermo
        funcs = list_module_functions(AeViz.simulation.methods.thermo)
        for name, obj in funcs:
            if name == 'entropy' and self.evolved_qts['entropie'] < 1:
                continue
            if name == 'temperature' and self.evolved_qts['tempratr'] < 1:
                continue
            if name == 'lorentz' and not self.relativistic:
                continue
            setattr(self, name, types.MethodType(obj, self))
    
    def __load_composition_methods(self):
        import AeViz.simulation.methods.composition
        funcs = list_module_functions(AeViz.simulation.methods.composition)
        for name, obj in funcs:
            if name == 'Ye' and self.evolved_qts['y_edim'] < 1:
                continue
            if self.evolved_qts['comp_dim'] < 1 and name in ['heavy_fraction',
                                                             'Abar', 'Zbar',
                                                             'proton_fraction',
                                                             'neutron_fraction',
                                                             'alpha_fraction']:
                continue
            if self.evolved_qts['cpot_dim'] < 1 and 'chemical_potential' in \
                name:
                continue
            setattr(self, name, types.MethodType(obj, self))
    
    def __load_neutrinos_methods(self):
        import AeViz.simulation.methods.neutrinos
        funcs = list_module_functions(AeViz.simulation.methods.neutrinos)
        for name, obj in funcs:
            setattr(self, name, types.MethodType(obj, self))
    
    def __load_magnetic_fields_methods(self):
        import AeViz.simulation.methods.magnetic_fields
        funcs = list_module_functions(AeViz.simulation.methods.magnetic_fields)
        for name, obj in funcs:
            if ('poloidal' in name or 'toroidal' in name) and self.GEOM != 2:
                continue
            setattr(self, name, types.MethodType(obj, self))
    
    def __load_GWs_methods(self):
        import AeViz.simulation.methods.GWs
        funcs = list_module_functions(AeViz.simulation.methods.GWs)
        for name, obj in funcs:
            setattr(self, name, types.MethodType(obj, self))
        try:
            import PyEMD
            import AeViz.simulation.methods.IMFs
            funcs = list_module_functions(AeViz.simulation.methods.IMFs)
            for name, obj in funcs:
                setattr(self, name, types.MethodType(obj, self))
        except:
            pass
    
    def __load_global_methods(self):
        import AeViz.simulation.methods.global_methods
        funcs = list_module_functions(AeViz.simulation.methods.global_methods)
        for name, obj in funcs:
            if 'neutrino' in name and self.evolved_qts['neudim'] < 1:
                continue
            if 'Ye' in name and self.evolved_qts['y_edim'] < 1:
                continue
            setattr(self, name, types.MethodType(obj, self))
    
    def __load_supernova_methods(self):
        import AeViz.simulation.methods.supernova
        funcs = list_module_functions(AeViz.simulation.methods.supernova)
        for name, obj in funcs:
            setattr(self, name, types.MethodType(obj, self))

    def __load_instabilities_methods(self):
        import AeViz.simulation.methods.instabilities
        funcs = list_module_functions(AeViz.simulation.methods.instabilities)
        for name, obj in funcs:
            setattr(self, name, types.MethodType(obj, self))

    def __load_spherical_harmonics_methods(self):
        import AeViz.simulation.methods.sph_harmonics
        funcs = list_module_functions(AeViz.simulation.methods.sph_harmonics)
        for name, obj in funcs:
            setattr(self, name, types.MethodType(obj, self))
    
    def __load_spherical_symmetry_methods(self):
        import AeViz.simulation.methods.other_spherical_sym
        funcs = list_module_functions(
            AeViz.simulation.methods.other_spherical_sym)
        for name, obj in funcs:
            setattr(self, name, types.MethodType(obj, self))