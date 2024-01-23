import numpy as np
import h5py, os
from utils.parfiles import get_indices_from_parfile, get_initial_parameters
from utils.path_utils import pltf, simulation_local_storage_folder, find_simulation
from utils.file_utils import load_file
from cell.cell import cell as cl
from cell.ghost import ghost as gh



class Simulation:
    def __init__(self, simulation_name, simulation_folder_path=None, dim = None):
        """
        Simulation class initialization. Attibutes:

        """
        ## Public attributes
        self.simulation_name = simulation_name
        self.path = find_simulation(self.simulation_name, pltf(), simulation_folder_path)
        self.hydroTHD_index, self.ghost_cells = get_indices_from_parfile('start.pars',
                                                        os.path.join(self.path, 'pars'))
        self.initial_parameters = get_initial_parameters(os.path.join(self.path, 'pars'))
        self.cell = cl(self.path, dim)
        self.dim = self.cell.simulation_dimension()
        self.ghost = gh(self.ghost_cells)
        self.storage_path = simulation_local_storage_folder(pltf(), self.simulation_name, self.dim)
        self.hdf_file_list = self.__get_hdf_file_list()
        ## Private attributes
        ## Paths to the main simulation folders
        self.__log_path = os.path.join(self.path, 'log')
        self.__hdf_path = os.path.join(self.path, 'outp-hdf')
        self.__grid_path = os.path.join(self.path, 'grid')
        ## Name of the files containing global simulation data
        self.__integrated_nu_path = 'neu.dat'
        self.__rho_max_path = 'rho.dat'
        self.__grw_path = 'grw.dat'
        self.__nu_E_grid_path = 'grid.e.dat'
        ## Opened file name
        self.__opened_h5_file = ''

    def __get_hdf_file_list(self):
        """
        List of all the 'timestep' files in the outp-hdf folder.
        """
        file_list = os.listdir(self.hdf_path)
        #remove x00 files
        file_list = [x for x in file_list if x.startswith('h')]
        #return the sorted list of files
        file_list.sort()
        return file_list
    
