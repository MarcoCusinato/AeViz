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
        Simulation class initialization.
        """
        self.simulation_name = simulation_name
        self.path = find_simulation(self.simulation_name, pltf(), simulation_folder_path)
        self.log_path = os.path.join(self.path, 'log')
        self.hdf_path = os.path.join(self.path, 'outp-hdf')
        self.grid_path = os.path.join(self.path, 'grid')
        self.par_path = os.path.join(self.path, 'pars')
        self.integrated_nu = 'neu.dat'
        self.rho_max_file = 'rho.dat'
        self.grw = 'grw.dat'
        self.nu_e_grid = 'grid.e.dat'
        self.hydroTHD_index, self.ghost_cells = get_indices_from_parfile('start.pars', self.par_path)
        self.cell = cl(self.path, dim)
        self.dim = self.cell.simulation_dimension()
        self.ghost = gh(self.ghost_cells)
        self.storage_path = simulation_local_storage_folder(pltf(), self.simulation_name, self.dim)
        self.initial_parameters = get_initial_parameters(self.par_path)