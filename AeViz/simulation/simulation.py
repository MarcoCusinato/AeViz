import numpy as np
import os
from AeViz.utils.parfiles import get_indices_from_parfile, get_initial_parameters
from AeViz.utils.path_utils import pltf, simulation_local_storage_folder, find_simulation
from AeViz.utils.file_utils import load_file
from AeViz.utils.decorators import hdf_isopen
from AeViz.cell.cell import cell as cl
from AeViz.cell.ghost import ghost as gh

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
        self.__data_h5 = None
        ## Opened file name
        self.__opened_hdf_file = ''
        self.hdf_file_list = self.__get_hdf_file_list()

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
    
    ## -----------------------------------------------------------------------------------------------------------------
    ## HYDRODYNAMICAL DATA
    ## -----------------------------------------------------------------------------------------------------------------
    
    @hdf_isopen
    def rho(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['hydro/data']) \
                                                        [..., self.hydroTHD_index['hydro']['I_RH']]), self.dim)
    
    @hdf_isopen
    def energy_MHD(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['hydro/data']) \
                                                        [..., self.hydroTHD_index['hydro']['I_EN']]), self.dim)
    
    ## VELOCITY
    @hdf_isopen
    def radial_velocity(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['thd/data']) \
                                                        [..., self.hydroTHD_index['thd']['I_VELX']]), self.dim)
    
    @hdf_isopen
    def theta_velocity(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['thd/data']) \
                                                        [..., self.hydroTHD_index['thd']['I_VELY']]), self.dim)
    
    @hdf_isopen
    def phi_velocity(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['thd/data']) \
                                                        [..., self.hydroTHD_index['thd']['I_VELZ']]), self.dim)
    
    @hdf_isopen
    def soundspeed(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['thd/data']) \
                                                        [..., self.hydroTHD_index['thd']['I_CSND']]), self.dim)
    ## -----------------------------------------------------------------------------------------------------------------
    ## THERMODYNAMICAL DATA
    ## -----------------------------------------------------------------------------------------------------------------

    ## THERMODYNAMICAL
    @hdf_isopen
    def gas_pressure(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['thd/data']) \
                                                        [..., self.hydroTHD_index['thd']['I_PGAS']]), self.dim)
    
    @hdf_isopen
    def temperature(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['thd/data']) \
                                                        [..., self.hydroTHD_index['thd']['I_TMPR']]), self.dim)
    
    @hdf_isopen
    def enthalpy(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['thd/data']) \
                                                        [..., self.hydroTHD_index['thd']['I_ENTH']]), self.dim)
    
    @hdf_isopen
    def entropy(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['thd/data']) \
                                                        [..., self.hydroTHD_index['thd']['I_ENTR']]), self.dim)
    
    @hdf_isopen
    def adiabatic_index(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['thd/data']) \
                                                        [..., self.hydroTHD_index['thd']['I_GAMM']]), self.dim)
    
    ## RELATIVISTIC
    @hdf_isopen
    def lorentz(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['thd/data']) \
                                                        [..., self.hydroTHD_index['thd']['I_LRTZ']]), self.dim)
    
    ## ENERGY
    @hdf_isopen
    def internal_energy(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['thd/data']) \
                                                        [..., self.hydroTHD_index['thd']['I_EINT']]), self.dim)
    
    @hdf_isopen
    def nu_heat(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['thd/data']) \
                                                        [..., self.hydroTHD_index['hydro']['I_HEAT']]), self.dim)
    
    ## COMPOSITION
    @hdf_isopen
    def Ye(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['hydro/data']) \
                                                        [..., self.hydroTHD_index['hydro']['I_YE']]), self.dim) / \
                self.rho(file_name)
    
    @hdf_isopen
    def neutron_fraction(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['thd/data']) \
                                                        [..., self.hydroTHD_index['thd']['I_COMP'][0]]), self.dim)
    
    @hdf_isopen
    def proton_fraction(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['thd/data']) \
                                                        [..., self.hydroTHD_index['thd']['I_COMP'][1]]), self.dim)
    
    @hdf_isopen
    def alpha_fraction(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['thd/data']) \
                                                        [..., self.hydroTHD_index['thd']['I_COMP'][2]]), self.dim)
    
    @hdf_isopen
    def heavy_fraction(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['thd/data']) \
                                                        [..., self.hydroTHD_index['thd']['I_COMP'][3]]), self.dim)
    
    @hdf_isopen
    def Abar(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['thd/data']) \
                                                        [..., self.hydroTHD_index['thd']['I_COMP'][4]]), self.dim)
    
    @hdf_isopen
    def Zbar(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['thd/data']) \
                                                        [..., self.hydroTHD_index['thd']['I_COMP'][5]]), self.dim)
    
    ## CHEMICAL POTENTIAL
    @hdf_isopen
    def electron_chemical_potential(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['thd/data']) \
                                                        [..., self.hydroTHD_index['thd']['I_CPOT'][0]]), self.dim)
    
    @hdf_isopen
    def neutron_chemical_potential(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['thd/data']) \
                                                        [..., self.hydroTHD_index['thd']['I_CPOT'][1]]), self.dim)
    
    @hdf_isopen
    def proton_chemical_potential(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['thd/data']) \
                                                        [..., self.hydroTHD_index['thd']['I_CPOT'][2]]), self.dim)
    
    @hdf_isopen
    def neutrino_chemical_potential(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(self.__data_h5['thd/data']) \
                                                        [..., self.hydroTHD_index['thd']['I_CPOT'][3]]), self.dim)
    
    