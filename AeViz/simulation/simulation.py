import numpy as np
from typing import Literal
import os, warnings
from AeViz.utils.parfiles import (get_indices_from_parfile,
                                  get_initial_parameters)
from AeViz.utils.path_utils import (pltf, simulation_local_storage_folder, 
                                    find_simulation)
from AeViz.utils.file_utils import load_file
from AeViz.utils.decorators import hdf_isopen
from AeViz.utils.math_utils import strfct2D
from AeViz.utils.file_utils import load_file, find_column_changing_line
from AeViz.utils.GW_utils import (GW_strain, calculate_AE220, GW_spectrogram,
                                  GWs_peak_indices, GWs_fourier_transform,
                                  GWs_frequency_peak_indices)
from AeViz.cell.cell import cell as cl
from AeViz.cell.ghost import ghost as gh
from AeViz.units.units import units

u = units()

class Simulation:
    def __init__(self, simulation_name, simulation_folder_path=None,
                 dim = None):
        """
        Simulation class initialization. Attibutes:

        """
        ## Public attributes
        self.simulation_name = simulation_name
        self.path = find_simulation(self.simulation_name, pltf(),
                                    simulation_folder_path)
        self.hydroTHD_index, \
            self.ghost_cells = get_indices_from_parfile('start.pars',
                                               os.path.join(self.path, 'pars'))
        self.initial_parameters = get_initial_parameters(os.path.join(self.path,
                                                                       'pars'))
        self.cell = cl(self.path, dim)
        self.dim = self.cell.simulation_dimension()
        self.ghost = gh(self.ghost_cells)
        self.storage_path = simulation_local_storage_folder(pltf(), 
                                                self.simulation_name, self.dim)
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
        self.__data_h5 = None
        ## Opened file name
        self.__opened_hdf_file = ''
        self.hdf_file_list = self.__get_hdf_file_list()
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

    def find_file_from_time(self, time_in_ms=True, return_index=False,
                            tob_corrected=True):
        """
        Returns the name of the file corresponding to the given time. If
        return_index is True, returns the index of the file in the 
        hdf_file_list. If time_in_ms is True, the time to giveis in ms,
        otherwise is in s.
        """
        if time_in_ms:
            time_to_find = u.convert_to_s(time_to_find)
        
        file_list = self.file_list_hdf()
        time = self.get_PNS_radius(ret_time=True, tob_corrected=tob_corrected)
        if return_index:
            index = np.argmax(time>=time_to_find)
            return file_list[index], index
        return file_list[np.argmax(time>=time_to_find)]

    ## TIME POINTS           
    def time_of_bounce(self):
        """
        Empirical criterion: time of bounce defined as the time at which
        the central density (max desity) raises above 2.5e14 g/cm^3 before a rapid fall.
        If we do not reach that density we lower the threshold to 2
        """
        rho_data = self.rho_max(False)
        rho_index = np.argmax(rho_data[:,1] > 1.4e14)
        if rho_index == 0 or rho_data[rho_index, 0] >= 0.6:
            rho_index = np.argmax(rho_data[:,1]>2e14)
        return rho_data[rho_index, 0]

    def time_of_BH(self):
        """
        Attempt to find when and if a BH is formed.
        For now is when the rho max raises above a certain limit (6e15 g/cm3)
        """
        rho_data = self.rho_max()
        rho_BH = 6e15
        return rho_data[np.argmax(rho_data[:,1]>=rho_BH),0]

    ## -----------------------------------------------------------------
    ## HYDRODYNAMICAL DATA
    ## -----------------------------------------------------------------
    
    @hdf_isopen
    def rho(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['hydro/data'])[..., self.hydroTHD_index['hydro']
                                          ['I_RH']]), self.dim)
    
    ## VELOCITY
    @hdf_isopen
    def radial_velocity(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['thd/data'])[..., self.hydroTHD_index['thd']
                                        ['I_VELX']]), self.dim)
    
    @hdf_isopen
    def theta_velocity(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['thd/data'])[..., self.hydroTHD_index['thd']
                                        ['I_VELY']]), self.dim)
    
    @hdf_isopen
    def phi_velocity(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['thd/data'])[..., self.hydroTHD_index['thd']
                                        ['I_VELZ']]), self.dim)
    
    @hdf_isopen
    def soundspeed(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['thd/data'])[..., self.hydroTHD_index['thd']
                                        ['I_CSND']]), self.dim)
    ## -----------------------------------------------------------------
    ## THERMODYNAMICAL DATA
    ## -----------------------------------------------------------------

    ## THERMODYNAMICAL
    @hdf_isopen
    def gas_pressure(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['thd/data'])[..., self.hydroTHD_index['thd']
                                        ['I_PGAS']]), self.dim)
    
    @hdf_isopen
    def temperature(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['thd/data'])[..., self.hydroTHD_index['thd']
                                        ['I_TMPR']]), self.dim)
    
    @hdf_isopen
    def enthalpy(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['thd/data'])[..., self.hydroTHD_index['thd']
                                        ['I_ENTH']]), self.dim)
    
    @hdf_isopen
    def entropy(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['thd/data'])[..., self.hydroTHD_index['thd']
                                        ['I_ENTR']]), self.dim)
    
    @hdf_isopen
    def adiabatic_index(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['thd/data'])[..., self.hydroTHD_index['thd']
                                        ['I_GAMM']]), self.dim)
    
    ## RELATIVITY AND GRAVITY
    @hdf_isopen
    def lorentz(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['thd/data'])[..., self.hydroTHD_index['thd']
                                        ['I_LRTZ']]), self.dim)
    
    @hdf_isopen
    def gravitational_potential(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['gravpot/data'])), self.dim)
    
    def gravitational_energy(self, file_name):
        return 0.5 * self.rho(file_name) * \
            self.gravitational_potential(file_name)
    
    ## ENERGY
    @hdf_isopen
    def internal_energy(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['thd/data'])[..., self.hydroTHD_index['thd']
                                        ['I_EINT']]), self.dim)
    
    @hdf_isopen
    def nu_heat(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['thd/data'])[..., self.hydroTHD_index['hydro']
                                        ['I_HEAT']]), self.dim)
    
    ## COMPOSITION
    @hdf_isopen
    def Ye(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['hydro/data'])[..., self.hydroTHD_index['hydro']
            ['I_YE']]), self.dim) / self.rho(file_name)
    
    @hdf_isopen
    def neutron_fraction(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['thd/data'])[..., self.hydroTHD_index['thd']
                                        ['I_COMP'][0]]), self.dim)
    
    @hdf_isopen
    def proton_fraction(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['thd/data'])[..., self.hydroTHD_index['thd']
                                        ['I_COMP'][1]]), self.dim)
    
    @hdf_isopen
    def alpha_fraction(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['thd/data'])[..., self.hydroTHD_index['thd']
                                        ['I_COMP'][2]]), self.dim)
    
    @hdf_isopen
    def heavy_fraction(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['thd/data'])[..., self.hydroTHD_index['thd']
                                        ['I_COMP'][3]]), self.dim)
    
    @hdf_isopen
    def Abar(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['thd/data'])[..., self.hydroTHD_index['thd']
                                        ['I_COMP'][4]]), self.dim)
    
    @hdf_isopen
    def Zbar(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['thd/data'])[..., self.hydroTHD_index['thd']
                                        ['I_COMP'][5]]), self.dim)
    
    ## CHEMICAL POTENTIAL
    @hdf_isopen
    def electron_chemical_potential(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['thd/data'])[..., self.hydroTHD_index['thd']
                                        ['I_CPOT'][0]]), self.dim)
    
    @hdf_isopen
    def neutron_chemical_potential(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['thd/data'])[..., self.hydroTHD_index['thd']
                                        ['I_CPOT'][1]]), self.dim)
    
    @hdf_isopen
    def proton_chemical_potential(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['thd/data'])[..., self.hydroTHD_index['thd']
                                        ['I_CPOT'][2]]), self.dim)
    
    @hdf_isopen
    def neutrino_chemical_potential(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['thd/data'])[..., self.hydroTHD_index['thd']
                                        ['I_CPOT'][3]]), self.dim)
    
    ## -----------------------------------------------------------------
    ## PARAMETERS
    ## -----------------------------------------------------------------

    ## ERROR
    @hdf_isopen
    def error(self, file_name):
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['hydro/data'])[..., self.hydroTHD_index['hydro']
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
    
    ## -----------------------------------------------------------------
    ## MAGNETIC FIELDS DATA
    ## -----------------------------------------------------------------

    @hdf_isopen
    def __CT_magnetic_fields(self, file_name):
        """
        Magnetic field at the cells border, use this ONLY to calculate
        streamlines. If you want to plot the actual magnetic fields use
        the 'magnetic_field' method.
        """
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['mag_CT/data'])), self.dim)
    
    @hdf_isopen
    def magnetic_fields(self, file_name):
        """
        Magnetic field at the cells center.
        """
        return self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['mag_vol/data'])), self.dim)
    
    def poloidal_magnetic_fields(self, file_name):
        data = self.magnetic_fields(file_name)
        return np.sqrt(data[..., 0] ** 2 + data[..., 1] ** 2)
    
    def toroidal_magnetic_fields(self, file_name):
        data = self.magnetic_fields(file_name)
        return data[..., 2]
    
    def magnetic_energy(self, file_name):
        """
        Magnetic energy density. Total, poloidal and toroidal.
        """
        data = self.magnetic_fields(file_name)
        return  0.5 * (data[..., 0] ** 2 + data[..., 1] ** 2 \
                       + data[..., 2] ** 2), \
                 0.5 * (data[..., 0] ** 2 + data[..., 1] ** 2), \
                 0.5 * data[..., 2] ** 2
    
    def stream_function(self, file_name, plane):
        return strfct2D(self.__CT_magnetic_fields(file_name), self.cell, 
                        self.ghost, plane)
    
    ## -----------------------------------------------------------------
    ## NEUTRINO DATA
    ## -----------------------------------------------------------------
    
    ## ENERGY DEPENDENT
    @hdf_isopen
    def neutrino_energy_density(self, file_name):
        nu_ene = self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['neutrino/e'])[..., 0]), self.dim)
        nu_ene[..., 2] /= 4
        return nu_ene
    
    @hdf_isopen
    def neutrino_momenta(self, file_name):
        """
        In the comoving rest frame of the fluid are equal to the
        neutrino energy fluxes
        """
        nu_flux = self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['neutrino/e'])[..., 1:]), self.dim)
        ## Insert a new axis to be consistent with the other dimensions
        if self.dim == 1:
            nu_flux = nu_flux[..., None]
        nu_flux[..., 2, :] /= 4
        return nu_flux
    
    @hdf_isopen
    def nutrino_momenta_opacities(self, file_name):
        nu_opac = self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['neutrino/oe'])[..., 1:]), self.dim)
        if self.dim == 1:
            nu_opac = nu_opac[..., None]
        return nu_opac
    
    def neutrino_number_density(self, file_name):
        return self.neutrino_energy_density(file_name) / \
            u.convert_to_erg(self.cell.E_nu()[:, None])
    
    def neutrino_mean_energy(self, file_name):
        """
        Average neutrino energy per cell so 
        <e> = sum_w E_nu(w) / sum_w N_nu(w),
          with w the center of the neutrino bin
        """
        return u.convert_to_MeV(
            self.neutrino_energy_density(file_name).sum(axis=-2) 
                        / self.neutrino_number_density(file_name).sum(axis=-2))
    ## GREY
    @hdf_isopen
    def neutrino_energy_density_grey(self, file_name):
        nu_ene = self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['/neutrinogrey/egrey'])[..., 0]), self.dim)
        nu_ene[..., 2] /= 4
        return nu_ene
    
    @hdf_isopen
    def neutrino_momenta_grey(self, file_name):
        """
        In the comoving rest frame of the fluid are equal to the 
        neutrino energy fluxes
        """
        nu_flux =  self.ghost.remove_ghost_cells(np.squeeze(np.array(
            self.__data_h5['/neutrinogrey/egrey'])[..., 1:]), self.dim)
        if self.dim == 1:
            nu_flux = nu_flux[..., None]
        nu_flux[..., 2, :] /= 4
        return nu_flux

    ## -----------------------------------------------------------------
    ## GRAVIATIONAL WAVES DATA
    ## -----------------------------------------------------------------

    def GW_Amplitudes(self, distance=1, tob_corrected=True, 
                      zero_correction=True, lower_refinement=True):
        """
        Params:
            distance: distance of the observer from the GW source
            zero_correction: shifts up or down the amplitude to make it
            centred with zero
            lower_refinement: on fine timestep refinements GWs are
            usually noisy, so we take 
            one every n points. n defined as the number to reach 0.1 ms
        Returns GWs amplitudes:
        1D:
            No GWs in spherical symmetry
        2D:
            Column 1: time
            Column 2: + ploarization
        3D:
            Column 1: time
            Column 2: + polarization equatorial plane
            Column 3: + polarization polar plane
            Column 4: x polarization equatorial plane
            Column 5: x polarization polar plane
        """
        data = load_file(self.__log_path, self.__grw_path)
        
        if lower_refinement:
            dt = data[1, 2] - data[0, 2]
            new_dt = dt
            n=1
            while new_dt < 5e-5:
                new_dt += dt
                n += 1
            data = data[::n, :]
        
        column_change = find_column_changing_line(self.__log_path,
                                                  self.__grw_path)
        if zero_correction:
            index = np.argmax((data[:, 2] - self.tob)  
            >= -0.01)
        else:
            index = None
        if tob_corrected and zero_correction:
            data[:,2] -= self.tob
        
        return GW_strain(self.dim, column_change, data, index) / distance
    
    def GW_spectrogram(self, distance=1, window_size=10, tob_corrected=True):
        """
        Parameters:
            window_size: value of the time window to use in ms
            distance: distance of the observer from the GW source
            tob_corrected: if the returned timeseries has to be 
            corrected for the tob
        Returns:
            time: timeseries in s
            frequency: aray of the frequencies in Hz
            Zxx: magnitude
        In 3D simulations:
            h_pl_e, h_pl_p, h_cr_e, h_cr_p
        """
        GW_strain = self.GW_Amplitudes(correct_for_tob=False)
        GW_strain[:, 1:] /= distance
        window = 0
        while (np.abs(GW_strain[window,0] - GW_strain[0,0]) < 
               u.convert_to_s(window_size)):
            window += 1
        
        time, frequency, Zxx = GW_spectrogram(self.dim, GW_strain, window)
        if tob_corrected:
            time -= self.tob
        return time, frequency, Zxx
    
    def Deltah(self, peak:Literal['bounce', 'highest']='bounce',
               interval=[None, None], min_time=1.75, max_time=2, distance=1, 
               coordinates=False, tob_corrected=True):
        """
        Returns the Delta h of the gravitational wave strain as defined
        in Richers et al. 2017 (https://arxiv.org/pdf/1701.02752.pdf).
        Basically the difference between the first maximum and minimum 
        postbounce, the first peak that appears is not considered.
        In case the highest peak is selected the amplitude returned is 
        the maxima between left and right.
        min and max time are windows in which to search the low and high
        peaks in the strain
        Return:
            amplitude in cm
            if coordinates
                time, h of the highest and lowest peak
        """
        GWs = self.GW_Amplitudes(distance, tob_corrected)
        indices = GWs_peak_indices(GWs, peak, interval, min_time, max_time)
        Deltah = np.abs(GWs[indices[1], 1] - GWs[indices[2], 1])
        if coordinates:
            x = [GWs[indices[1], 0], GWs[indices[2], 0]]
            y = [GWs[indices[1], 1], GWs[indices[2], 1]]
            return Deltah, np.array(x), np.array(y)
        return Deltah
    
    def GWs_peak_frequencies(self, peak:Literal['bounce', 'highest']='bounce',
                              min_time=1.75, max_time=2, interval=[None, None],
                              return_intensities=False, return_fourier=False):
        """
        Calculates the dominant and the second dominant frequency of a GWs peak
        Return:
            frequencies: dominat, second dominant
            if return intensities
                intensities: dominant, second dominant
            if return return fourier
                frequencies array
                htilde
        """
        GWs = self.GW_Amplitudes()
        indices = GWs_peak_indices(GWs, peak, interval, min_time, max_time)
        frequency, htilde = GWs_fourier_transform(GWs, indices)
        indices = GWs_frequency_peak_indices(frequency, htilde)
        return_list = [frequency[indices]]
        if return_intensities:
            return_list.append(htilde[indices])
        if return_fourier:
            return_list.append([frequency, htilde])
        if len(return_list) == 1:
            return return_list[0]
        return return_list
    
    def GWs_dE_dt(self, tob_corrected=True):
        """
        Returns the energy carried away by the GWs in erg/s
        """
        GWs = self.GW_Amplitudes(tob_corrected)
        GWs[:, 1] = u.speed_light ** 3 / u.G * 2 / 15 * GWs[:,1] ** 2
        return GWs


    ## FIX ME!!!!!
    def AE220(self, tob_corrected):
        """
        Calculates the AE220 from density and velocities for a
        2D simulation.
        ONLY 2D
        Returns
            radius
            time: array of time step
            AE220: len(radius), len(time) array
            full_strain: GWs strain from the full star
            inner_strain: GWs strain from the PNS nucleus
            convection_strain: GWs strain from the convection region
            outer_strain: GWs strain from the outer layers
        """
        AE220_file = os.path.join(self.storage_path,'AE220.h5')
        const = -0.125 *  np.sqrt(15/np.pi)
        if not os.path.exists(AE220_file):
            warnings.warn("AE220 file not found. Creating one. \nPlease wait...")
            time, AE220, full_strain, inner_strain, \
            convection_strain, outer_strain = calculate_AE220(self)
            AE220_hdf = h5py.File(AE220_file, 'w')
            AE220_hdf.create_dataset('time', data = time)
            AE220_hdf.create_dataset('AE220', data = AE220)
            AE220_hdf.create_dataset('h_full_strain', data = full_strain * const)
            AE220_hdf.create_dataset('h_inner_strain', data = inner_strain * const)
            AE220_hdf.create_dataset('h_convection_strain', data = convection_strain * const)
            AE220_hdf.create_dataset('h_outer_strain', data = outer_strain * const)
            AE220_hdf.close()
        data = self.open_h5(AE220_file)
        time = data["time"][...]
        AE220 = data["AE220"][...]
        full_strain = data["h_full_strain"][...]
        inner_strain = data["h_inner_strain"][...]
        convection_strain = data["h_convection_strain"][...]
        outer_strain = data["h_outer_strain"][...]
        self.close_h5(data)
        
        if not tob_corrected:
            time += self.time_of_bounce_rho()
        return self.cell.radius(self.ghost), time, AE220 * const, \
                full_strain, inner_strain, convection_strain, outer_strain
    
    ## -----------------------------------------------------------------
    ## GLOBAL DATA
    ## -----------------------------------------------------------------
    def global_neutrino_luminosity(self, tob_corrected=True):
        """
        indices
        1: time
        2: luminosity flux nue  5: number luminosity flux nue
        3: luminosity flux nua  6: number luminosity flux nua
        4: luminosity flux nux  7: number luminosity flux nux
        """
        nu_tmp = load_file(self.__log_path, self.__integrated_nu_path)
        if tob_corrected:
            nu_tmp[:, 2] -= self.tob
        return np.stack((nu_tmp[:, 2], nu_tmp[:, 38], nu_tmp[:, 39],
                         0.25 * nu_tmp[:, 40], nu_tmp[:, 35],
                         nu_tmp[:, 36], 0.25 * nu_tmp[:, 37]), axis=1)
    
    def global_neutrino_mean_energies(self, tob_corrected=True):
        """
        indices
        1: time
        2: mean energy nue  
        3: mean energy nua  
        4: mean energy nux  
        """
        nu = self.global_neutrino_luminosity(tob_corrected)
        ene_mean = np.zeros((nu.shape[0],4))
        ene_mean[:,0] = nu[:,0]
        ene_mean[:,1] = u.convert_to_MeV(nu[:,1]/nu[:,4])
        ene_mean[:,2] = u.convert_to_MeV(nu[:,2]/nu[:,5])
        ene_mean[:,3] = u.convert_to_MeV(nu[:,3]/nu[:,6])
        return np.stack((nu[:, 0], u.convert_to_MeV(nu[:, 1]/nu[:, 4]),
                         u.convert_to_MeV(nu[:, 2]/nu[:, 5]),
                         u.convert_to_MeV(nu[:, 3]/nu[:, 6])), axis=1)
    
    def total_mass(self, tob_corrected=True):
        """
        Returns the total mass of the star at every timestep
        indices
        1: time
        2: total mass
        """
        M = load_file(self.__log_path, self.__rho_max_path)
        if tob_corrected:
            M[:,2] -= self.tob
        return np.stack((M[:, 2], u.convert_to_solar_masses(M[:, 4])), axis=1)
    
    def rho_max(self, correct_for_tob=True):
        """
        indices
        1: time
        2: rho max
        """
        rho = load_file(self.__log_path, self.__rho_max_path)
        if correct_for_tob:
            rho[:,2] -= self.tob
        return np.stack((rho[:, 2], rho[:, 3]), axis = 1)

    def global_Ye(self, tob_corrected=True):
        """
        indices
        1: time
        2: Ye max
        3: Ye min
        4: Ye cent
        """
        Ye = load_file(self.__log_path, self.__rho_max_path)
        if tob_corrected:
            Ye[:,2] -= self.tob
        return np.stack((Ye[:,2], Ye[:,6], Ye[:,7], Ye[:,8]), axis=1)

    def global_rotational_energy(self, tob_corrected=True):
        """
        1: time
        2: total rotational energy
        """
        en = load_file(self.__log_path, self.__mag_data)
        if tob_corrected:
            en[:,2] -= self.tob
        return np.stack((en[:, 2], en[:, 3]), axis=1)

    ## -----------------------------------------------------------------
    ## PNS DATA
    ## -----------------------------------------------------------------