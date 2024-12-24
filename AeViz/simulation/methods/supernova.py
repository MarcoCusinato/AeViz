from AeViz.simulation.methods import *
from AeViz.utils.kick_vel_utils import calculate_kick
from AeViz.utils.load_save_mass_ene_utils import calculate_masses_energies
from AeViz.utils.load_save_radii_utils import calculate_radius
from AeViz.utils.PNS_ang_mom_nu_utils import calculate_angular_mom_PNS_nu
import os

"""
Methods to calculate and handle data from a supernova simulation.
So far any simulation in spherical coordinates is a supernova
simulation.
These functions are not meant to be used standalone, but rather to be
imported into the Simulation class.
"""

## TIME POINTS           
def time_of_bounce(self):
    """
    Empirical criterion: time of bounce defined as the time at which
    the central density (max desity) raises above 2.5e14 g/cm^3 
    before a rapid fall. If we do not reach that density we lower 
    the threshold to 2
    """
    if os.path.exists(os.path.join(self.storage_path, 'tob.dat')):
        tob = np.loadtxt(os.path.join(self.storage_path, 'tob.dat'))
        return tob
    rho_data = self.rho_max(False)
    rho_index = np.argmax(rho_data[:,1] > 1.4e14)
    if rho_index == 0 or rho_data[rho_index, 0] >= 0.6:
        rho_index = np.argmax(rho_data[:,1]>2e14)
    return rho_data[rho_index, 0]

def time_of_explosion(self):
    """
    Empirical criterion: time of explosion defined as the time at
    which the explosion energy raises above 5e48 erg
    """
    time, _, ene, *_ = self.explosion_mass_ene()
    _, _, shock_max, *_ = self.shock_radius()
    index = np.where((ene > 1e48) & (shock_max > 3e7))[0][0]
    return time[index]

def time_of_BH(self):
    """
    Attempt to find when and if a BH is formed.
    For now is when the rho max raises above a certain limit
    (6e15 g/cm3)
    """
    rho_data = self.rho_max()
    rho_BH = 6e15
    return rho_data[np.argmax(rho_data[:,1]>=rho_BH),0]

## BOUNCE COMPACTNESS
def bounce_compactness(self, mass = None):
    """
    Compactness parameter as defined in O'Connor 2011  
    `10.1088/0004-637X/730/2/72`
    It returns the compacteness at bounce (maximum compactness of
    the core) for the model. 
    If a mass at which the compactness has to be provided is 
    specified it will return just a value
    """
    bounce_file = self.find_file_from_time(0, True)
    rho = self.rho(bounce_file)
    radius = u.convert_to_km(self.cell.radius(self.ghost)) / 1000
    dV = self.cell.dVolume_integration(self.ghost)
    if self.dim > 1:
        enclosed_mass = np.sum(rho * dV, axis = tuple(range(self.dim - 1)))
    enclosed_mass = u.convert_to_solar_masses(np.cumsum(enclosed_mass))
    compactness = enclosed_mass / radius
    if mass is not None:
        return compactness[ np.argmax( enclosed_mass >= mass ) ]
    else:
        return compactness

## -----------------------------------------------------------------
## RADII DATA
## -----------------------------------------------------------------

@smooth
@derive
def PNS_radius(self, tob_corrected=True, save_checkpoints=True, **kwargs):
    """
    Returns the radius of the PNS at every timestep.
    If tob_corrected is True, the time is corrected for the time of
    bounce. If save_checkpoints is True, the checkpoints are saved
    during the calculation.
    Returns: time, radius(phi, theta), max_radius, min_radius,
                average_radius, number of ghost cells
    """
    data = calculate_radius(self, 'PNS', save_checkpoints)
    if not tob_corrected:
        data[0] += self.tob
    return data

@smooth
@derive
def shock_radius(self, tob_corrected=True, save_checkpoints=True, **kwargs):
    """
    Returns the shock radius at every timestep.
    If tob_corrected is True, the time is corrected for the time of
    bounce. If save_checkpoints is True, the checkpoints are saved
    during the calculation.
    Returns: time, radius(phi, theta), max_radius, min_radius,
                average_radius, number of ghost cells
    """
    data = calculate_radius(self, 'shock', save_checkpoints)
    if not tob_corrected:
        data[0] += self.tob
    return data

@smooth
@derive
def neutrino_spheres(self, tob_corrected=True, save_checkpoints=True, **kwargs):
    """
    Returns the neutrino spheres at every timestep.
    If tob_corrected is True, the time is corrected for the time of
    bounce. If save_checkpoints is True, the checkpoints are saved
    during the calculation.
    Returns: time, radius(phi, theta), max_radius, min_radius,
                average_radius, number of ghost cells
    Radii are returned as disctionaries of nue, nua and nux
    """
    data = calculate_radius(self, 'neutrino', save_checkpoints)
    if not tob_corrected:
        data[0] += self.tob
    return data

@smooth
@derive
def gain_radius(self, tob_corrected=True, save_checkpoints=True, **kwargs):
    """
    Returns the gain radius at every timestep.
    If tob_corrected is True, the time is corrected for the time of
    bounce. If save_checkpoints is True, the checkpoints are saved
    during the calculation.
    Returns: time, radius(phi, theta), max_radius, min_radius,
                average_radius, number of ghost cells
    """
    data = calculate_radius(self, 'gain', save_checkpoints)
    if not tob_corrected:
        data[0] += self.tob
    return data

@smooth
@derive
def PNS_nucleus_radius(self, tob_corrected=True, save_checkpoints=True, **kwargs):
    """
    Returns the PNS nucleus at every timestep.
    If tob_corrected is True, the time is corrected for the time of
    bounce. If save_checkpoints is True, the checkpoints are saved
    during the calculation.
    Returns: time, radius(phi, theta), max_radius, min_radius,
                average_radius, number of ghost cells
    """
    data = calculate_radius(self, 'nucleus', save_checkpoints)
    if not tob_corrected:
        data[0] += self.tob
    return data

@smooth
@derive
def innercore_radius(self, tob_corrected=True, save_checkpoints=True, **kwargs):
    """
    Returns the inner core radius at every timestep.
    If tob_corrected is True, the time is corrected for the time of
    bounce. If save_checkpoints is True, the checkpoints are saved
    during the calculation.
    Returns: time, radius(phi, theta), max_radius, min_radius,
                average_radius, number of ghost cells
    """
    data = calculate_radius(self, 'innercore', save_checkpoints)
    if not tob_corrected:
        data[0] += self.tob
    return data

## -----------------------------------------------------------------
## MASS AND ENERGY DATA
## -----------------------------------------------------------------

@smooth
@derive
def PNS_mass_ene(self, tob_corrected=True, save_checkpoints=True, **kwargs):
    """
    Returns the PNS mass and energy at every timestep.
    If tob_corrected is True, the time is corrected for the time of
    bounce. If save_checkpoints is True, the checkpoints are saved
    during the calculation.
    Returns: time, mass, kinetic energy, magnetic energy,
                rotational energy, gravitational energy, total energy,
                convective energy, T/W
    """
    time, _, _, _, data, *_= \
        calculate_masses_energies(self, save_checkpoints)
    if not tob_corrected:
        time += self.tob
    return [time, data['mass'], data['kinetic_ene'], data['magnetic_ene'],
            data['rotational_ene'], data['grav_ene'], data['total_ene'],
            data['convective_ene'], data['rotational_ene'] / \
                np.abs(data['grav_ene'])]

@smooth
@derive
def PNS_angular_mom(self, tob_corrected=True, save_checkpoints=True,
                    **kwargs):
    """
    Returns the PNS angular momentum at every timestep.
    If tob_corrected is True, the time is corrected for the time of
    bounce. If save_checkpoints is True, the checkpoints are saved
    during the calculation.
    Returns: time, Lx, Ly, Lz, L_tot
    """
    time, _, _, _, data, *_ = \
        calculate_masses_energies(self, save_checkpoints)
    if not tob_corrected:
        time += self.tob
    return [time, data['L']['Lx'], data['L']['Ly'],
            data['L']['Lz'], data['L']['L_tot']]

@smooth
@derive
def PNS_angular_momentum_neutrinos(self, tob_corrected=True,
                                    save_checkpoints=True, 
                                    **kwargs):
    """
    Returns the PNS angular momentum carried away by neutrinos at every
    timestep.
    If tob_corrected is True, the time is corrected for the time of
    bounce. If save_checkpoints is True, the checkpoints are saved
    during the calculation.
    Returns: time, Lx, Ly, Lz, Lx_tot, Ly_tot, Lz_tot, L_tot
        The momentum Lx, Ly and Lz are dictionaries
    """
    time, Lx, Ly, Lz, Lx_tot, Ly_tot, Lz_tot, L_tot = \
        calculate_angular_mom_PNS_nu(self, save_checkpoints)
    if not tob_corrected:
        time += self.tob
    return [time, Lx, Ly, Lz, Lx_tot, Ly_tot, Lz_tot, L_tot]

@smooth
@derive
def explosion_mass_ene(self, tob_corrected=True, save_checkpoints=True, **kwargs):
    """
    Returns the explosion mass and energy at every timestep.
    If tob_corrected is True, the time is corrected for the time of
    bounce. If save_checkpoints is True, the checkpoints are saved
    during the calculation.
    Returns: time, unbound mass, energy, and 
        kinetic energy, magnetic energy of unbounded material
    """
    time, _, _, _, _, data, *_ = \
        calculate_masses_energies(self, save_checkpoints)
    if not tob_corrected:
        time += self.tob
    return [time, data['mass'], data['energy'], data['kinetic_ene'],
            data['magnetic_ene']]

@smooth
@derive
def gain_mass_nu_heat(self, tob_corrected=True, save_checkpoints=True, **kwargs):
    """
    Returns the gain mass and neutrino heating at every timestep.
    If tob_corrected is True, the time is corrected for the time of
    bounce. If save_checkpoints is True, the checkpoints are saved
    during the calculation.
    Returns: time, mass, neutrino heating
    """
    time, _, _, data, *_ = \
        calculate_masses_energies(self, save_checkpoints)
    if not tob_corrected:
        time += self.tob
    return [time, data['mass'], data['heating_ene']]

@smooth
@derive
def innercore_mass_ene(self, tob_corrected=True, save_checkpoints=True, **kwargs):
    """
    Returns the inner core mass and energy at every timestep.
    If tob_corrected is True, the time is corrected for the time of
    bounce. If save_checkpoints is True, the checkpoints are saved
    during the calculation.
    Returns: time, mass, kinetic energy, magnetic energy,
                rotational energy, gravitational energy, total energy,
                T/W
    """
    time, _, data, *_ = \
        calculate_masses_energies(self, save_checkpoints)
    if not tob_corrected:
        time += self.tob
    return [time, data['mass'], data['kinetic_ene'], data['magnetic_ene'],
            data['rotational_ene'], data['grav_ene'], data['total_ene'],
            data['T_W']]

@smooth
@derive
def PNS_core_mass_ene(self, tob_corrected=True, save_checkpoints=True, **kwargs):
    """
    Returns the PNS core mass and energy at every timestep.
    If tob_corrected is True, the time is corrected for the time of
    bounce. If save_checkpoints is True, the checkpoints are saved
    during the calculation.
    Returns: time, mass, kinetic energy, magnetic energy,
                rotational energy, gravitational energy, total energy,
                T/W
    """
    time, _, _, _, _, _, data = \
        calculate_masses_energies(self, save_checkpoints)
    if not tob_corrected:
        time += self.tob
    return [time, data['mass'], data['kinetic_ene'], data['magnetic_ene'],
            data['rotational_ene'], data['grav_ene'], data['total_ene'],
            data['T_W']]

@smooth
@derive
def mass_accretion_500km(self, tob_corrected=True, save_checkpoints=True,
                            **kwargs):
    """
    Returns the mass accretion rate at 500 km from the center of the
    star at every timestep.
    If tob_corrected is True, the time is corrected for the time of
    bounce.
    Returns: time, mass accretion rate
    """
    time, data, *_ = \
        calculate_masses_energies(self, save_checkpoints)
    if not tob_corrected:
        time += self.tob
    return [time, data]

## -----------------------------------------------------------------
## VELOCITIES DATA
## -----------------------------------------------------------------

@smooth
@derive
def PNS_kick_velocity(self, tob_corrected=True, save_checkpoints=True,
                        **kwargs):
    """
    Returns the modeule of the PNS kick velocity at every timestep.
    If tob_corrected is True, the time is corrected for the time of
    bounce. If save_checkpoints is True, the checkpoints are saved
    during the calculation.
    Returns: time, kick velocity, hydro kick velocity, 
                nu kick velocity
    """
    def modulus(v):
        vtot = 0
        for comp in v:
            vtot += comp ** 2
        return vtot ** 0.5
    def sum_components(vcomp):
        vtot = 0
        for comp in vcomp:
            vtot += comp
        return vtot
        
    time, hydro, vnue, vnua, vnux = \
                    self.PNS_kick_velocity_components(tob_corrected,
                                                    save_checkpoints)
    vkick = modulus([sum_components([hydro[0], vnue[0], vnua[0], vnux[0]]),
                    sum_components([hydro[1], vnue[1], vnua[1], vnux[1]]),
                    sum_components([hydro[2], vnue[2], vnua[2], vnux[2]])])
    vkick_hydro = modulus(hydro)
    vkick_nue = modulus([sum_components([vnue[0], vnua[0], vnux[0]]),
                        sum_components([vnue[1], vnua[1], vnux[1]]),
                        sum_components([vnue[2], vnua[2], vnux[2]])])
    return time, vkick, vkick_hydro, vkick_nue

@smooth
@derive
def PNS_kick_velocity_components(self, tob_corrected=True,
                                    save_checkpoints=True, **kwargs):
    """
    Returns the components of the PNS kick velocity at every timestep.
    If tob_corrected is True, the time is corrected for the time of
    bounce. If save_checkpoints is True, the checkpoints are saved
    during the calculation.
    Returns: time, vr, vtheta, vphi, v_nu
            v_nu is returned for each neutrino species (nue, nuebar, nux)
    """
    time, hydro, nu_flux = \
        calculate_kick(self, save_checkpoints)
    if not tob_corrected:
        time += self.tob
    dt = np.zeros(time.shape[0])
    dt[1:] = time[1:] - time[:-1]
    dt[0] = dt[1]
    _, PNSmass,*_ = self.PNS_mass_ene()
    vnue = [np.cumsum(comp * dt) / PNSmass for comp in [nu_flux['nue']['x'],
                                                    nu_flux['nue']['y'],
                                                    nu_flux['nue']['z']]]
    vnua = [np.cumsum(comp * dt) / PNSmass for comp in [nu_flux['nua']['x'],
                                                    nu_flux['nua']['y'],
                                                    nu_flux['nua']['z']]]
    vnux = [np.cumsum(comp * dt) / PNSmass for comp in [nu_flux['nux']['x'],
                                                    nu_flux['nux']['y'],
                                                    nu_flux['nux']['z']]]
    hydro = [comp / PNSmass for comp in [hydro['x'], hydro['y'],
                                            hydro['z']]]
    return time, hydro, vnue, vnua, vnux, 
