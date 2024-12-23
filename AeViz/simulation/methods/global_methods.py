from AeViz.simulation.methods import *
from AeViz.utils.file_utils import load_file

## -----------------------------------------------------------------
## GLOBAL DATA
## -----------------------------------------------------------------

@smooth
@derive
def global_neutrino_luminosity(self, tob_corrected=True, **kwargs):
    """
    indices
    1: time
    2: luminosity flux nue  5: number luminosity flux nue
    3: luminosity flux nua  6: number luminosity flux nua
    4: luminosity flux nux  7: number luminosity flux nux
    """
    nu_tmp = load_file(self._Simulation__log_path,
                       self._Simulation__integrated_nu_path)
    if tob_corrected:
        nu_tmp[:, 2] -= self.tob
    return np.stack((nu_tmp[:, 2], nu_tmp[:, 38], nu_tmp[:, 39],
                        0.25 * nu_tmp[:, 40], nu_tmp[:, 35],
                        nu_tmp[:, 36], 0.25 * nu_tmp[:, 37]), axis=1)

@smooth
@derive
def global_neutrino_mean_energies(self, tob_corrected=True, **kwargs):
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

@smooth
@derive
def total_mass(self, tob_corrected=True, **kwargs):
    """
    Returns the total mass of the star at every timestep
    indices
    1: time
    2: total mass
    """
    M = load_file(self._Simulation__log_path, self._Simulation__rho_max_path)
    if tob_corrected:
        M[:,2] -= self.tob
    return np.stack((M[:, 2], u.convert_to_solar_masses(M[:, 4])), axis=1)

@smooth
@derive
def rho_max(self, correct_for_tob=True, **kwargs):
    """
    indices
    1: time
    2: rho max
    """
    rho = load_file(self._Simulation__log_path, self._Simulation__rho_max_path)
    if correct_for_tob:
        rho[:,2] -= self.tob
    return np.stack((rho[:, 2], rho[:, 3]), axis = 1)

@smooth
@derive
def global_Ye(self, tob_corrected=True, **kwargs):
    """
    indices
    1: time
    2: Ye max
    3: Ye min
    4: Ye cent
    """
    Ye = load_file(self._Simulation__log_path, self._Simulation__rho_max_path)
    if tob_corrected:
        Ye[:,2] -= self.tob
    return np.stack((Ye[:,2], Ye[:,6], Ye[:,7], Ye[:,8]), axis=1)

@smooth
@derive
def global_rotational_energy(self, tob_corrected=True, **kwargs):
    """
    1: time
    2: total rotational energy
    """
    en = load_file(self._Simulation__log_path, self._Simulation__mag_data)
    if tob_corrected:
        en[:,2] -= self.tob
    return np.stack((en[:, 2], en[:, 3]), axis=1)