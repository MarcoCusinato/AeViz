from AeViz.simulation.methods import *
from AeViz.utils.file_utils import load_file
from AeViz.units.aeseries import aeseries
from AeViz.units.aerray import aerray
from AeViz.units import u

"""
Function to load and process log files from a simulation.
These functions are not meant to be used standalone, but rather to be
imported into the Simulation class.
"""

## -----------------------------------------------------------------
## GLOBAL DATA
## -----------------------------------------------------------------

@smooth
@derive
@subtract_tob
def global_neutrino_luminosity(self, tob_corrected=True, **kwargs):
    """
    aeseries
    luminosity flux nue  5: number luminosity flux nue
    luminosity flux nua  6: number luminosity flux nua
    luminosity flux nux  7: number luminosity flux nux
    """
    nu_tmp = load_file(self._Simulation__log_path,
                       self._Simulation__integrated_nu_path)
    time = aerray(nu_tmp[:,2], u.s, 'time', r'$t$', None, [0, nu_tmp[-1, 2]], False)
    return [aeseries(
                    aerray(nu_tmp[:, 38], u.erg / u.s, 'Lnue',
                           r'$L_\mathrm{\nu_e}$', None, [0, 1e52]),
                           time=time),
            aeseries(
                    aerray(nu_tmp[:, 39], u.erg / u.s, 'Lnua',
                           r'$L_\mathrm{\overline{\nu}_e}$', None, [0, 1e52]),
                           time=time.copy()),
            aeseries(
                    aerray(nu_tmp[:, 40] / 4, u.erg / u.s, 'Lnux',
                           r'$L_\mathrm{\nu_x}$', None, [0, 1e52]),
                           time=time.copy()),
            ]

@subtract_tob
def global_neutrino_number_luminosity(self, tob_corrected=True, **kwargs):
    """
    aeseries
    number luminosity flux nue
    number luminosity flux nua
    number luminosity flux nux
    """
    nu_tmp = load_file(self._Simulation__log_path,
                       self._Simulation__integrated_nu_path)
    time = aerray(nu_tmp[:,2], u.s, 'time', r'$t$', None, [0, nu_tmp[-1, 2]], False)
    return [aeseries(
                    aerray(nu_tmp[:, 35], u.dimensionless_unscaled / u.s,
                           'Lnue',
                           r'$N_\mathrm{\nu_e}$', None, [0, 1e55]),
                           time=time),
            aeseries(
                    aerray(nu_tmp[:, 36], u.dimensionless_unscaled / u.s,
                           'Lnumnua',
                           r'$N_\mathrm{\overline{\nu}_e}$', None, [0, 1e55]),
                           time=time.copy()),
            aeseries(
                    aerray(nu_tmp[:, 37] / 4, u.dimensionless_unscaled / u.s,
                           'Lnumnux',
                           r'$N_\mathrm{\nu_x}$', None, [0, 1e55]),
                           time=time.copy()),
    ]

@smooth
@derive
def global_neutrino_mean_energies(self, tob_corrected=True, **kwargs):
    """
    list of aeseries
    1: mean energy nue  
    2: mean energy nua  
    3: mean energy nux  
    """
    nu_lum = self.global_neutrino_luminosity(tob_corrected)
    nu_num = self.global_neutrino_number_luminosity(tob_corrected)
    mean_ene = [(nu_lum[i] / nu_num[i]).to(u.MeV) for i in range(len(nu_lum))]
    for me, nm, lb in zip(mean_ene, ['Enue', 'Enua', 'Enux'],
                          [r'$\langle E_{\nu_e}\rangle$',
                           r'$\langle E_{\overline{\nu}_e}\rangle$',
                           r'$\langle E_{\nu_x}\rangle$']):
        me.data.set(name=nm, limits=[0,30], label=lb) 
    return mean_ene

@smooth
@derive
@subtract_tob
def total_mass(self, tob_corrected=True, **kwargs):
    """
    Aeseries with mass in solar masses
    """
    file = load_file(self._Simulation__log_path, self._Simulation__rho_max_path)
    M = aerray(file[:, 4], u.g, 'mtot', r'$M_\mathrm{tot}$', limits=[0, 10]).to(u.M_sun)
    time = aerray(file[:, 2], u.s, 'time', r'$t$', None, [0, file[-1, 2]], False)
    return aeseries(M, time=time)

@smooth
@derive
@subtract_tob
def rho_max(self, tob_corrected=True, **kwargs):
    """
    indices
    1: time
    2: rho max
    """
    rho = load_file(self._Simulation__log_path, self._Simulation__rho_max_path)
    time = aerray(rho[:,2], u.s, 'time', r'$t$', None, [0, rho[-1, 2]], False)       
    return aeseries(
        aerray(rho[:, 3], u.g / u.cm ** 3, 'rho_max', r'$\rho_\mathrm{max}$',
               'viridis', [1e4, 1e15], True),
        time=time
    )

@smooth
@derive
@subtract_tob
def Ye_max(self, tob_corrected=True, **kwargs):
    Ye = load_file(self._Simulation__log_path, self._Simulation__rho_max_path)
    return aeseries(
        aerray(Ye[:, 6], u.dimensionless_unscaled, 'Yemax', r'$Y_\mathrm{e,max}$',
               None, [0, 0.5], False),
        aerray(Ye[:,2], u.s, 'time', r'$t$', None, [0, Ye[-1, 2]], False)
    )

@smooth
@derive
@subtract_tob
def Ye_min(self, tob_corrected=True, **kwargs):
    Ye = load_file(self._Simulation__log_path, self._Simulation__rho_max_path)
    return aeseries(
        aerray(Ye[:, 7], u.dimensionless_unscaled, 'Yemin', r'$Y_\mathrm{e,min}$',
               None, [0, 0.5], False),
        aerray(Ye[:,2], u.s, 'time', r'$t$', None, [0, Ye[-1, 2]], False)
    )

@smooth
@derive
@subtract_tob
def Ye_min(self, tob_corrected=True, **kwargs):
    Ye = load_file(self._Simulation__log_path, self._Simulation__rho_max_path)
    return aeseries(
        aerray(Ye[:, 8], u.dimensionless_unscaled, 'Yecent', r'$Y_\mathrm{e,cent}$',
               None, [0, 0.5], False),
        aerray(Ye[:,2], u.s, 'time', r'$t$', None, [0, Ye[-1, 2]], False)
    )

@smooth
@derive
@subtract_tob
def rotational_energy_total(self, tob_corrected=True, **kwargs):
    """
    1: time
    2: total rotational energy
    """
    en = load_file(self._Simulation__log_path, self._Simulation__mag_data)
    return aeseries(
        aerray(en[:, 3], u.dimensionless_unscaled, 'Erottot', r'$E_\mathrm{rot,tot}$',
               None, [0, 1e53], False),
        aerray(en[:, 2], u.s, 'time', r'$t$', None, [0, en[-1, 2]], False)
    )