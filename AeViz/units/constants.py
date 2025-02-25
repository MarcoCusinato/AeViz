from AeViz.units.aerray import aerray
from astropy import units as u
import numpy as np

class constants:
    c = aerray(2.99792458e10, u.cm / u.s, 'speed_of_light', label=r'$c$')
    G = aerray(6.6743e-8, u.cm ** 3 / u.g / u.s ** 2, 'gravitational_constant',
               label=r'$G$')
    h = aerray(6.62607015e-27, u.erg * u.s, 'planck_constant', label='$h$')
    hbar = aerray(1.0545718e-27, u.erg * u.s, 'reduced_planck_constant',
                  label=r'$\hbar$')
    kB = aerray(1.3807e-16, u.cm ** 2 / u.g / u.s ** 3 / u.K,
                'boltzmann_constant', label=r'$k_\mathrm{B}$')
    sSB = aerray(5.6704e-5, u.g / u.s ** 3 / u.K ** 4,
                 name='stefan_boltzmann_constant',
                 label=r'$\sigma_\mathrm{SB}$')
    lambda_cosm = aerray(1.089e-56, u.cm ** (-2), name='cosmological_constant',
                         label=r'$\Lambda$')
    e = aerray(4.8032e-10, u.cm ** (3/2) * u.g ** (1/2) / u.s, name='e',
               label='$e$')
    me = aerray(9.1093837139e-28, u.g, name='electron_mass',
                label=r'$m_\mathrm{e}$')
    mp = aerray(1.67262192595e-24, u.g, name='proton_mass',
                label=r'$m_\mathrm{p}$')
    mn = aerray(1.67492750056-24, u.g, name='neutron_mass',
                label=r'$m_\mathrm{n}$')
    Msol = aerray(1.989e33, u.g, name='solar_mass', label=r'$M_\odot$')
    Rsol = aerray(6.955e10, u.cm, name='solar_radius', label=r'$R_\odot$')
    mu0 = aerray(1, u.cm * u.s **2 * u.G ** 2 / u.g, name='mu0', 
                 label=r'$\mu_0$')