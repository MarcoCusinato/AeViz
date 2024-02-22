import numpy as np
## PLOTTING UTILS DICTIONARIES

## GW limits
def GW_limit(GW_data):
       lim = np.max(np.abs(GW_data))
       lim *= 1.1
       if lim > 150:
              lim = 150
       return (-lim, lim)

## PLOT LABELS, COLORMAPS, AND Y-LIMITS or COLORBAR LIMITS and LIMITS
## SCALES (LOG, LINEAR, or SYMLOG)
plot_labels = {
        'rho': {'log': True,
                'lim': (1e4, 1e15),
                'cmap': 'viridis',
                'label': r'$\rho$ [g$\cdot$cm$^{-3}$]'},
        'rho_max': {'log': True,
                'lim': (1e4, 1e15),
                'cmap': 'viridis',
                'label': r'$\rho$ [g$\cdot$cm$^{-3}$]'},
        'MHD_energy': {'log': True,
                'lim': (1e24, 1e35),
                'cmap': 'nipy_spectral',
                'label': r'$E$ [erg]'},
        'VX': {'log': True,
               'lim': (-1e10, 5e10),
               'cmap': 'Spectral_r',
               'label': r'$v_r$ [cm$\cdot$s$^{-1}$]'},
        'VY': {'log': True,
               'lim': (-1e10, 1e10),
               'cmap': 'Spectral_r',
               'label': r'$v_\theta$ [cm$\cdot$s$^{-1}$]'},
        'VZ': {'log': True,
               'lim': (1e7, 1e10),
               'cmap': 'cividis',
               'label': r'$v_\phi$ [cm$\cdot$s$^{-1}$]'},
        'Ye': {'log': False,
               'lim': (0.0, 0.5),
               'cmap': 'gist_rainbow',
               'label': r'Y$_e$'},
        'YN': {'log': False,
               'lim': (0, 1.0),
               'cmap': 'gist_rainbow',
               'label': r'Y$_N$'},
        'error': {'log': False,
                'lim': (0, 1),
                'cmap': 'seismic',
                'label': r'Error'},
        'lorentz': {'log': False,
                 'lim': (1, 1.1),
                 'cmap': 'gist_rainbow',
                 'label': r'$\gamma$'},
        'DENS': {'log': True,
                 'lim': (1e4, 1e15),
                 'cmap': 'viridis',
                 'label': r'$\rho$ [g cm$^{-3}$]'},
        'internal_energy': {'log': True,
                 'lim': (1e24, 1e35),
                 'cmap': 'nipy_spectral',
                 'label': r'$E_{\mathrm{int}}$ [erg]'},
        'enthalphy': {'log': True,
                 'lim': (1e25, 1e36),
                 'cmap': 'gist_stern',
                 'label': r'$H$ [erg'},
        'PELE': {'log': True,
                 'lim': (1e4, 1e15),
                 'cmap': 'viridis',
                 'label': r'$\rho$ [g cm$^{-3}$]'},
        'TELE': {'log': True,
                 'lim': (1e4, 1e15),
                 'cmap': 'viridis',
                 'label': r'$\rho$ [g cm$^{-3}$]'},
        'NELE': {'log': True,
                 'lim': (1e4, 1e15),
                 'cmap': 'viridis',
                 'label': r'$\rho$ [g cm$^{-3}$]'},
        'PION': {'log': True,
                 'lim': (1e4, 1e15),
                 'cmap': 'viridis',
                 'label': r'$\rho$ [g cm$^{-3}$]'},
        'TION': {'log': True,
                 'lim': (1e4, 1e15),
                 'cmap': 'viridis',
                 'label': r'$\rho$ [g cm$^{-3}$]'},
        'NION': {'log': True,
                 'lim': (1e4, 1e15),
                 'cmap': 'viridis',
                 'label': r'$\rho$ [g cm$^{-3}$]'},
        'radial_velocity': {'log': True,
                 'lim': (-1e10, 5e10),
                 'cmap': 'Spectral_r',
                 'label': r'$v_r$ [cm$\cdot$s$^{-1}$]'},
        'theta_velocity':  {'log': True,
                  'lim': (-1e10, 1e10),
                  'cmap': 'Spectral_r',
                  'label': r'$v_\theta$ [cm$\cdot$s$^{-1}$]'},
        'phi_velocity':  {'log': True,
                  'lim': (1e7, 1e10),
                  'cmap': 'cividis',
                  'label': r'$v_\phi$ [cm$\cdot$s$^{-1}$]'},
        'temperature': {'log': False,
              'lim': (0, 40),
              'cmap': 'inferno',
              'label': r'T [MeV]'},
        'entropy': {'log': False,
                 'lim': (1.5, 15),
                 'cmap': 'gist_rainbow_r',
                 'label': r'S [k$_b$/bry]'},
        'adiabatic_index': {'log': False,
                 'lim': (0.5, 3.5),
                 'cmap': 'cividis',
                 'label': r'$\Gamma$'},
        'HEAT': {'log': True,
                 'lim': (-1e31, 1e32),
                 'cmap': 'Spectral_r',
                 'label': r'$Q_\nu$ [erg]'},
        'DELP': {'log': True,
                 'lim': (1e4, 1e15),
                 'cmap': 'viridis',
                 'label': r''},
        'JX': {'log': True,
               'lim': (-1e24, 1e24),
               'cmap': 'coolwarm',
               'label': r'$j_r$ [g$\cdot$cm$^{2}\cdot$s$^{-2}$]'},
        'JY': {'log': True,
               'lim': (-1e20, 1e20),
               'cmap': 'Spectral_r',
               'label': r'$j_\theta$ [g$\cdot$cm$^{2}\cdot$s$^{-2}$]'},
        'JZ': {'log': True,
               'lim': (-1e23, 1e23),
               'cmap': 'seismic',
               'label': r'$j_\phi$ [g$\cdot$cm$^{2}\cdot$s$^{-2}]$'},
        'gas_pressure': {'log': True,
                 'lim': (1e25, 1e34),
                 'cmap': 'gist_rainbow_r',
                 'label': r'P$_\mathrm{gas}$ [dyn cm$^{-2}$]'},
        'soundspeed': {'log': True,
                   'lim': (1e8, 1e10),
                   'cmap': 'nipy_spectral',
                   'label': r'v$_\mathrm{sound}$ [cm$\cdot$s$^{-1}$]'},
        'alfven_velocity': {'log': True,
                   'lim': (1e4, 1e7),
                   'cmap': 'gnuplot',
                   'label': r'v$_\mathrm{a}$ [cm$\cdot$s$^{-1}$]'},
        'convective_velocity': {'log': True,
                   'lim': (-1e9, 1e9),
                   'cmap': 'Spectral_r',
                   'label': r'v$_\mathrm{conv}$ [cm$\cdot$s$^{-1}$]'},
        'turbulent_velocity': {'log': True,
                   'lim': (1e8, 1e10),
                   'cmap': 'cividis',
                   'label': r'v$_\mathrm{turb}$ [cm$\cdot$s$^{-1}$]'},
        'omega': {'log': True,
                   'lim': (1e0, 1e4),
                   'cmap': 'magma',
                   'label': r'$\Omega$ [s$^{-1}$]'},
        'gravitational_potential': {'log': True,
                  'lim': (-1e22, -1e15),
                   'cmap': 'magma',
                   'label': r'$\Phi$ [erg$\cdot$g$^{-1}$]'},
        'gravitational_energy': {'log': True,
                  'lim': (-1e22, -1e15),
                   'cmap': 'magma',
                   'label': r'E$_mathrm{grav}$ [erg]'},
        'neutron_fraction': {'log': False,
                'lim': (0, 1),
                'cmap': 'cividis',
                'label': r'X$_n$'},
        'proton_fraction': {'log': False,
                'lim': (0, 1),
                'cmap': 'viridis', 'label': r'X$_p$'},
        'alpha_fraction': {'log': False,
                    'lim': (0, 1), 
                    'cmap': 'plasma',
                    'label': r'X$_\alpha$'},
        'heavy_fraction': {'log': False,
                'lim': (0, 1),
                'cmap': 'magma',
                'label': r'X$_h$'},
        'Abar': {'log': True,
                 'lim': (4, 80),
                 'cmap': 'gist_rainbow_r',
                 'label': r'$\bar{A}$'},
        'Zbar': {'log': False,
                 'lim': (1, 34), 
                 'cmap': 'nipy_spectral', 
                 'label': r'$\bar{Z}$'},
        'electron_chemical_potential': {'log': True, 
                   'lim': (0.1, 300), 
                   'cmap': 'coolwarm', 
                   'label': r'$\mu_e$ [erg$\cdot$g$^{-1}$]'},
        'neutron_chemical_potential': {'log': True, 
                   'lim': (-2e2, 3e3), 
                   'cmap': 'bwr', 
                   'label': r'$\mu_n$ [erg$\cdot$g$^{-1}$]'},
        'proton_chemical_potential': {'log': True, 
                   'lim': (-2e2, 3e3), 
                   'cmap': 'seismic', 
                   'label': r'$\mu_p$ [erg$\cdot$g$^{-1}$]'},
        'neutrino_chemical_potential': {'log': True, 
                    'lim': (-2e2, 3e3), 
                    'cmap': 'Spectral_r', 
                    'label': r'$\mu_\nu$ [erg$\cdot$g$^{-1}$]'},
        'BHEX': {'log': True, 
                 'lim': (1e4, 1e15), 
                 'cmap': 'viridis', 
                 'label': r'$\rho$ [g cm$^{-3}$]'},
        'NUEX': {'log': True, 
                 'lim': (-1e39, 1e39), 
                 'cmap': 'Spectral_r', 
                 'label': r'F$_{\nu_e,r}$ [$\#_\nu\cdot$ s$^{-1}$]'},
        'NUEY': {'log': True, 
                 'lim': (-1e38, 1e38), 
                 'cmap': 'Spectral_r', 
                 'label': r'F$_{\nu_e,\theta}$ [$\#_\nu\cdot$ s$^{-1}$]'},
        'NUEZ': {'log': True, 
                 'lim': (-1e38, 1e38), 
                 'cmap': 'Spectral_r', 
                 'label': r'F$_{\nu_e,\phi}$ [$\#_\nu\cdot$ s$^{-1}$]'},
        'NUEE': {'log': True, 
                 'lim': (1e25, 1e34), 
                 'cmap': 'viridis', 
                 'label': r'E$_{\nu_e}$ [erg$\cdot$ s$^{-1}$]'},
        'NUAX': {'log': True, 
                 'lim': (-1e39, 1e39), 
                 'cmap': 'Spectral_r', 
                 'label': r'F$_{\overline{\nu}_e,r}$' \
                     r' [$\#_\nu\cdot$ s$^{-1}$]'},
        'NUAY': {'log': True, 'lim': (-1e38, 1e38),
                 'cmap': 'Spectral_r',
                 'label': r'F$_{\overline{\nu}_e,\theta}$' \
                     r' [$\#_\nu\cdot$ s$^{-1}$]'},
        'NUAZ': {'log': True,
                 'lim': (-1e38, 1e38),
                 'cmap': 'Spectral_r',
                 'label': r'F$_{\overline{\nu}_e,\phi}$' \
                     r' [$\#_\nu\cdot$ s$^{-1}$]'},
        'NUAE': {'log': True, 
                 'lim': (1e25, 1e34),
                 'cmap': 'viridis',
                 'label': r'E$_{\overline{\nu}_e}$ [erg$\cdot$ s$^{-1}$]'},
        'NUXX': {'log': True,
                 'lim': (-1e39, 1e39),
                 'cmap': 'Spectral_r',
                 'label': r'F$_{\nu_x,r}$ [$\#_\nu\cdot$ s$^{-1}$]'},
        'NUXY': {'log': True,
                 'lim': (-1e38, 1e38),
                 'cmap': 'Spectral_r', 
                 'label': r'F$_{\nu_x,\theta}$ [$\#_\nu\cdot$ s$^{-1}$]'},
        'NUXZ': {'log': True,
                 'lim': (-1e38, 1e38),
                 'cmap': 'Spectral_r',
                 'label': r'F$_{\nu_x,\phi}$ [$\#_\nu\cdot$ s$^{-1}$]'},
        'NUXE': {'log': True,
                 'lim': (1e25, 1e34),
                 'cmap': 'viridis',
                 'label': r'E$_{\nu_x}$ [erg$\cdot$ s$^{-1}$]'},
        'BX': {'log': True,
               'lim': (-1e15, 1e15),
               'cmap': 'coolwarm',
               'label': r'B$_r$ [G]'},
        'BY': {'log': True,
               'lim': (-1e15, 1e15),
               'cmap': 'coolwarm',
               'label': r'B$_\theta$ [G]'},
        'BZ': {'log': True,
               'lim': (-1e15, 1e15),
               'cmap': 'coolwarm',
               'label': r'B$_\phi$ [G]'},
        'poloidal_magnetic_fields': {'log': True,
               'lim': (1e11, 1e15),
               'cmap': 'inferno',
               'label': r'B$_\mathrm{pol}$ [G]'},
        'toroidal_magnetic_fields': {'log': True,
               'lim': (-1e12, 1e12),
               'cmap': 'coolwarm',
               'label': r'B$_\mathrm{tor}$ [G]'},
        'toroidal_magnetic_energy': {'log': True,
               'lim': (1e20, 1e28),
               'cmap': 'magma',
               'label': r'E$_\mathrm{mag, tor}$ [erg]'},
        'poloidal_magnetic_energy': {'log': True,
               'lim': (1e20, 1e28),
               'cmap': 'magma',
               'label': r'E$_\mathrm{mag, pol}$ [erg]'},
        'total_magnetic_energy': {'log': True,
               'lim': (1e20, 1e28),
               'cmap': 'magma',
               'label': r'E$_\mathrm{mag}$ [erg]'},
        'shock_radius': {'log': True,
                         'lim': (1e0, 1e4),
                         'label': r'R$_\mathrm{shock}$ [km]'},
        'PNS_radius': {'log': False,
                         'lim': (0, 150),
                         'label': r'R$_\mathrm{PNS}$ [km]'},
        'neutrino_spheres_nua': {'log': False,
                         'lim': (0,150),
                         'label': r'R$_{\overline{\nu}_e}$ [km]'},
        'neutrino_spheres_nue': {'log': False,
                         'lim': (0,150),
                         'label': r'R$_{\nu_e}$ [km]'},
        'neutrino_spheres_nux': {'log': False,
                         'lim': (0,150),
                         'label': r'R$_{\nu_x}$ [km]'},
        'gain_radius': {'log': True,
                         'lim': (1e0, 1e4),
                         'label': r'R$_\mathrm{gain}$ [km]'},
        'innercore_radius': {'log': True,
                         'lim': (1e0,1e4),
                         'label': r'R$_\mathrm{innercore}$ [km]'},
        'PNS_nucleus_radius': {'log': False,
                         'lim': (0, 40),
                         'label': r'R$_\mathrm{PNS, nuc}$ [km]'},
        'convective_flux': {'log': True,
                         'lim': (-1e40, 1e40),
                         'cmap': 'RdYlGn_r',
                         'label': r'F$_\mathrm{conv}$ [erg$\cdot$s$^{-1}$]'},
        'BV_frequency': {'log': True,
                         'lim': (-1e5, 1e5),
                         'cmap': 'RdYlBu_r',
                         'label': r'$\omega_\mathrm{BV}$ [s$^{-1}$]'},            
        'Rossby_number': {'log': True,
                         'lim': (-1e-4, 1e-4),
                         'cmap': 'RdYlBu_r',
                         'label': r'Ro'},
        'GW_Amplitudes_h+eq': {'log': False,
                         'lim': GW_limit,
                         'label': r'$\mathcal{D}h_+^{eq}$ [cm]'},
        'GW_Amplitudes_hxeq': {'log': False,
                         'lim': GW_limit, 
                         'label': r'$\mathcal{D}h_\times^{eq}$ [cm]'},
        'GW_Amplitudes_h+pol': {'log': False,
                         'lim': GW_limit, 
                         'label': r'$\mathcal{D}h_+^{pol}$ [cm]'},
        'GW_Amplitudes_hxpol': {'log': False,
                         'lim': GW_limit, 
                         'label': r'$\mathcal{D}h_\times^{pol}$ [cm]'},
        'mass_accretion_500km': {'log': False,
                         'lim': (0, 3),
                         'label': r'$\dot{M}_{500km}$ '\
                                r'[M$_\odot\cdot$s$^{-1}$]'},
         'explosion_mass': {'log': True,
                         'lim': (1e-4, 1e0),
                         'label': r'M$_\mathrm{exp}$ [M$_\odot$]'},
         'explosion_ene': {'log': True,
                         'lim': (1e47, 3e51),
                         'label': r'E$_\mathrm{exp}$ [erg]'},
         'explosion_kin' : {'log': True,
                         'lim': (1e47, 3e51),
                         'label': r'E$_\mathrm{kin}$ [erg]'},
         'explosion_mag': {'log': True,
                          'lim': (1e47, 3e51),
                          'label': r'E$_\mathrm{mag}$ [erg]'},
         'gain_mass': {'log': False,
                       'lim': (0, 1e0),
                       'label': r'M$_\mathrm{gain}$ [M$_\odot$]'},
         'gain_ene': {'log': False,
                      'lim': (-1e53, 3e51),
                      'label': r'Q$_\nu$ [erg]'},
         'innercore_mass': {'log': False,
                            'lim': (0, 2),
                            'label': r'M$_\mathrm{innercore}$ [M$_\odot$]'},
         'innercore_ene': {'log': True,
                            'lim': (1e47, 1e53),
                            'label': r'E$_\mathrm{innercore}$ [erg]'},
         'innercore_kin': {'log': True,
                            'lim': (1e47, 1e53),
                            'label': r'E$_\mathrm{kin}$ [erg]'},
         'innercore_mag': {'log': True,
                           'lim': (1e47, 1e52),
                           'label': r'E$_\mathrm{mag}$ [erg]'},
         'innercore_rot': {'log': True,
                           'lim': (1e47, 1e52),
                           'label': r'E$_\mathrm{rot}$ [erg]'},
         'innercore_grav': {'log': True,
                            'lim': (-1e54, -1e48),
                            'label': r'E$_\mathrm{grav}$ [erg]'},
         'innercore_T/W': {'log': False,
                           'lim': (0, 0.05),
                           'label': r'$T/|W|$'},
         'PNS_mass': {'log': False,
                      'lim': (0, 2),
                      'label': r'M$_\mathrm{PNS}$ [M$_\odot$]'},
         'PNS_ene': {'log': True,
                     'lim': (1e47, 1e52),
                     'label': r'E$_\mathrm{PNS}$ [erg]'},
         'PNS_kin': {'log': True,
                     'lim': (1e47, 1e52),
                     'label': r'E$_\mathrm{kin}$ [erg]'},
         'PNS_mag': {'log': True,
                     'lim': (1e47, 1e51),
                     'label': r'E$_\mathrm{mag}$ [erg]'},
         'PNS_rot': {'log': True,
                     'lim': (1e47, 1e51),
                     'label': r'E$_\mathrm{rot}$ [erg]'},
         'PNS_grav': {'log': True,
                      'lim': (-1e54, -1e49),
                      'label': r'E$_\mathrm{grav}$ [erg]'},
         'PNS_conv': {'log': False,
                      'lim': (1e47, 1e51),
                      'label': r'E$_\mathrm{conv}$ [erg]'},
         'nue_moment_x': {'log': True,
                      'lim': (-1e40, 1e40),
                      'cmap': 'PiYG_r',
                      'label': r'F$_{\nu_{e, r}}$ [erg$\cdot$s$^{-1}$]'},
         'nue_moment_y': {'log': True,
                      'lim': (-1e39, 1e39),
                      'cmap': 'PRGn_r',
                      'label': r'F$_{\nu_{e, \theta}}$ [erg$\cdot$s$^{-1}$]'},
         'nue_moment_z': {'log': True,
                      'lim': (-1e39, 1e39),
                      'cmap': 'BrBG_r',
                      'label': r'F$_{\nu_{e, \phi}}$ [erg$\cdot$s$^{-1}$]'},
         'nue_moment_e': {'log': True,
                       'lim': (1e28, 1e32),
                       'label': r'E$_{\nu_e}$ [erg]',
                       'cmap': 'viridis'},
         'nua_moment_x': {'log': True,
                      'lim': (-1e40, 1e40),
                      'cmap': 'PuOr_r',
                      'label': r'F$_{\overline{\nu}_{e, r}}$ '\
                             r'[erg$\cdot$s$^{-1}$]'},
         'nua_moment_y': {'log': True,
                      'lim': (-1e39, 1e39),
                      'cmap': 'RdGy_r',
                      'label':r'F$_{\{overline{\nu}_{e, \theta}}$ '\
                             r'[erg$\cdot$s$^{-1}$]'},
         'nua_moment_z': {'log': True,
                      'lim': (-1e39, 1e39),
                      'cmap': 'RdBu_r',
                      'label': r'F$_{\overline{\nu}_{e, \phi}}$ '\
                             '[erg$\cdot$s$^{-1}$]'},
         'nua_moment_e': {'log': True,
                       'lim': (1e28, 1e32),
                       'label': r'E$_{\{overline{\nu}_e}$ [erg]',
                       'cmap': 'magma'},
         'nux_moment_x': {'log': True,
                      'lim': (-1e40, 1e40),
                      'cmap': 'RdYlBu_r',
                      'label': r'F$_{\nu_{x, r}}$ [erg$\cdot$s$^{-1}$]'},
         'nux_moment_y': {'log': True,
                      'lim': (-1e39, 1e39),
                      'cmap': 'RdYlGn_r',
                      'label': r'F$_{\nu_{x, \theta}}$ [erg$\cdot$s$^{-1}$]'},
         'nux_moment_z': {'log': True,
                      'lim': (-1e39, 1e39),
                      'cmap': 'Spectral_r',
                      'label': r'F$_{\nu_{x, \phi}}$ [erg$\cdot$s$^{-1}$]'},
         'nux_moment_e': {'log': True,
                       'lim': (1e28, 1e32),
                       'label': r'E$_{\nu_x}$ [erg]',
                       'cmap': 'plasma'},
        'nue_mean_ene': {'log': False,
                         'lim': (10, 60),
                         'cmap': 'CMRmap',
                         'label': r'$\langle E_{\nu_e} \rangle$ [MeV]'},
        'nua_mean_ene': {'log': False,
                         'lim': (10, 50),
                         'cmap': 'terrain',
                         'label':\
                            r'$\langle E_{\overline{\nu}_e} \rangle$ [MeV]'},
        'nux_mean_ene': {'log': False,
                         'lim': (15, 60),
                         'cmap': 'brg',
                         'label': r'$\langle E_{\nu_x} \rangle$ [MeV]'},
        'nu_integrated_ene_all': {'log': False,
                                  'lim': (0, 30),
                                  'label': r'$\langle E_\nu\rangle$ [MeV]'},
        'nu_integrated_ene_nue': {'log': False,
                                  'lim': (0, 30),
                                  'label': r'$\langle E_{\nu_e}\rangle$ [MeV]'},
         'nu_integrated_ene_nua': {'log': False,
                                  'lim': (0, 30),
                                  'label': r'$\langle E_{\overline{\nu}_e}\rangle$ [MeV]'},
         'nu_integrated_ene_nux': {'log': False,
                                  'lim': (0, 30),
                                  'label': r'$\langle E_{\nu_x}\rangle$ [MeV]'},
         'nu_integrated_lum_all': {'log': False,
                                   'lim': (0, 1),
                                   'label': r'$L_\nu$ [$10^{53}$ erg$\cdot$s$^{-1}$]'},
         'nu_integrated_lum_nue': {'log': False,
                                   'lim': (0, 1),
                                   'label': r'$L_{\nu_e}$ [$10^{53}$ erg$\cdot$s$^{-1}$]'},
         'nu_integrated_lum_nua': {'log': False,
                                   'lim': (0, 1),
                                   'label': r'$L_{\overline{\nu}_e}$ [$10^{53}$ erg$\cdot$s$^{-1}$]'},
         'nu_integrated_lum_nux': {'log': False,
                                   'lim': (0, 1),
                                   'label': r'$L_{\nu_x}$ [$10^{53}$ erg$\cdot$s$^{-1}$]'}
}

## PLOT LABELS FOR X-AXIS
xaxis_labels = {'radius': 'R [km]',
                'theta': r'$\theta$ [rad]',
                'phi': r'$\phi$ [rad]',
                'time': 't [s]'
                }