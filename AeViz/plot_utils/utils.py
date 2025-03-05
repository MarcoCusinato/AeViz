from AeViz.plot_utils.limits_utils import GW_limit, max_min
## PLOTTING UTILS DICTIONARIES

## PLOT LABELS, COLORMAPS, AND Y-LIMITS or COLORBAR LIMITS and LIMITS
## SCALES (LOG, LINEAR, or SYMLOG)
plot_labels = {
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
        'YN': {'log': False,
               'lim': (0, 1.0),
               'cmap': 'gist_rainbow',
               'label': r'Y$_N$'},

        'DENS': {'log': True,
                 'lim': (1e4, 1e15),
                 'cmap': 'viridis',
                 'label': r'$\rho$ [g cm$^{-3}$]'},
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
         'rho_spherical_harmonics': {'log': False,
                                    'lim': max_min,
                                    'label': r'$\tilde{\rho}_{XXYY}$',
                                    'cmap': 'YlGn'},
         'love_number': {'log': False,
                             'lim': (0, 0.002),
                             'label': r'$\kappa_2$'},
          'tidal_deformability': {'log': True,
                                       'lim': (1, 10000),
                                       'label': r'$\Lambda$'},
}

keys = list(plot_labels.keys())
for key in keys:
       if 'cmap' in plot_labels[key]:
              plot_labels['dr_' + key] = {'log': plot_labels[key]['log'],
                                           'lim': (None, None),
                                           'cmap': plot_labels[key]['cmap'],
                                           'label':r'$\partial_r$ ' + \
                                            plot_labels[key]['label'].split(']')[0] + '/cm]',
                                            }
              plot_labels['dtheta_' + key] = {'log': plot_labels[key]['log'],
                                              'lim': (None, None),
                                              'cmap': plot_labels[key]['cmap'],
                                              'label': r'$\partial_\theta$ ' + \
                                              plot_labels[key]['label'],
                                              }
              plot_labels['dphi_' + key] = {'log': plot_labels[key]['log'],
                                            'lim': (None, None),
                                            'cmap': plot_labels[key]['cmap'],
                                            'label': r'$\partial_\phi$ ' + \
                                            plot_labels[key]['label'],
                                            }
              plot_labels['dt_' + key] = {'log': plot_labels[key]['log'],
                                          'lim': (None, None),
                                          'cmap': plot_labels[key]['cmap'],
                                          'label': r'$\partial_t$ ' + \
                                          plot_labels[key]['label'].split(']')[0] + '/s]',
                                          }
       else:
              plot_labels['dr_' + key] = {'log': plot_labels[key]['log'],
                                          'lim': (None, None),
                                          'label':r'$\partial_r$ ' + \
                                          plot_labels[key]['label'].split(']')[0] + '/cm]',
                                          }
              plot_labels['dtheta_' + key] = {'log': plot_labels[key]['log'],
                                              'lim': (None, None),
                                              'label': r'$\partial_\theta$ ' + \
                                              plot_labels[key]['label'],
                                             }
              plot_labels['dphi_' + key] = {'log': plot_labels[key]['log'],
                                            'lim': (None, None),
                                            'label': r'$\partial_\phi$ ' + \
                                                   plot_labels[key]['label'],
                                           }
              plot_labels['dt_' + key] = {'log': plot_labels[key]['log'],
                                          'lim': (None, None),
                                          'label': r'$\partial_t$ ' + \
                                                  plot_labels[key]['label'].split(']')[0] + '/s]',
                                          }

## PLOT LABELS FOR X-AXIS
xaxis_labels = {'radius': 'R [km]',
                'theta': r'$\theta$ [rad]',
                'phi': r'$\phi$ [rad]',
                'time': 't [s]'
                }