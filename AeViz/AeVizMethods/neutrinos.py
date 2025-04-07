from AeViz.AeVizMethods import *

## ---------------------------------------------------------------------
## NEUTRINO DATA
## ---------------------------------------------------------------------

"""
Module containing all the neutrino quantities plotting methods.
These functions are not meant to be used standalone, but rather to be
imported into the AeViz class.
"""

@fig_window_open
def neutrino_energy_density(self, file=None, projection:Literal['1D', '2D']='1D',
                            plane='radius', comp:Literal['nue', 'nua', 'nux']='nue',
                            **kwargs):
    kwargs['comp'] = comp
    AeViz_plot_panel(self, 'neutrino_energy_density_grey', file, projection,
                     plane, **kwargs)

@fig_window_open
def neutrino_mean_energy(self, file=None, projection:Literal['1D', '2D']='1D',
                        plane='radius', comp:Literal['nue', 'nua', 'nux']='nue',
                        **kwargs):
    kwargs['comp'] = comp
    AeViz_plot_panel(self, 'neutrino_energy_density_grey', file, projection,
                     plane, **kwargs)
    
@fig_window_open
def neutrino_number_density(self, file=None, projection:Literal['1D', '2D']='1D',
                        plane='radius', comp:Literal['nue', 'nua', 'nux']='nue',
                        **kwargs):
    kwargs['comp'] = comp
    AeViz_plot_panel(self, 'neutrino_number_density_grey', file, projection,
                     plane, **kwargs)

@fig_window_open
def nue_moment(self, file=None, projection:Literal['1D', '2D']='1D',
                comp:Literal['x', 'y', 'z', 'e', 'all']='all', index1=None,
                index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                'theta', 'phi', 'time']='radius', **kwargs):
    if comp == 'all':
        qt = ['nue_moment_x', 'nue_moment_y', 'nue_moment_z',
                'nue_moment_e']
    else:
        qt = ['nue_moment_' + comp]
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
@fig_window_open
def nua_moment(self, file=None, projection:Literal['1D', '2D']='1D',
                comp:Literal['x', 'y', 'z', 'e', 'all']='all', index1=None,
                index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                'theta', 'phi', 'time']='radius', **kwargs):
    if comp == 'all':
        qt = ['nua_moment_x', 'nua_moment_y', 'nua_moment_z',
                'nua_moment_e']
    else:
        qt = ['nua_moment_' + comp]
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def nux_moment(self, file=None, projection:Literal['1D', '2D']='1D',
                comp:Literal['x', 'y', 'z', 'e', 'all']='all', index1=None,
                index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                'theta', 'phi', 'time']='radius', **kwargs):
    if comp == 'all':
        qt = ['nux_moment_x', 'nux_moment_y', 'nux_moment_z',
                'nux_moment_e']
    else:
        qt = ['nux_moment_' + comp]
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
@fig_window_open
def nu_mean_ene(self, file=None, projection:Literal['1D', '2D']='1D',
                comp:Literal['nue', 'nua', 'nux', 'all']='all',
                index1=None, index2=None,
                plane:Literal['xz', 'yz', 'xy', 'radius', 'theta', 'phi',
                                'time']='radius', **kwargs):
    if comp == 'all':
        qt = ['nue_mean_ene', 'nua_mean_ene', 'nux_mean_ene']
    else:
        qt = [comp + '_mean_ene']
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def nu_integrated(self, flavour:Literal['all', 'nue', 'nua', 'nux']='all',
                    comp: Literal['ene', 'lum']='lum', **kwargs):
    qt = 'nu_integrated_' + comp + '_' + flavour
    self.plot1D(None, qt, 'time', None, None, **kwargs)