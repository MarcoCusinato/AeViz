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
    AeViz_plot_panel(self, 'neutrino_mean_energy', file, projection,
                     plane, **kwargs)
    
@fig_window_open
def neutrino_number_density(self, file=None, projection:Literal['1D', '2D']='1D',
                        plane='radius', comp:Literal['nue', 'nua', 'nux']='nue',
                        **kwargs):
    kwargs['comp'] = comp
    AeViz_plot_panel(self, 'neutrino_number_density_grey', file, projection,
                     plane, **kwargs)
    
@fig_window_open
def neutrino_energy_density(self, file=None, projection:Literal['1D', '2D']='1D',
                plane='radius', comp:Literal['nue', 'nua', 'nux']='nue',
                **kwargs):
    kwargs['comp'] = comp
    AeViz_plot_panel(self, 'neutrino_momenta_grey', file, projection,
                     plane, **kwargs)

@fig_window_open
def global_neutrino_luminosities(self, projection:Literal['1D', '2D']='1D',
                            comp:Literal['nue', 'nua', 'nux', 'all']='all',
                            **kwargs):
    kwargs['comp'] = comp
    if comp == 'all':
        self.set_simple_labelling()
        for cm in ['nue', 'nua', 'nux']:
            kwargs['comp'] = cm
            AeViz_plot_panel(self, 'global_neutrino_luminosity', None, projection,
                            'time', **kwargs)
        self.set_simple_labelling()
    else:
        kwargs['comp'] = comp
        AeViz_plot_panel(self, 'global_neutrino_luminosity', None, projection,
                        'time', **kwargs)
    
@fig_window_open
def global_neutrino_number_luminosities(self, projection:Literal['1D', '2D']='1D',
                            comp:Literal['nue', 'nua', 'nux', 'all']='all',
                            **kwargs):
    kwargs['comp'] = comp
    if comp == 'all':
        self.set_simple_labelling()
        for cm in ['nue', 'nua', 'nux']:
            kwargs['comp'] = cm
            AeViz_plot_panel(self, 'global_neutrino_number_luminosity', None, projection,
                            'time', **kwargs)
        self.set_simple_labelling()
    else:
        kwargs['comp'] = comp
        AeViz_plot_panel(self, 'global_neutrino_number_luminosity', None, projection,
                        'time', **kwargs)
        
@fig_window_open
def global_neutrino_mean_energies(self, projection:Literal['1D', '2D']='1D',
                            comp:Literal['nue', 'nua', 'nux', 'all']='all',
                            **kwargs):
    kwargs['comp'] = comp
    if comp == 'all':
        self.set_simple_labelling()
        for cm in ['nue', 'nua', 'nux']:
            kwargs['comp'] = cm
            AeViz_plot_panel(self, 'global_neutrino_mean_energies', None, projection,
                            'time', **kwargs)
        self.set_simple_labelling()
    else:
        kwargs['comp'] = comp
        AeViz_plot_panel(self, 'global_neutrino_mean_energies', None, projection,
                        'time', **kwargs)

@fig_window_open
def neutrino_event_rate(self, projection:Literal['1D', '2D']='1D', **kwargs):
    AeViz_plot_panel(self, 'neutrino_event_rate', None, projection,
                        'time', **kwargs)