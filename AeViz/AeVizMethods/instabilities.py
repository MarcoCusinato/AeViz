from AeViz.AeVizMethods import *

## ---------------------------------------------------------------------
## INSTABILITIES
## ---------------------------------------------------------------------

"""
Module containing all the instabilities plotting methods.
These functions are not meant to be used standalone, but rather to be
imported into the AeViz class.
"""

@fig_window_open
def convective_velocity(self, file=None, projection:Literal['1D', '2D']='1D',
                 plane='radius', **kwargs):
    AeViz_plot_panel(self, 'convective_velocity', file, projection, plane,
                     **kwargs)

@fig_window_open
def turbulent_velocity(self, file=None, projection:Literal['1D', '2D']='1D',
                 plane='radius', **kwargs):
    AeViz_plot_panel(self, 'turbulent_velocity', file, projection, plane,
                     **kwargs)

@fig_window_open
def convective_flux(self, file=None, projection:Literal['1D', '2D']='1D',
                 plane='time', **kwargs):
    AeViz_plot_panel(self, 'convective_flux', file, projection, plane,
                     **kwargs)

@fig_window_open
def Rossby_number(self, file=None, projection:Literal['1D', '2D']='1D',
                  plane='time', **kwargs):
    AeViz_plot_panel(self, 'Rossby_number', file, projection, plane, **kwargs)

@fig_window_open
def BV_frequency(self, file=None, projection:Literal['1D', '2D']='1D',
                 plane='time', **kwargs):
    AeViz_plot_panel(self, 'BV_frequency', file, projection, plane, 
                     **kwargs)

@fig_window_open
def epicyclic_frequency(self, file=None, projection:Literal['1D', '2D']='1D',
                        plane='time', **kwargs):
    AeViz_plot_panel(self, 'epicyclic_frequency', file, projection, plane,
                     **kwargs)