from AeViz.AeVizMethods import *

## ---------------------------------------------------------------------
## AUXILIARY METHODS
## ---------------------------------------------------------------------

"""
Module containing all the auxiliary methods plotting methods.
These functions are not meant to be used standalone, but rather to be
imported into the AeViz class.
"""

@fig_window_open
def add_field(self, plot, comp: Literal['velocity', 'Bfield'], **kwargs):
    """
    Adds velocity or magnetic field to the selected plot.
    """
    self.add_2Dfield(plot, comp)

@fig_window_open
def error(self, file=None, projection:Literal['1D', '2D']='1D',
          plane='time', **kwargs):
    AeViz_plot_panel(self, 'error', file, projection, plane, **kwargs)