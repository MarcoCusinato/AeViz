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
def add_field(self, file, plot, comp: Literal['velocity', 'Bfield'],
                plane: Literal['xy', 'yz', 'xz']='xz', **kwargs):
    """
    Adds velocity or magnetic field to the selected plot.
    """
    self.add_2Dfield(file, plot, comp, plane)

@fig_window_open
def error(self, file=None, projection:Literal['1D', '2D']='1D',
            index1=None, index2=None,
            plane:Literal['xz', 'yz', 'xy', 'radius', 'theta', 'phi',
                        'time']='radius', **kwargs):
    qt = 'error'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)