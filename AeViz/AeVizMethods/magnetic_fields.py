from AeViz.AeVizMethods import *

## ---------------------------------------------------------------------
## MAGNETIC FIELD DATA
## ---------------------------------------------------------------------

"""
Module containing all the magnetic field quantities plotting methods.
These functions are not meant to be used standalone, but rather to be
imported into the AeViz class.
"""

@fig_window_open
def alfven_velocity(self, file=None, projection:Literal['1D', '2D']='1D',
                    index1=None, index2=None, plane:Literal['xz', 'yz', 'xy',
                            'radius', 'theta', 'phi', 'time']='radius',
                    **kwargs):
    qt = 'alfven_velocity'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def magnetic_fields(self, file=None, projection:Literal['1D', '2D']='1D',
                    index1=None, index2=None,
                    comp:Literal['all', 'r', 'th', 'ph', 'torpol',
                                    'poloidal', 'toroidal']='all',
                    plane:Literal['xz', 'yz', 'xy', 'radius', 'theta',
                                    'phi', 'time']='radius', **kwargs):
    if comp == 'all':
        qt = ['BX', 'BY', 'BZ']
    elif comp == 'r':
        qt = ['BX']
    elif comp == 'th':
        qt = ['BY']
    elif comp == 'ph':
        qt = ['BZ']
    elif comp == 'poloidal':
        qt = ['poloidal_magnetic_fields']
    elif comp == 'toroidal':
        qt = ['toroidal_magnetic_fields']
    elif comp == 'torpol':
        qt = ['poloidal_magnetic_fields', 'toroidal_magnetic_fields']
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def magnetic_energy(self, file=None, projection:Literal['1D', '2D']='1D', 
                    index1=None, index2=None,
                    comp:Literal['all', 'tot', 'poloidal',
                                    'toroidal']='all', 
                    plane:Literal['xz', 'yz', 'xy', 'radius', 'theta',
                                    'phi', 'time']='radius', **kwargs):
    if comp == 'all':
        qt = ['total_magnetic_energy', 'poloidal_magnetic_energy',
                'toroidal_magnetic_energy']
    elif comp == 'tot':
        qt = ['total_magnetic_energy']
    elif comp == 'poloidal':
        qt = ['poloidal_magnetic_energy']
    elif comp == 'toroidal':
        qt = ['toroidal_magnetic_energy']
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

