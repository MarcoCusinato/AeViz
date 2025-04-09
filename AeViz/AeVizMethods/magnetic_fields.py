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
                    plane='radius', **kwargs):
    AeViz_plot_panel(self, 'alfven_velocity', file, projection, plane, **kwargs)

@fig_window_open
def magnetic_fields(self, file=None, projection:Literal['1D', '2D']='1D',
                    comp: Literal['all', 'pol', 'tor', 'r', 'th', 'ph']='all',
                    plane='radius', **kwargs):
    if comp == 'all':
        AeViz_plot_panel(self, 'poloidal_magnetic_fields', file, projection,
                         plane, **kwargs)
        AeViz_plot_panel(self, 'toroidal_magnetic_fields', file, projection,
                         plane, **kwargs)
    elif comp == 'pol':
        AeViz_plot_panel(self, 'poloidal_magnetic_fields', file, projection,
                         plane, **kwargs)
    elif comp == 'tor':
        AeViz_plot_panel(self, 'toroidal_magnetic_fields', file, projection,
                         plane, **kwargs)
    else:
        kwargs['comp'] = comp
        AeViz_plot_panel(self, 'poloidal_magnetic_fields', file, projection,
                         plane, **kwargs)

@fig_window_open
def magnetic_energy(self, file=None, projection:Literal['1D', '2D']='1D',
                    comp: Literal['all', 'pol', 'tor', 'tot']='tot',
                    plane='radius', **kwargs):
    if comp == 'all':
        kwargs['comp'] = 'pol'
        AeViz_plot_panel(self, 'magnetic_energy', file, projection,
                         plane, **kwargs)
        kwargs['comp'] = 'tor'
        AeViz_plot_panel(self, 'magnetic_energy', file, projection,
                         plane, **kwargs)
    else:
        kwargs['comp'] = comp
        AeViz_plot_panel(self, 'magnetic_energy', file, projection,
                         plane, **kwargs)