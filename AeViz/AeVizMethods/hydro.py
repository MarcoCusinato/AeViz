from AeViz.AeVizMethods import *

## ---------------------------------------------------------------------
## HYDRODYNAMICAL DATA
## ---------------------------------------------------------------------

"""
Module containing all the hydrodynamical quantities plotting methods.
These functions are not meant to be used standalone, but rather to be
imported into the AeViz class.
"""

@fig_window_open
def rho(self, file=None, projection:Literal['1D', '2D']='1D', plane='time',
        **kwargs):
    AeViz_plot_panel(self, 'rho', file, projection, plane, **kwargs)

@fig_window_open
def MHD_energy(self, file=None, projection:Literal['1D', '2D']='1D',
                    plane='time', **kwargs):
    AeViz_plot_panel(self, 'MHD_energy', file, projection, plane, **kwargs)

@fig_window_open
def radial_velocity(self, file=None, projection:Literal['1D', '2D']='1D',
                    plane='time', **kwargs):
    AeViz_plot_panel(self, 'radial_velocity', file, projection, plane, **kwargs)

@fig_window_open
def theta_velocity(self, file=None, projection:Literal['1D', '2D']='1D',
                    plane='time', **kwargs):
    AeViz_plot_panel(self, 'theta_velocity', file, projection, plane, **kwargs)
    
@fig_window_open
def phi_velocity(self, file=None, projection:Literal['1D', '2D']='1D',
                    plane='time', **kwargs):
    AeViz_plot_panel(self, 'phi_velocity', file, projection, plane, **kwargs)

@fig_window_open
def soundspeed(self, file=None, projection:Literal['1D', '2D']='1D', plane='time',
            **kwargs):
    AeViz_plot_panel(self, 'soundspeed', file, projection, plane, **kwargs)

@fig_window_open
def omega(self, file=None, projection:Literal['1D', '2D']='1D', plane='time',
            **kwargs):
    AeViz_plot_panel(self, 'omega', file, projection, plane, **kwargs)