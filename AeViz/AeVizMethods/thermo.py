from AeViz.AeVizMethods import *
## ---------------------------------------------------------------------
## THERMODYNAMICAL DATA
## ---------------------------------------------------------------------

"""
Module containing all the thermodynamical quantities plotting methods.
These functions are not meant to be used standalone, but rather to be
imported into the AeViz class.
"""

@fig_window_open
def gas_pressure(self, file=None, projection:Literal['1D', '2D']='1D', plane='time',
        **kwargs):
   AeViz_plot_panel(self, 'gas_pressure', file, projection, plane, **kwargs)

@fig_window_open
def temperature(self, file=None, projection:Literal['1D', '2D']='1D', plane='time',
        **kwargs):
    AeViz_plot_panel(self, 'temperature', file, projection, plane, **kwargs)

@fig_window_open
def enthalphy(self, file=None, projection:Literal['1D', '2D']='1D', plane='time',
        **kwargs):
    AeViz_plot_panel(self, 'enthalphy', file, projection, plane, **kwargs)

@fig_window_open
def entropy(self, file=None, projection:Literal['1D', '2D']='1D', plane='time',
        **kwargs):
    AeViz_plot_panel(self, 'entropy', file, projection, plane, **kwargs)

@fig_window_open
def adiabatic_index(self, file=None, projection:Literal['1D', '2D']='1D', plane='time',
        **kwargs):
    AeViz_plot_panel(self, 'adiabatic_index', file, projection, plane, **kwargs)

@fig_window_open
def lorentz_factor(self, file=None, projection:Literal['1D', '2D']='1D', plane='time',
        **kwargs):
    AeViz_plot_panel(self, 'lorentz', file, projection, plane, **kwargs)

@fig_window_open
def gravitational_potential(self, file=None, projection:Literal['1D', '2D']='1D',
                            plane='time', **kwargs):
    AeViz_plot_panel(self, 'gravitational_potential', file, projection, plane, **kwargs)

@fig_window_open
def gravitational_energy(self, file=None, projection:Literal['1D', '2D']='1D',
                            plane='time', **kwargs):
    AeViz_plot_panel(self, 'gravitational_energy', file, projection, plane, **kwargs)

@fig_window_open
def int_ene(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
        index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                        'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'internal_energy'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def nu_heat(self, file=None, projection:Literal['1D', '2D']='1D', plane='time',
            **kwargs):
    AeViz_plot_panel(self, 'nu_heat', file, projection, plane, **kwargs)