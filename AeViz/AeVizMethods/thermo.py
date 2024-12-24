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
def gas_pressure(self, file=None, projection:Literal['1D', '2D']='1D',
                    index1=None, index2=None, plane:Literal['xz', 'yz', 'xy',
                    'radius', 'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'gas_pressure'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def temperature(self, file=None, projection:Literal['1D', '2D']='1D',
                index1=None, index2=None, plane:Literal['xz', 'yz', 'xy',
                'radius', 'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'temperature'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def enthalphy(self, file=None, projection:Literal['1D', '2D']='1D',
                index1=None, index2=None, plane:Literal['xz', 'yz', 'xy',
                'radius', 'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'enthalphy'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def entropy(self, file=None, projection:Literal['1D', '2D']='1D',
            index1=None, index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
            'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'entropy'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def adiabatic_index(self, file=None, projection:Literal['1D', '2D']='1D',
                    index1=None, index2=None, plane:Literal['xz', 'yz', 'xy',
                    'radius', 'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'adiabatic_index'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def lorentz_factor(self, file=None, projection:Literal['1D', '2D']='1D',
                    index1=None, index2=None,
                    plane:Literal['xz', 'yz', 'xy', 'radius', 'theta', 
                                    'phi', 'time']='radius', **kwargs):
    qt = 'lorentz'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def grav_pot(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
        index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                        'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'gravitational_potential'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def grav_ene(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
        index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                        'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'gravitational_energy'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def int_ene(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
        index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                        'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'internal_energy'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def nu_heat(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
        index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                        'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'nu_heat'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)