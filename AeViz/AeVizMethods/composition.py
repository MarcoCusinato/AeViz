from AeViz.AeVizMethods import *

## ---------------------------------------------------------------------
## COMPOSITION DATA
## ---------------------------------------------------------------------

"""
Module containing all the composition quantities plotting methods.
These functions are not meant to be used standalone, but rather to be
imported into the AeViz class.
"""

@fig_window_open
def Ye(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
        index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                        'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'Ye'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def Xn(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
        index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                        'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'neutron_fraction'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def Xp(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
        index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                        'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'proton_fraction'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def Xalpha(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
        index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                        'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'alpha_fraction'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def Xheavy(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
        index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                        'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'heavy_fraction'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def Abar(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
        index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                        'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'Abar'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def Zbar(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
        index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                        'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'Zbar'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def cpot_e(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
        index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                        'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'electron_chemical_potential'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def cpot_n(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
        index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                        'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'neutron_chemical_potential'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def cpot_p(self, file=None, projection:Literal['1D', '2D']='1D',
            index1=None, index2=None,
            plane:Literal['xz', 'yz', 'xy', 'radius', 'theta', 'phi',
                            'time']='radius', **kwargs):
    qt = 'proton_chemical_potential'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def cpot_nu(self, file=None, projection:Literal['1D', '2D']='1D',
            index1=None, index2=None,
            plane:Literal['xz', 'yz', 'xy', 'radius', 'theta', 'phi',
                            'time']='radius', **kwargs):
    qt = 'neutrino_chemical_potential'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)