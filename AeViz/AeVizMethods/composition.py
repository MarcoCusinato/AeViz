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
def Ye(self, file=None, projection:Literal['1D', '2D']='1D', plane='time',
       **kwargs):
    AeViz_plot_panel(self, 'Ye', file, projection, plane, **kwargs)

@fig_window_open
def Xn(self, file=None, projection:Literal['1D', '2D']='1D', plane='time',
       **kwargs):
    AeViz_plot_panel(self, 'neutron_fraction', file, projection, plane,
                     **kwargs)

@fig_window_open
def Xp(self, file=None, projection:Literal['1D', '2D']='1D', plane='time',
       **kwargs):
    AeViz_plot_panel(self, 'proton_fraction', file, projection, plane,
                     **kwargs)

@fig_window_open
def Xalpha(self, file=None, projection:Literal['1D', '2D']='1D', plane='time',
       **kwargs):
    AeViz_plot_panel(self, 'alpha_fraction', file, projection, plane, **kwargs)

@fig_window_open
def Xheavy(self, file=None, projection:Literal['1D', '2D']='1D', plane='time',
       **kwargs):
    AeViz_plot_panel(self, 'heavy_fraction', file, projection, plane, **kwargs)

@fig_window_open
def Abar(self, file=None, projection:Literal['1D', '2D']='1D', plane='time',
       **kwargs):
    AeViz_plot_panel(self, 'Abar', file, projection, plane, **kwargs)

@fig_window_open
def Zbar(self, file=None, projection:Literal['1D', '2D']='1D', plane='time',
       **kwargs):
    AeViz_plot_panel(self, 'Zbar', file, projection, plane, **kwargs)

@fig_window_open
def cpot_e(self, file=None, projection:Literal['1D', '2D']='1D', plane='time',
       **kwargs):
    AeViz_plot_panel(self, 'electron_chemical_potential', file, projection,
                     plane, **kwargs)

@fig_window_open
def cpot_n(self, file=None, projection:Literal['1D', '2D']='1D', plane='time',
       **kwargs):
    AeViz_plot_panel(self, 'neutron_chemical_potential', file, projection,
                     plane, **kwargs)

@fig_window_open
def cpot_p(self, file=None, projection:Literal['1D', '2D']='1D', plane='time',
       **kwargs):
    AeViz_plot_panel(self, 'proton_chemical_potential', file, projection,
                     plane, **kwargs)

@fig_window_open
def cpot_nu(self, file=None, projection:Literal['1D', '2D']='1D', plane='time',
       **kwargs):
    AeViz_plot_panel(self, 'neutrino_chemical_potential', file, projection,
                     plane, **kwargs)