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
def convective_velocity(self, file=None,
                        comp:Literal['all', 'convective', 'turbulent']='all',
                        **kwargs):
    if comp == 'all':
        qt = ['convective_velocity', 'turbulent_velocity']
    elif comp == 'convective':
        qt = ['convective_velocity']
    elif comp == 'turbulent':
        qt = ['turbulent_velocity']

    for q in qt:
        self.plot1D(file, q, 'radius', None, None, **kwargs)


@fig_window_open
def convective_flux(self, file=None,
                    plane:Literal['time', 'radius']='radius', **kwargs):
    qt = 'convective_flux'
    if plane == 'time':
        self.plotProfile(qt)
    else:
        self.plot1D(file, qt, plane, None, None, **kwargs)

@fig_window_open
def Rossby_number(self, file=None, projection:Literal['1D', '2D']='1D',
                    index1=None, index2=None, plane:Literal['xz', 'yz', 'xy',
                    'radius', 'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'Rossby_number'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def BV_frequency(self, file=None, projection:Literal['1D', '2D']='1D',
                    index1=None, index2=None,
                    plane:Literal['xz', 'yz', 'xy', 'radius', 'theta',
                                    'phi', 'time']='radius', **kwargs):
    qt = 'BV_frequency'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def epicyclic_frequency(self, file=None, projection:Literal['1D', '2D']='1D',
                        index1=None, index2=None,
                        plane:Literal['xz', 'yz', 'xy', 'radius', 'theta',
                                        'phi', 'time']='radius', **kwargs):
    qt = 'epicyclic_frequency'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)