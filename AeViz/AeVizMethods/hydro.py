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
def rho(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
        index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                        'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'rho'
    if projection == '1D':
        if plane == 'time':
            qt = 'rho_max'
        self.plot1D(file, qt, plane, index1, index2, **kwargs)
    elif projection == '2D':
        if plane == 'time':
            self.plotProfile(qt)
        else:
            self.plot2D(file, plane, qt, **kwargs)

@fig_window_open
def MHD_energy(self, file=None, projection:Literal['1D', '2D']='1D',
                index1=None, index2=None, plane:Literal['xz', 'yz', 'xy',
                    'radius', 'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'MHD_energy'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def velocity(self, file=None, projection:Literal['1D', '2D']='1D', 
                index1=None, index2=None, comp:Literal['all', 'r', 'th',
                'ph']='all', plane:Literal['xz', 'yz', 'xy', 'radius', 'theta',
                                        'phi', 'time']='radius', **kwargs):
    if comp == 'all':
        qt = ['radial_velocity', 'theta_velocity', 'phi_velocity']
    elif comp == 'r':
        qt = ['radial_velocity']
    elif comp == 'th':
        qt = ['theta_velocity']
    elif comp == 'ph':
        qt = ['phi_velocity']
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

@fig_window_open
def soundspeed(self, file=None, projection:Literal['1D', '2D']='1D',
                index1=None, index2=None, plane:Literal['xz', 'yz', 'xy',
                'radius', 'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'soundspeed'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)


@fig_window_open
def omega(self, file=None, projection:Literal['1D', '2D']='1D', index1=None,
          index2=None, plane:Literal['xz', 'yz', 'xy',
                                                    'radius',
                        'theta', 'phi', 'time']='radius', **kwargs):
    qt = 'omega'
    plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)