from AeViz.AeVizMethods import *

## ---------------------------------------------------------------------
## SUPERNOVA
## ---------------------------------------------------------------------

"""
Module containing all the supernova plotting methods.
These functions are not meant to be used standalone, but rather to be
imported into the AeViz class.
"""

## RADII

@fig_window_open
def PNS_radius(self, rad:Literal['all', 'min', 'max', 'avg']='avg', **kwargs):
    self.plot1D(None, 'PNS_radius_' + comp, 'time', None, None, **kwargs)

@fig_window_open
def shock_radius(self, comp:Literal['all', 'min', 'max','avg']='avg',
                    **kwargs):
    self.plot1D(None, 'shock_radius_' + comp, 'time', None, None, **kwargs)

@fig_window_open
def neutrino_spheres(self, comp:Literal['all', 'min', 'max', 'avg']='avg',
                        **kwargs):
    for flavour in ['nue', 'nua', 'nux']:
        self.plot1D(None, 'neutrino_spheres_' + flavour + '_' + comp,
                    'time', None, None, **kwargs)

@fig_window_open
def gain_radius(self,  comp:Literal['all', 'min', 'max', 'avg']='avg',
                **kwargs):
    self.plot1D(None, 'gain_radius_' + comp, 'time', None, None, **kwargs)

@fig_window_open
def innercore_radius(self, comp:Literal['all', 'min', 'max','avg']='avg',
                        **kwargs):
    self.plot1D(None, 'innercore_radius_' + comp, 'time', None, None, **kwargs)

@fig_window_open
def PNS_nucleus_radius(self, comp:Literal['all', 'min', 'max',
                                            'avg']='avg', **kwargs):
    self.plot1D(None, 'PNS_nucleus_radius_' + comp, 'time', None, None,
                **kwargs)
    
## ENERGIES
@fig_window_open
def explosion(self, comp: Literal['all', 'mass', 'ene', 'kin', 'mag'],
                **kwargs):
    if comp == 'all':
        qt = ['mass', 'ene', 'kin', 'mag']
    else:
        qt = [comp]
    for q in qt:
        self.plot1D(None, 'explosion_' + q, 'time', None, None, **kwargs)

@fig_window_open
def gain(self, comp: Literal['all', 'mass', 'ene'], **kwargs):
    if comp == 'all':
        qt = ['mass', 'ene']
    else:
        qt = [comp]
    for q in qt:
        self.plot1D(None, 'gain_' + q, 'time', None, None, **kwargs)

@fig_window_open
def innercore(self, comp: Literal['mass', 'ene', 'kin', 'mag', 'rot',
                                    'grav', 'T/W'], **kwargs):
    self.plot1D(None, 'innercore_' + comp, 'time', None, None, **kwargs)

@fig_window_open
def PNS_nucleus(self, comp: Literal['mass', 'ene', 'kin', 'mag', 'rot',
                                    'grav', 'T/W'], **kwargs):
    self.plot1D(None, 'PNS_core_' + comp, 'time', None, None, **kwargs)

@fig_window_open
def PNS(self, comp: Literal['mass', 'ene', 'kin', 'mag', 'rot', 'grav',
                            'conv'], **kwargs):
    self.plot1D(None, 'PNS_' + comp, 'time', None, None, **kwargs)

## MASSES, ACCRETION AND VELOCITIES

@fig_window_open
def kick_velocity(self,
                    comp:Literal['modulus', 'nu', 'hydro', 'all']='all',
                    **kwargs):
    plot_qt(self, None, 'kick_velocity_'+comp, '1D', None, None, 'time')

@fig_window_open
def PNS_angular_mom(self, comp: Literal['all', 'Lx', 'Ly', 'Lz', 'Ltot'],
                    **kwargs):
    
    self.plot1D(None, 'PNS_angular_mom_' + comp, 'time', None, None,
                **kwargs)

@fig_window_open
def mass_accretion(self, **kwargs):
    self.plot1D(None, 'mass_accretion_500km', 'time', None, None, **kwargs)




