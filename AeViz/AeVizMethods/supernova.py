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
def PNS_radius(self, rad:Literal['all', 'min', 'max', 'avg']='avg',
               projection:Literal['1D', '2D']='1D', **kwargs):
    AeViz_plot_radius_panel(self, 'PNS_radius', projection, rad, **kwargs)

@fig_window_open
def shock_radius(self, rad:Literal['all', 'min', 'max', 'avg']='avg',
                 projection:Literal['1D', '2D']='1D', **kwargs):
    AeViz_plot_radius_panel(self, 'shock_radius', projection, rad, **kwargs)

@fig_window_open
def neutrino_spheres(self, rad:Literal['all', 'min', 'max', 'avg']='avg',
                     comp:Literal['nue', 'nua', 'nux']='nue',
                     projection:Literal['1D', '2D']='1D', **kwargs):
    kwargs['comp'] = comp
    AeViz_plot_radius_panel(self, 'neutrino_spheres', projection, rad, **kwargs)

@fig_window_open
def gain_radius(self, rad:Literal['all', 'min', 'max', 'avg']='avg',
                 projection:Literal['1D', '2D']='1D', **kwargs):
    AeViz_plot_radius_panel(self, 'gain_radius', projection, rad, **kwargs)

@fig_window_open
def innercore_radius(self, rad:Literal['all', 'min', 'max', 'avg']='avg',
                 projection:Literal['1D', '2D']='1D', **kwargs):
    AeViz_plot_radius_panel(self, 'innercore_radius', projection, rad, **kwargs)

@fig_window_open
def PNS_nucleus_radius(self, rad:Literal['all', 'min', 'max', 'avg']='avg',
                 projection:Literal['1D', '2D']='1D', **kwargs):
    AeViz_plot_radius_panel(self, 'PNS_nucleus_radius', projection, rad,
                            **kwargs)

@fig_window_open
def isodensities_lines(self, rad:Literal['all', 'min', 'max', 'avg']='avg',
                     comp:Literal['1e+14', '1e+13', '1e+12', '1e+11',
                                    '1e+10', '1e+09', '1e+08']='1e+14',
                     projection:Literal['1D', '2D']='1D', **kwargs):
    kwargs['comp'] = comp
    AeViz_plot_radius_panel(self, 'isodensities_lines', projection, rad,
                            **kwargs)
    
## ENERGIES
@fig_window_open
def explosion(self, projection:Literal['1D', '2D']='1D',
              comp:Literal['mass', 'tot', 'kin', 'mag', 'ratio']='mass',
              **kwargs):
    kwargs['comp'] = comp
    AeViz_plot_panel(self, 'explosion_mass_ene', None, projection, 'time',
                     **kwargs)

@fig_window_open
def gain(self,  projection:Literal['1D', '2D']='1D',
         comp: Literal['mass', 'nu_heath']='mass', **kwargs):
    kwargs['comp'] = comp
    AeViz_plot_panel(self, 'gain_mass_nu_heat', None, projection, 'time',
                     **kwargs)

@fig_window_open
def innercore(self,  projection:Literal['1D', '2D']='1D',
              comp:Literal['mass', 'kin', 'mag', 'rot', 'grav', 'tot', 'T/W']='mass',
              **kwargs):
    kwargs['comp'] = comp
    AeViz_plot_panel(self, 'innercore_mass_ene', None, projection, 'time',
                     **kwargs)

@fig_window_open
def PNS_nucleus(self,  projection:Literal['1D', '2D']='1D',
                comp: Literal['mass', 'ene', 'kin', 'mag',
                              'rot', 'grav', 'T/W']='mass', **kwargs):
    kwargs['comp'] = comp
    AeViz_plot_panel(self, 'PNS_core_mass_ene', None, projection, 'time',
                     **kwargs)

@fig_window_open
def PNS(self, projection:Literal['1D', '2D']='1D',
        comp:Literal['mass', 'kin', 'mag', 'rot',
                                    'conv', 'grav', 'tot', 'T/W']='mass',
        **kwargs):
    kwargs['comp'] = comp
    AeViz_plot_panel(self, 'PNS_mass_ene', None, projection, 'time',
                     **kwargs)

## MASSES, ACCRETION AND VELOCITIES

@fig_window_open
def kick_velocity(self, projection:Literal['1D', '2D']='1D',
                 comp:Literal['x', 'y', 'z', 'tot']='tot',
                 flavour:Literal['hydro', 'nu', 'nue', 'nua', 'nux', 'all']='all',
                 **kwargs):
    kwargs['comp'] = comp
    kwargs['flavour'] = flavour
    AeViz_plot_panel(self, 'PNS_kick_velocity', None, projection, 'time',
                     **kwargs)

@fig_window_open
def PNS_angular_mom(self, projection:Literal['1D', '2D']='1D',
                    comp: Literal['Lx', 'Ly', 'Lz', 'Ltot']='Ltot',
                    **kwargs):
    kwargs['comp'] = comp
    AeViz_plot_panel(self, 'PNS_angular_mom', None, projection, 'time',
                     **kwargs)

@fig_window_open
def PNS_angular_mom_nu(self, projection:Literal['1D', '2D']='1D',
                    comp: Literal['Lx', 'Ly', 'Lz', 'Ltot']='Ltot',
                    flavour:Literal['nue', 'nua', 'nux']='nue',
                    **kwargs):
    kwargs['comp'] = comp
    kwargs['flavour'] = flavour
    AeViz_plot_panel(self, 'PNS_angular_momentum_neutrinos', None, projection,
                     'time', **kwargs)

@fig_window_open
def mass_accretion(self, projection:Literal['1D', '2D']='1D', **kwargs):
    AeViz_plot_panel(self, 'mass_accretion_500km', None, projection, 'time',
                     **kwargs)

@fig_window_open
def mass_flux(self, file=None, projection:Literal['1D', '2D']='1D', plane='time',
        **kwargs):
    AeViz_plot_panel(self, 'mass_flux', file, projection, plane, **kwargs)


## NEUTRINO LUMINOSITY

@fig_window_open
def neutrino_luminosity(self, file=None, projection:Literal['1D', '2D']='1D', \
                        comp:Literal['nue', 'nua', 'nux', 'all']='all', \
                        plane='time', **kwargs):
    kwargs['comp'] = comp
    if comp == 'all':
        self.set_simple_labelling()
        for cm in ['nue', 'nua', 'nux']:
            kwargs['comp'] = cm
            AeViz_plot_panel(self, 'total_neutrino_luminosity', file, projection, \
                        plane, **kwargs)
        self.set_simple_labelling()
    else:
        kwargs['comp'] = comp
        AeViz_plot_panel(self, 'total_neutrino_luminosity', file, projection, \
                        plane, **kwargs)

