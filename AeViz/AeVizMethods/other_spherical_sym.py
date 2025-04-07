from AeViz.AeVizMethods import *

## ---------------------------------------------------------------------
## OTHER METHODS
## ---------------------------------------------------------------------

"""
Module containing all the other methods plotting methods.
These functions are not meant to be used standalone, but rather to be
imported into the AeViz class.
"""

## SPHERICAL HARMONICS

@fig_window_open
def rho_decomposition(self, l=0, m=None, correct_zero=True, **kwargs):
    loc = locals()
    for q in loc.keys():
        if q not in ['self', 'kwargs']:
            kwargs[q] = loc[q]
    qt = 'rho_spherical_harmonics'
    self.plotProfile(qt, **kwargs)

@fig_window_open
def rho_decomposition_barcode(self, lmin=None, lmax=None, msum=False,
                                r=None, rhomin=None, rhomax=None,
                                zero_norm=True, **kwargs):
    loc = locals()
    for q in loc.keys():
        if q not in ['self', 'kwargs']:
            kwargs[q] = loc[q]
    self.barcode(**kwargs)

## TIDAL DATA

@fig_window_open
def tidal_number(self, projection:Literal['1D', '2D']='1D',
                 type:Literal['love', 'lambda']='love',
                 comp:Literal['PNS_core', 'PNS']='PNS',
                 **kwargs):
    kwargs['comp'] = comp
    if type == 'love':
        qt = 'love_number'
    else:
        qt = 'tidal_deformability'
    AeViz_plot_panel(self, qt, None, projection, 'time', **kwargs)