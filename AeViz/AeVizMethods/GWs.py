from AeViz.AeVizMethods import *
from AeViz.simulation.methods.GWs import GW_Amplitudes

## ---------------------------------------------------------------------
## GRAVITATIONAL WAVES
## ---------------------------------------------------------------------

"""
Module containing all the gravitational waves plotting methods.
These functions are not meant to be used standalone, but rather to be
imported into the AeViz class.
"""

@fig_window_open
def GWs(self, projection:Literal['1D', '2D']='1D',
        comp:Literal['all', 'h+eq', 'h+pol', 'hxeq', 'hxpol']='h+eq',
        decomposition=False, spectrogram=False, **kwargs):
    if self.sim_dim == 1:
        pass
    kwargs['spectrogram'] = spectrogram
    if self.sim_dim == 2:
        kwargs['comp'] = 'h+eq'
        AeViz_plot_panel(self, 'GW_Amplitudes', None, projection, 'time',
                         **kwargs)
    else:
        if comp == 'all':
            for cm in ['h+eq', 'h+pol', 'hxeq', 'hxpol']:
                kwargs['comp'] = cm
                AeViz_plot_panel(self, 'GW_Amplitudes', None, projection, 
                                 'time', **kwargs)
        else:
            kwargs['comp'] = comp
            AeViz_plot_panel(self, 'GW_Amplitudes', None, projection, 
                                 'time', **kwargs)
        
    kwargs['comp'] = comp

@fig_window_open
def IMFs(self, strain:Literal['h+eq', 'hxeq', 'h+pol', 'hxpol']='h+eq',
            mode:Literal['EMD', 'EEMD']='EMD', min_imfs=0, max_imfs=10,
            spectro=False, **kwargs):
    loc = locals()
    for q in loc.keys():
        if q not in ['self', 'kwargs']:
            kwargs[q] = loc[q]
    if spectro:
        self.plotHHT(**kwargs)
    else:
        raise NotImplementedError('IMFs plotting is not implemented yet.')
        self.plotIMFs(**kwargs)
