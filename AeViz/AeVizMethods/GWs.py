from AeViz.AeVizMethods import *

## ---------------------------------------------------------------------
## GRAVITATIONAL WAVES
## ---------------------------------------------------------------------

"""
Module containing all the gravitational waves plotting methods.
These functions are not meant to be used standalone, but rather to be
imported into the AeViz class.
"""

@fig_window_open
def GWs(self, comp:Literal['all', 'h+eq', 'h+pol', 'hxeq', 'hxpol']='h+eq',
        decomposition=False, spectrogram=False, **kwargs):
    if self.sim_dim == 1:
        pass
    else:
        if comp == 'all' and self.sim_dim == 3:
            qt = ['h+eq', 'h+pol', 'hxeq', 'hxpol']
        elif self.sim_dim == 2:
            qt = ['h+eq']
        else:
            qt = [comp]
        if decomposition:
            for q in qt:
                self.plotGWDecomposition(q, **kwargs)
                
        elif spectrogram:
            for q in qt:
                self.plotGWspectrogram(q)
        else:
            for q in qt:
                self.plot1D(None, 'GW_Amplitudes_' + q, 'time', None, None,
                            **kwargs)

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
