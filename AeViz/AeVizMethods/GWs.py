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
    kwargs.setdefault('mode', None)
    if self.sim_dim == 2:
        kwargs['comp'] = 'h+eq'
        if decomposition:
            self.plotGWDecomposition('hydro_strain', **kwargs)
        elif kwargs['mode'] is not None and spectrogram:
            kwargs['spectrogram'] = False
            kwargs['return_letter'] = True
            kwargs['min_legend_len'] = 2
            kwargs['legend'] = [self.loaded_data.simulation_name]
            letter = AeViz_plot_panel(self, 'hchar', None, '1D', 'frequency',
                                      **kwargs)
            leg = [self.loaded_data.simulation_name]
            if 'detectors' in kwargs:
                kwargs.pop('legend', None)
                detectors = kwargs['detectors']
                if not isinstance(detectors, list):
                    detectors = [detectors]
                kwargs.pop('color', None)
                kwargs.pop('c', None)
                kwargs['overplot'] = True
                kwargs['plot'] = letter
                for det in detectors:
                    kwargs['detector'] = det
                    leg.append(det)
                    AeViz_plot_panel(self, 'ASD', None, '1D', 'frequency',
                                     **kwargs)
                if letter in self.legend or len(leg) > 1:
                    self.update_legend(leg, letter)
        else:
            AeViz_plot_panel(self, 'GW_Amplitudes', None, projection, 'time',
                             **kwargs)
    else:
        if comp == 'all':
            for cm in ['h+eq', 'h+pol', 'hxeq', 'hxpol']:
                kwargs['comp'] = cm
                if decomposition:
                    self.plotGWDecomposition('hydro_strain', **kwargs)
                elif kwargs['mode'] is not None and spectrogram:
                    kwargs['spectrogram'] = False
                    kwargs['return_letter'] = True
                    kwargs['min_legend_len'] = 2
                    kwargs['legend'] = [self.loaded_data.simulation_name]
                    letter = AeViz_plot_panel(self, 'hchar', None, '1D',
                                              'frequency', **kwargs)
                    leg = [self.loaded_data.simulation_name]
                    if 'detectors' in kwargs:
                        leg = [self.loaded_data.simulation_name]
                        detectors = kwargs['detectors']
                        if not isinstance(detectors, list):
                            detectors = [detectors]
                        kwargs.pop('color', None)
                        kwargs.pop('c', None)
                        kwargs['overplot'] = True
                        kwargs['plot'] = letter
                        for det in detectors:
                            kwargs['detector'] = det
                            leg.append(det)
                            AeViz_plot_panel(self, 'ASD', None, '1D',
                                             'frequency', **kwargs)
                        if letter in self.legend or len(leg) > 1:
                            self.update_legend(leg, letter)
                else:
                    AeViz_plot_panel(self, 'GW_Amplitudes', None, projection, 
                                     'time', **kwargs)
        else:
            kwargs['comp'] = comp
            if decomposition:
                self.plotGWDecomposition('hydro_strain', **kwargs)
            elif kwargs['mode'] is not None and spectrogram:
                kwargs['spectrogram'] = False
                kwargs['return_letter'] = True
                kwargs['min_legend_len'] = 2
                kwargs['legend'] = [self.loaded_data.simulation_name]
                letter = AeViz_plot_panel(self, 'hchar', None, '1D',
                                            'frequency', **kwargs)
                if 'detectors' in kwargs:
                    leg = [self.loaded_data.simulation_name]
                    detectors = kwargs['detectors']
                    if not isinstance(detectors, list):
                        detectors = [detectors]
                    kwargs.pop('color', None)
                    kwargs.pop('c', None)
                    kwargs['overplot'] = True
                    kwargs['plot'] = letter
                    for det in detectors:
                        kwargs['detector'] = det
                        leg.append(det)
                        AeViz_plot_panel(self, 'ASD', None, '1D', 'frequency',
                                         **kwargs)
                    if letter in self.legend or len(leg) > 1:
                        self.update_legend(leg, letter)
            else:
                AeViz_plot_panel(self, 'GW_Amplitudes', None, projection, 
                                 'time', **kwargs)
        
    kwargs['comp'] = comp

@fig_window_open
def IMFs(self, projection:Literal['1D', '2D']='1D',
         strain:Literal['h+eq', 'hxeq', 'h+pol', 'hxpol']='h+eq',
         mode:Literal['EMD', 'EEMD']='EMD', min_imfs=0, max_imfs=10,
         spectrogram=False, **kwargs):
    loc = locals()
    for q in loc.keys():
        if q not in ['self', 'kwargs']:
            kwargs[q] = loc[q]
    if 'plot_f' not in kwargs:
        kwargs['plot_f'] = False
    if projection == '2D':
        if spectrogram:
            self.plotHHT('HH_spectrum', **kwargs)
    else:
        self.plotIMFs('IMFs', **kwargs)
