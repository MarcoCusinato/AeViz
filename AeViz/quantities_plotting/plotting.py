import os
import numpy as np
from AeViz.load_utils.data_load_utils import Data
from AeViz.plot_utils.plotting_utils import PlottingUtils
import matplotlib.pyplot as plt
from AeViz.units import u, aeseries, aerray
from AeViz.quantities_plotting.plotting_helpers import (recognize_quantity,
                                                        setup_cbars,
                                                        setup_cbars_profile,
                                                        setup_cbars_spectrogram,
                                                        setup_cbars_HHT,
                                                        show_figure,
                                                        plot_panel,
                                                        plot_profile_panel,
                                                        plot_hammer_panel,
                                                        remove_labelling)
from AeViz.utils.decorators.grid import _get_plane_indices
import cv2
from AeViz.plot_utils.utils import xaxis_labels

class Plotting(PlottingUtils, Data):
    def __init__(self):
        PlottingUtils.__init__(self)
        Data.__init__(self)
        self.__simple_labelling = False
        self.__no_nu = False
        
    def set_simple_labelling(self, no_nu=False):
        if self.__simple_labelling:
            self.__simple_labelling = False
        else:
            self.__simple_labelling = True
        self.__no_nu = no_nu
    
    def plot1D(self, file, qt, plane, **kwargs):
        """
        Plots a line for the quantity in the xaxis. This can be either a
        radial average, angular or a single radius.
        """
        a = kwargs.pop('a', 1.0)
        exp = kwargs.pop('exp', 1.0)
        axd_letters = ['A', 'B', 'C', 'D']
        legend = None
        d_kwargs = kwargs.copy()
        if plane != 'time':
            d_kwargs['plane'] = plane
        data = self._Data__get_data_from_name(name=qt, file=file, **d_kwargs)
        if a != 1.0 or exp != 1.0:
            if a != 1.0:
                data *= a
            if exp != 1.0:
                data = data ** exp
        # Fix the label if a != 1.0 or exp != 1.0
        ## CHECK IF ALL THE PLOTS ARE 1D
        if 'overplot' in kwargs:
            overplot = kwargs['overplot']
        else:
            overplot = False
            if self.axd is not None:
                for ax_letter in self.axd:
                    if ax_letter.islower():
                        continue
                    if any([pdim != 1 for pdim in self.plot_dim[ax_letter]]):
                        overplot = True
                        break
        if not overplot:
        ##PLOT CREATION
            if type(plane) == tuple:
                if plane[0] is None or type(plane[0]) == int:
                    if isinstance(data, list):
                        plane = data[0].return_axis_names()[0]
                    else:
                        plane = data.return_axis_names()[0]
            if self.__simple_labelling:
                data, label = remove_labelling(data, self.__no_nu)
                legend = [label]
            number = self.__check_axd_1D(data.data.label, getattr(data, plane))
            self._PlottingUtils__update_params(
                                                file=file,
                                                ax_letter=axd_letters[number],
                                                plane=plane,
                                                data=data,
                                                cbar_position=None,
                                                dim=1,
                                                sim_dim=self.sim_dim,
                                                **kwargs
                                                )
            self._PlottingUtils__plot1D(axd_letters[number])
            ## SET THE LIMITS
            if axd_letters[number] not in self.xlims:
                self.xlim(getattr(data, plane).limits, axd_letters[number])
                self.ylim(data.data.limits, axd_letters[number])
            self._PlottingUtils__save_labels(axd_letters[number])
            self.Xscale(getattr(data, plane).log, axd_letters[number])
            self.Yscale(data.data.log, axd_letters[number])
            self.update_legend(legend, axd_letters[number])
        else:
            if 'plot' in kwargs:
                ax_letter = kwargs['plot']
            else:
                # USE last active plot
                ax_letter = list(self.axd.keys())[-1]
            if ax_letter not in self.axd:
                ax_letter = list(self.axd.keys())[-1]
            self._PlottingUtils__update_params(
                                                file=file,
                                                ax_letter=ax_letter,
                                                plane=plane,
                                                data=data,
                                                cbar_position=None,
                                                dim=1,
                                                sim_dim=self.sim_dim,
                                                **kwargs
                                                )
            self._PlottingUtils__plot1D(ax_letter)

    def plot1D_2Dradius(self, qt, **kwargs):
        ax_letter = kwargs['plot']
        index = self.plot_dim[ax_letter].index(2)
        plane = self.plane[ax_letter][index]
        file = self.file[ax_letter][index]
        d_kwargs = kwargs.copy()
        d_kwargs['plane'] = plane
        d_kwargs['time'] = file
        data = self._Data__get_data_from_name(name=qt, file=None, **d_kwargs)
        self._PlottingUtils__update_params(
                                            file=file,
                                            ax_letter=ax_letter,
                                            plane='X',
                                            data=data,
                                            cbar_position=None,
                                            dim=1,
                                            sim_dim=self.sim_dim,
                                            **kwargs
                                            )
        self._PlottingUtils__plot1D(ax_letter)
 
    def plot2D(self, file, plane, qt1=None, qt2=None, qt3=None, qt4=None,
               **kwargs):
        """
        Plot a contourf plot of the quantity in the plane.
        """
        ## Set the ghost cells
        self.ghost.update_ghost_cells(t_l=3, t_r=3, p_l=3, p_r=3)
        ## Get the number of quantities to plot
        qt1, qt2, qt3, qt4 = recognize_quantity(qt1, qt2, qt3, qt4, True)
        number_of_quantities = sum(x is not None for x in [qt1, qt2, qt3, qt4])
        ## Create the plot
        redo = False
        if self.axd is not None:
            redo = True
            if 'C' in self.axd and 'D' in self.axd:
                idxc = self.plot_dim['C'].index(2)
                idxd = self.plot_dim['D'].index(2)
                if np.all(self.data['C'][idxc] == self.data['D'][idxd]):
                    del self.axd['D']
                    self._PlottingUtils__clear_param_key('D')
            
            if (number_of_quantities == 4) or \
                (number_of_quantities == 3 and 'B' in self.axd) or \
                (number_of_quantities == 2 and 'C' in self.axd) or \
                (number_of_quantities == 1 and 'D' in self.axd) or \
                (self.form_factor not in [None, 2]):
                self.Close()
                number, form_factor, cbars = setup_cbars(qt1, qt2, qt3, qt4)
            elif number_of_quantities == 1:
                if 'C' in self.axd:
                    qt4 = qt1
                    number, form_factor, cbars = setup_cbars(True, True, True,
                                                             qt4)
                elif 'B' in self.axd:
                    qt3 = qt1
                    number, form_factor, cbars = setup_cbars(True, True, qt3,
                                                             qt4)
                elif 'A' in self.axd:
                    qt2 = qt1
                    number, form_factor, cbars = setup_cbars(True, qt2, qt3,
                                                             qt4)
                qt1 = None
            elif number_of_quantities == 2:
                if 'B' in self.axd:
                    qt4 = qt2
                    qt3 = qt1
                    number, form_factor, cbars = setup_cbars(True, True, qt3,
                                                             qt4)
                elif 'A' in self.axd:
                    qt2 = qt1
                    qt3 = qt2
                    number, form_factor, cbars = setup_cbars(True, qt2, qt3,
                                                             qt4)
                qt1 = None
                qt2 = None
            elif number_of_quantities == 3:
                qt2 = qt1
                qt3 = qt2
                qt4 = qt3
                number, form_factor, cbars = setup_cbars(True, qt2, qt3, qt4)
                qt1 = None
        else:
            number, form_factor, cbars = setup_cbars(qt1, qt2, qt3, qt4)
        
        self._PlotCreation__setup_axd(number, form_factor)
        
        if qt1 is not None:
            plot_panel(self, 'A', file, qt1, cbars, plane, **kwargs)
        if qt2 is not None:
            plot_panel(self, 'B', file, qt2, cbars, plane, **kwargs)
        if qt3 is not None and qt4 is None:
            plot_panel(self, 'C', file, qt3, cbars, plane, **kwargs)
            plot_panel(self, 'D', file, qt3, cbars, plane, **kwargs)
        elif qt3 is not None:
            plot_panel(self, 'C', file, qt3, cbars, plane, **kwargs)
        if qt4 is not None:
            plot_panel(self, 'D', file, qt4, cbars, plane, **kwargs)
        self.ghost.restore_default()
        if redo:
            for ax_letter in self.axd:
                if ax_letter.islower():
                    continue
                self._PlottingUtils__update_cbar_position(ax_letter,
                                                          cbars[ax_letter]) 
            self._PlottingUtils__redo_plot()
        
        show_figure()

    def plotProfile(self, qt1=None, qt2=None, qt3=None, qt4=None, **kwargs):
        """
        Plot the time profile of the quantity.
        """
        number_of_quantities = sum(x is not None for x in [qt1, qt2, qt3, qt4])
        redo = False
        
        if self.axd is not None:
            redo = True            
            if (number_of_quantities == 4) or \
                (number_of_quantities == 3 and 'B' in self.axd) or \
                (number_of_quantities == 2 and 'C' in self.axd) or \
                (number_of_quantities == 1 and 'D' in self.axd) or \
                (self.form_factor not in [None, 4]):
                self.Close()
                number, form_factor, cbars = setup_cbars_profile(qt1, qt2,
                                                                 qt3, qt4)
            elif number_of_quantities == 1:
                if 'C' in self.axd:
                    qt4 = qt1
                    number, form_factor, cbars = setup_cbars_profile(True,
                                                            True, True, qt4)
                elif 'B' in self.axd:
                    qt3 = qt1
                    number, form_factor, cbars = setup_cbars_profile(True, True, qt3,
                                                             qt4)
                elif 'A' in self.axd:
                    qt2 = qt1
                    number, form_factor, cbars = setup_cbars_profile(True, qt2,
                                                                     qt3, qt4)
                qt1 = None
            elif number_of_quantities == 2:
                if 'B' in self.axd:
                    qt4 = qt2
                    qt3 = qt1
                    number, form_factor, cbars = setup_cbars_profile(True, 
                                                                True, qt3, qt4)
                elif 'A' in self.axd:
                    qt2 = qt1
                    qt3 = qt2
                    number, form_factor, cbars = setup_cbars_profile(True, qt2,
                                                                     qt3, qt4)
                qt1 = None
                qt2 = None
            elif number_of_quantities == 3:
                qt2 = qt1
                qt3 = qt2
                qt4 = qt3
                number, form_factor, cbars = setup_cbars_profile(True, qt2,
                                                                 qt3, qt4)
                qt1 = None
        else:
            number, form_factor, cbars = setup_cbars_profile(qt1, qt2, qt3,
                                                             qt4)
        self._PlotCreation__setup_axd(number, form_factor)
        if qt1 is not None:
            plot_profile_panel(self, 'A', qt1, cbars, **kwargs)
        if qt2 is not None:
            plot_profile_panel(self, 'B', qt2, cbars, **kwargs)
        if qt3 is not None:
            plot_profile_panel(self, 'C', qt3, cbars, **kwargs)
        if qt4 is not None:
            plot_profile_panel(self, 'D', qt4, cbars, **kwargs)
        if redo:
            for ax_letter in self.axd:
                if ax_letter.islower():
                    continue
                self._PlottingUtils__update_cbar_position(ax_letter,
                                                          cbars[ax_letter]) 
            self._PlottingUtils__redo_plot()
        show_figure()

    def plotHammer(self, file, plane, qt1=None, qt2=None, qt3=None, qt4=None, **kwargs):
        number_of_quantities = sum(x is not None for x in [qt1, qt2, qt3, qt4])
        redo = False
        
        if self.axd is not None:
            redo = True            
            if (number_of_quantities == 4) or \
                (number_of_quantities == 3 and 'B' in self.axd) or \
                (number_of_quantities == 2 and 'C' in self.axd) or \
                (number_of_quantities == 1 and 'D' in self.axd) or \
                (self.form_factor not in [None, 4]):
                self.Close()
                number, form_factor, cbars = setup_cbars_profile(qt1, qt2,
                                                                 qt3, qt4)
            elif number_of_quantities == 1:
                if 'C' in self.axd:
                    qt4 = qt1
                    number, form_factor, cbars = setup_cbars_profile(True,
                                                            True, True, qt4)
                elif 'B' in self.axd:
                    qt3 = qt1
                    number, form_factor, cbars = setup_cbars_profile(True, True, qt3,
                                                             qt4)
                elif 'A' in self.axd:
                    qt2 = qt1
                    number, form_factor, cbars = setup_cbars_profile(True, qt2,
                                                                     qt3, qt4)
                qt1 = None
            elif number_of_quantities == 2:
                if 'B' in self.axd:
                    qt4 = qt2
                    qt3 = qt1
                    number, form_factor, cbars = setup_cbars_profile(True, 
                                                                True, qt3, qt4)
                elif 'A' in self.axd:
                    qt2 = qt1
                    qt3 = qt2
                    number, form_factor, cbars = setup_cbars_profile(True, qt2,
                                                                     qt3, qt4)
                qt1 = None
                qt2 = None
            elif number_of_quantities == 3:
                qt2 = qt1
                qt3 = qt2
                qt4 = qt3
                number, form_factor, cbars = setup_cbars_profile(True, qt2,
                                                                 qt3, qt4)
                qt1 = None
        else:
            number, form_factor, cbars = setup_cbars_profile(qt1, qt2, qt3,
                                                             qt4)
        self._PlotCreation__setup_axd(number, form_factor, prj='hammer')
        if qt1 is not None:
            plot_hammer_panel(self, 'A', file, qt1, cbars, plane, **kwargs)
        if qt2 is not None:
            plot_hammer_panel(self, 'B', file, qt2, cbars, plane, **kwargs)
        if qt3 is not None:
            plot_hammer_panel(self, 'C', file, qt3, cbars, plane, **kwargs)
        if qt4 is not None:
            plot_hammer_panel(self, 'D', file, qt4, cbars, plane, **kwargs)
        if redo:
            for ax_letter in self.axd:
                if ax_letter.islower():
                    continue
                self._PlottingUtils__update_cbar_position(ax_letter,
                                                          cbars[ax_letter]) 
            self._PlottingUtils__redo_plot()
        show_figure()

    def plotGWDecomposition(self, qt, **kwargs):
        """
        Plot the GW decomposition. ONLY in 2D.
        """
        self.Close()
        self._PlotCreation__setup_axd(5, 1)
        AE220, oth  = self._Data__get_data_from_name(name=qt, file=None, **kwargs)
        f_h, nuc_h, conv_h, out_h = oth
        kwargs['color'] = 'C0'
        self._PlottingUtils__update_params(
                                            ax_letter='A',
                                            plane='time',
                                            data=f_h,
                                            cbar_position=None,
                                            dim=1,
                                            sim_dim=self.sim_dim,
                                            **kwargs
                                            )
        kwargs['color'] = 'C1'
        self._PlottingUtils__update_params(
                                            ax_letter='B',
                                            plane='time',
                                            data=nuc_h,
                                            cbar_position=None,
                                            dim=1,
                                            sim_dim=self.sim_dim,
                                            **kwargs
                                            )
        kwargs['color'] = 'C2'
        self._PlottingUtils__update_params(
                                            ax_letter='C',
                                            plane='time',
                                            data=conv_h,
                                            cbar_position=None,
                                            dim=1,
                                            sim_dim=self.sim_dim,
                                            **kwargs
                                            )
        kwargs['color'] = 'C3'
        self._PlottingUtils__update_params(
                                            ax_letter='D',
                                            plane='time',
                                            data=out_h,
                                            cbar_position=None,
                                            dim=1,
                                            sim_dim=self.sim_dim,
                                            **kwargs
                                            )
        kwargs.pop('color')
        for (axd_letter, dd) in zip(['A', 'B', 'C', 'D'],
                                       [f_h, nuc_h, conv_h, out_h]
                                       ):
            self._PlottingUtils__plot1D(axd_letter)
            self.update_legend([dd.data.label], axd_letter)
            self.ylim(dd.data.limits, axd_letter)
            self.Yscale(dd.data.log, axd_letter)
            self.Xscale(dd.time.log, axd_letter)
        self._PlottingUtils__update_params(
                                            ax_letter='E',
                                            plane=('time', 'radius'),
                                            data=AE220,
                                            cbar_position='R',
                                            dim=-1,
                                            sim_dim=self.sim_dim,
                                            **kwargs
                                            )
        self._PlottingUtils__plot2D('E')
        self.Xscale(AE220.time.log, 'E')
        self.Yscale(AE220.radius.log, 'E')
        self.xlim(AE220.time.limits, 'E')
        self.ylim(AE220.radius.limits, 'E')
        self.plot1D(None, 'innercore_radius', 'time', rad='avg', plot='E',
                    color='black', ls='dashed', lw=0.75)
        self.plot1D(None, 'PNS_nucleus_radius', 'time', rad='avg', plot='E',
                    color='black', lw=0.75)        
        show_figure()

    def plotSpectrogram(self, qt, **kwargs):
        """
        Plots the GW spectrogram and the GW signal.
        """
        redo = False
        if self.axd is not None:
            number_spect = sum([-2 in self.plot_dim[ax_letter] 
                                 for ax_letter in self.axd if ax_letter 
                                 in self.plot_dim])
            number_curve = sum([((1 in self.plot_dim[ax_letter]) and 
                             (-2 not in self.plot_dim[ax_letter]))
                             for ax_letter in self.axd if ax_letter in 
                             self.plot_dim])
            if number_spect != number_curve:
                self.Close()
                number = 1
            else:
                number = number_curve + 1
                redo = True
        else:
            number = 1
        cbars, plots = setup_cbars_spectrogram(number)
        self._PlotCreation__setup_axd(number, 5)
        data = self._Data__get_data_from_name(qt, **kwargs)
        keep_kwargs = {}
        for nm in ['window_size', 'check_spacing', 'time_range', 'scale_to',
                   'windowing', 'overlap']:
            if nm in kwargs:
                keep_kwargs[nm] = kwargs[nm]
        spectrogram = data.stft(**keep_kwargs)
    
        ## 1D plot of GWs
        self._PlottingUtils__update_params(ax_letter=plots[0],
                                           plane='time',
                                           data=data,
                                           cbar_position=None,
                                           dim=1,
                                           sim_dim=self.sim_dim,
                                           **kwargs)
        self._PlottingUtils__plot1D(plots[0])
        self.Xscale(data.time.log, plots[0])
        self.Yscale(data.data.log, plots[0])
        self.xlim(data.time.limits, plots[0])
        self.ylim(data.data.limits, plots[0])
        ## 2D plot of spectrogram
        self._PlottingUtils__update_params(ax_letter=plots[1],
                                           plane=('time', 'frequency'),
                                           data=spectrogram,
                                           cbar_position=cbars[plots[1]],
                                           dim=-2,
                                           sim_dim=self.sim_dim,
                                           **kwargs)
        self._PlottingUtils__plot2Dmesh(plots[1])
        self.Xscale(spectrogram.time.log, plots[1])
        self.Yscale(spectrogram.frequency.log, plots[1])
        self.xlim(spectrogram.time.limits, plots[1])
        self.ylim(spectrogram.frequency.limits, plots[1])
        if redo:
            for ax_letter in self.axd:
                if ax_letter.islower() or -2 not in self.plot_dim[ax_letter]:
                    continue
                self._PlottingUtils__update_cbar_position(ax_letter,
                                                          cbars[ax_letter]) 
            self._PlottingUtils__redo_plot()
        show_figure()
    
    def plot1DSpectrogram(self, qt, **kwargs):
        """
        Plots the GW spectrogram and the GW signal.
        """
        redo = False
        if self.axd is not None:
            number_spect = sum([-4 in self.plot_dim[ax_letter] 
                                 for ax_letter in self.axd if ax_letter 
                                 in self.plot_dim])
            number_curve = sum([((1 in self.plot_dim[ax_letter]) and 
                             (-4 not in self.plot_dim[ax_letter]))
                             for ax_letter in self.axd if ax_letter in 
                             self.plot_dim])
            if number_spect != number_curve:
                self.Close()
                number = 1
            else:
                number = number_curve + 1
                redo = True
        else:
            number = 1
        if number == 1:
            plots = ["A", "B"]
        elif number == 2:
            plots = ["C", "D"]
        elif number == 3:
            plots = ["E", "F"]
        elif number == 4:
            plots = ["G", "H"]
        
        self._PlotCreation__setup_axd(number, 7)
        data = self._Data__get_data_from_name(qt, **kwargs)
        keep_kwargs = {}
        for nm in ['norm', 'norm_by_max', 'time_range', 'windowing',
                   'check_spacing', 'axis']:
            if nm in kwargs:
                keep_kwargs[nm] = kwargs[nm]
        spectrogram = data.rfft(**keep_kwargs)

        ## 1D plot of quantity
        kwargs_1D = kwargs.copy()
        kwargs_1D['color'] = 'gainsboro'
        self._PlottingUtils__update_params(ax_letter=plots[0],
                                           plane='time',
                                           data=data,
                                           cbar_position=None,
                                           dim=1,
                                           sim_dim=self.sim_dim,
                                         **kwargs_1D)
        self._PlottingUtils__plot1D(plots[0])
        i_start = 0
        i_end = len(data.time)
        if 'time_range' in kwargs:
            if len(kwargs['time_range']) == 2:
                #cut the signal!!
                i_start = np.argmax(data.time >= kwargs['time_range'][0])
                i_end = np.argmax(data.time >= kwargs['time_range'][1])

            
        kwargs_1D['color'] = 'k'
        self._PlottingUtils__update_params(ax_letter=plots[0],
                                           plane='time',
                                           data=data[i_start:i_end],
                                           cbar_position=None,
                                           dim=1,
                                           sim_dim=self.sim_dim,
                                           **kwargs_1D)
        
        self._PlottingUtils__plot1D(plots[0])
        self.Xscale(data.time.log, plots[0])
        self.Yscale(data.data.log, plots[0])
        self.xlim(data.time.limits, plots[0])
        self.ylim(data.data.limits, plots[0])
        ## 1D plot of spectrogram
        self._PlottingUtils__update_params(ax_letter=plots[1],
                                           plane='frequency',
                                           data=spectrogram,
                                           cbar_position=None,
                                           dim=-4,
                                           sim_dim=self.sim_dim,
                                           **kwargs)
        self._PlottingUtils__plot1D(plots[1])
        self.Xscale(spectrogram.frequency.log, plots[1])
        self.Yscale(spectrogram.data.log, plots[1])
        self.xlim(spectrogram.frequency.limits, plots[1])
        self.ylim(spectrogram.data.limits, plots[1])
        if redo:
            self._PlottingUtils__redo_plot()
        show_figure()
        
    ## IMFs Stuff
    def plotIMFs(self, qt, **kwargs):
        ## Plots the IMFs in a top down fashion
        data = self._Data__get_data_from_name(qt, **kwargs)
        self.Close()
        self._PlotCreation__setup_axd(6, len(data))
        full = aeseries(
            aerray(np.zeros(len(data[0].data)), u.cm, 'full_h', label='GWs'),
            time=data[0].time.copy()
        )
        for dd in data:
            full += dd
        full.data.set(name='full_h', label=r'GWs')
        self._PlottingUtils__update_params(
                                           ax_letter='full',
                                           plane='time',
                                           data=full,
                                           cbar_position=None,
                                           dim=1,
                                           sim_dim=self.sim_dim,
                                           **kwargs
                                           )
        self._PlottingUtils__plot1D('full')
        for i in range(len(data)):
            self._PlottingUtils__update_params(
                                               ax_letter=f'IMF{i+1}',
                                               plane='time',
                                               data=data[i],
                                               cbar_position=None,
                                               dim=1,
                                               sim_dim=self.sim_dim,
                                               **kwargs
                                               )
            self._PlottingUtils__plot1D(f'IMF{i+1}')
        show_figure()

    def plotHHT(self, qt, **kwargs):
        """
        Plots the Hilbert-Huang spectrum of the GW signal.
        """
        redo = False
        number_spect = 0
        number = 0
        if self.axd is not None:
            number_spect = sum([-3 in self.plot_dim[ax_letter] 
                                 for ax_letter in self.axd if ax_letter 
                                 in self.plot_dim])
            if len(self.plot_dim) != number_spect:
                self.Close()
                number = 0
            else:
                number = number_spect
                redo = True
        
        if number == 4:
            self.Close()
            number = 1
        else:
            number += 1
            redo = True
        plot, cbars = setup_cbars_HHT(number)
        self._PlotCreation__setup_axd(number, 4)
        if kwargs['plot_f']:
            kwargs['IMFs'] = self._Data__get_data_from_name('IMFs', **kwargs)
            IFs = self._Data__get_data_from_name('instantaneous_frequency',
                                                **kwargs)
        spectrogram = self._Data__get_data_from_name(qt, **kwargs)
        ## 2D plot of spectrogram
        self._PlottingUtils__update_params(ax_letter=plot,
                                           plane=('time', 'frequency'),
                                           data=spectrogram,
                                           cbar_position=cbars[plot],
                                           dim=-3,
                                           sim_dim=self.sim_dim,
                                           **kwargs)
        self._PlottingUtils__plot2Dmesh(plot)
        self.Xscale(spectrogram.time.log, plot)
        self.Yscale(spectrogram.frequency.log, plot)
        self.xlim(spectrogram.time.limits, plot)
        self.ylim(spectrogram.frequency.limits, plot)
        
        ## Overplot the instantaneous frequency
        if kwargs['plot_f']:
            if not 'alpha' in kwargs:
                kwargs['alpha'] = 0.5
            if not 'color' in kwargs:
                kwargs['color'] = 'gainsboro'
            for IF in IFs:
                self._PlottingUtils__update_params(
                                            ax_letter=plot,
                                            plane='time',
                                            data=IF,
                                            cbar_position=None,
                                            dim=1,
                                            sim_dim=self.sim_dim,
                                            **kwargs
                                            )
        if redo:
            for ax_letter in self.axd:
                if ax_letter.islower():
                    continue
                self._PlottingUtils__update_cbar_position(ax_letter,
                                                          cbars[ax_letter]) 
            self._PlottingUtils__redo_plot()
        show_figure()
        
    def barcode(self, **kwargs):
        if self.fig is not None:
            self.Close()
        self._PlotCreation__setup_axd(1, 4)
        plot, cbars = setup_cbars_HHT(1)
        data = self._Data__get_data_from_name('data_for_barcode', **kwargs)
        self._PlottingUtils__update_params(
                                    ax_letter=plot,
                                    plane=('time', 'Y'),
                                    data=data,
                                    cbar_position=cbars[plot],
                                    dim=-3,
                                    sim_dim=self.sim_dim,
                                    **kwargs
                                    )
               
        self._PlottingUtils__plot2Dmesh(plot, **kwargs)
        self.Xscale(data.time.log, plot)
        self.Yscale(data.Y.log, plot)
        self.xlim(data.time.limits, plot)
        self.ylim(data.Y.limits, plot)
        show_figure()

    def add_2Dfield(self, axd_letter, comp):
        """
        Adds a 2D field to the figure.
        """
        if axd_letter not in self.axd:
            raise ValueError('The axis letter is not in the figure.')
        if 2 not in self.plot_dim[axd_letter]:
            raise ValueError('The axis letter is not a 2D plot.')
        number = self.plot_dim[axd_letter].index(2)
        file = self.file[axd_letter][number]
        plane = self.plane[axd_letter][number]
        self.ghost.update_ghost_cells(t_l=3, t_r=3, p_l=3, p_r=3)
        if comp == 'velocity':
            self.__addVelocity_field(file, axd_letter, plane, number)
        elif comp == 'Bfield':
            self.__addBfield(file, axd_letter, plane, number)
    
    def make_movie(self, qt1=None, qt2=None, qt3=None, qt4=None, top=None,
              plane='xz', start_time=None, end_time=None,
              vfield=False, Bfield=False, top_time=False, lims=None):
        """
        Makes a movie with the quantities specified by the user.
        """
        TMN = globals()["TERMINAL"]
        globals()["TERMINAL"] = False
        number_of_quantities = sum(x is not None for x in [qt1, qt2, qt3, qt4])
        if number_of_quantities != 2:
            top = None
        if start_time is not None:
            start_time_ind = self.loaded_data.find_file_from_time(start_time)
            start_time_ind = self.loaded_data.hdf_file_list.index(
                start_time_ind)
            start_time = u.convert_to_s(start_time)
        else:
            start_time_ind = None
        if end_time is not None:
            end_time_ind = self.loaded_data.find_file_from_time(end_time)
            end_time_ind = self.loaded_data.hdf_file_list.index(end_time_ind)
            end_time = u.convert_to_s(end_time)
        else:
            end_time_ind = None
        if number_of_quantities > 1:
            save_name = '_'.join(filter(None, [qt1, qt2, qt3, qt4]))
        else:
            save_name = qt1
        frame_folder = os.path.join(self.save_path, save_name)
        if not os.path.exists(frame_folder):
            os.mkdir(frame_folder)
        frame_list = []
        for (fi, f) in enumerate(
            self.loaded_data.hdf_file_list[start_time_ind:end_time_ind]):
            self.Close()
            self.plot2D(f, plane, None, qt1, qt2, qt3, qt4)
            if lims is not None:
                self.xlim(lims, 'A')
            if top is not None:
                self.plot1D(None, top, 'time', None, None)
                if top_time:
                    self.axd['A'].axvline(
                        self.loaded_data.time(f),
                        color='black', lw=0.75, ls='dashed')
                    self.xlim((start_time, end_time), 'A')
            if vfield:
                if 'C' in self.axd and 'D' in self.axd:
                    self.add_2Dfield(f, 'C', 'velocity', plane, None)
                    self.add_2Dfield(f, 'D', 'velocity', plane, None)
                elif 'C' in self.axd:
                    self.add_2Dfield(f, 'C', 'velocity', plane, None)
                elif 'B' in self.axd:
                    self.add_2Dfield(f, 'B', 'velocity', plane, None)
                else:
                    self.add_2Dfield(f, 'A', 'velocity', plane, None)
            if Bfield:
                if top is not None:
                    self.add_2Dfield(f, 'B', 'Bfield', plane, None)
                else:
                    if 'C' in self.axd and 'D' in self.axd:
                        self.add_2Dfield(f, 'A', 'Bfield', plane, None)
                        self.add_2Dfield(f, 'B', 'Bfield', plane, None)
                    else:
                        self.add_2Dfield(f, 'A', 'Bfield', plane, None)
            save_path = os.path.join(frame_folder,
                                     'frame_{:04d}.png'.format(fi))
            frame_list.append(save_path)
            self.fig.savefig(save_path, dpi=300)
        globals()["TERMINAL"] = TMN
        ## RENDER MOVIE
        frame = cv2.imread(frame_list[0])
        height, width, channels = frame.shape

        # Define the video codec and create a VideoWriter object
        fourcc = cv2.VideoWriter_fourcc(*'vp80')
        video_writer = cv2.VideoWriter(
            os.path.join(self.save_path, save_name + '.webm'), 
            fourcc, 10.0, (width, height))

        # Write each frame to the video
        for frame_name in frame_list:
            frame = cv2.imread(frame_name)
            video_writer.write(frame)

        # Release the video writer and print a success message
        video_writer.release()

    def Close(self):
        """
        Close the figure, and deletes the dictionaries.
        """
        self._PlotCreation__close_figure()
        self._PlottingUtils__reset_params()

    def __addVelocity_field(self, file, axd_letter, plane, grid_number):
        """
        Plots the 2D velocity field in the plane specified by the user.
        """
        ## Get the plane indices
        index_theta, index_phi = _get_plane_indices(self, plane)
        vr = self._Data__get_data_from_name('radial_velocity', file)
        vr = self.__plane_cut(self.sim_dim, vr, index_theta, index_phi)
        if self.sim_dim == 2:
            va = self._Data__get_data_from_name('theta_velocity', file)
            va = self.__plane_cut(self.sim_dim, va, index_theta, index_phi)
            angle = self.cell.theta(self.ghost)
            vx  = (vr * np.sin(angle)[:, None] +  va * np.cos(angle)[:, None])
            vy = (vr * np.cos(angle)[:, None] -  va * np.sin(angle)[:, None])
        else:
            theta = self.cell.theta(self.ghost)
            phi = self.cell.phi(self.ghost)
            vtheta = self._Data__get_data_from_name('theta_velocity', file)
            vtheta = self.__plane_cut(self.sim_dim, vtheta, index_theta, index_phi)
            vphi = self._Data__get_data_from_name('phi_velocity', file)
            vphi = self.__plane_cut(self.sim_dim, vphi, index_theta, index_phi)
            
            if plane == 'xy':
                theta = theta[index_theta]
                vx = vr * np.sin(theta) * np.cos(phi)[:, None] + \
                     vtheta * np.cos(theta) * np.cos(phi)[:, None] - \
                     vphi * np.sin(phi)[:, None]
                vy = vr * np.sin(theta) * np.sin(phi)[:, None] + \
                    vtheta * np.cos(theta) * np.sin(phi)[:, None] + \
                    vphi * np.cos(phi)[:, None]
            elif plane == 'xz':
                ## We drop vphi because it is moduled by sin(phi)
                ## that is zero in the xz plane
                phi = phi[index_phi]
                ## Theta angle is wrong, but gives the correct
                ## modulation of phi
                theta = np.concatenate((theta, theta + np.pi))
                vx = vr * np.sin(theta)[:, None] * np.cos(phi) + \
                     vtheta * np.cos(theta)[:, None] * np.cos(phi)
                vy = vr * np.cos(theta)[:, None] - \
                      vtheta * np.sin(theta)[:, None]
            elif plane == 'yz':
                ## We drop vphi because it is moduled by cos(phi)
                ## that is zero in the yz plane
                phi = phi[index_phi]
                ## Theta angle is wrong, but gives the correct
                ## modulation of phi
                theta = np.concatenate((theta, theta + np.pi))
                ## Minus are necessary to get the correct result
                vx = -vr * np.sin(theta)[:, None] * np.sin(phi) + \
                     -vtheta * np.cos(theta)[:, None] * np.sin(phi)
                vy = vr * np.cos(theta)[:, None] - \
                      vtheta * np.sin(theta)[:, None]
        self._PlottingUtils__update_fields_params(axd_letter, (vx, vy), 'v')
        self._PlottingUtils__plot2Dfield(axd_letter, grid_number)

    def __addBfield(self, file, axd_letter, plane, grid_number):
        """
        Plots the 2D magnetic field in the plane specified by the user.
        """
        streamlines = self.loaded_data.stream_function(file, plane.lower())
        self._PlottingUtils__update_fields_params(axd_letter,
                                                  streamlines, 'B')
        self._PlottingUtils__plot2Dfield(axd_letter, grid_number)
             
    def __check_axd_1D(self, qt, xaxis):
        """
        Helper method to check if the quantity that we want to plot is
        already in the figure. If it is, it returns the number of the
        axis where the quantity is plotted. If it is not, it resets the
        figure and returns the number of the axis where the quantity
        will be plotted.
        """
        number = 1 
        if not self.fig_is_open():
            show_figure()
            self._PlotCreation__setup_axd(number, 1)
        elif self.number == 2 and self.form_factor == 2:
            self.number = 3
            self._PlottingUtils__copy_param_key('C', 'B')
            self._PlottingUtils__copy_param_key('B', 'A')
            self._PlottingUtils__clear_param_key('A')
            self._PlottingUtils__redo_plot()
            return number - 1
        elif self.number == 3 and self.form_factor == 2:
            if (qt != self.ylabels['A']) or (xaxis.label != self.xlabels['A']):
                self.Close()
                show_figure()
                self._PlotCreation__setup_axd(number, 1)
        else:
            for axd_letter in self.axd:
                if (qt == self.ylabels[axd_letter]) and \
                    (xaxis.label == self.xlabels[axd_letter]):
                    return number - 1
                number += 1
            if number > 4:
                raise ValueError('No more axes available.')
            self.number = number
            self._PlottingUtils__redo_plot()
        return number - 1
    
    def __plane_cut(self, dim, data, indextheta=None, indexphi=None):
        if dim == 1:
            return np.tile(data, (64, 1))
        elif dim == 2:
            if (indexphi == None and indextheta == None):
                return data
            elif indextheta is not None:
                return np.tile(data[data.shape[0] // 2, :], (64, 1))
        elif dim == 3:
            if indexphi is not None:
                return np.concatenate([np.flip(data[indexphi, :, :], axis=0),
                                    data[(indexphi + data.shape[0] // 2) % 
                                            data.shape[0], :, :]], axis=0)
            elif indextheta is not None:
                
                return data[:, indextheta, :]
        else:
            raise ValueError('The simulation dimension is not supported.')
