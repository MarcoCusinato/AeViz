import os
import numpy as np
from AeViz.load_utils.data_load_utils import Data
from AeViz.plot_utils.plotting_utils import PlottingUtils
import matplotlib.pyplot as plt
from AeViz.units import u
from AeViz.quantities_plotting.plotting_helpers import (recognize_quantity,
                                                        setup_cbars,
                                                        setup_cbars_profile,
                                                        setup_cbars_spectrogram,
                                                        setup_cbars_HHT,
                                                        get_plane_indices,
                                                        show_figure,
                                                        plot_panel,
                                                        plot_profile_panel)
from AeViz.plot_utils.utils import GW_limit
import cv2

class Plotting(PlottingUtils, Data):
    def __init__(self):
        PlottingUtils.__init__(self)
        Data.__init__(self)
    
    def plot1D(self, file, qt, plane, **kwargs):
        """
        Plots a line for the quantity in the xaxis. This can be either a
        radial average, angular or a single radius.
        """
        axd_letters = ['A', 'B', 'C', 'D']
        legend = None
        d_kwargs = kwargs.copy()
        if plane != 'time':
            d_kwargs['plane'] = plane
        data = self._Data__get_data_from_name(name=qt, file=file, **d_kwargs)
        ## CHECK IF ALL THE PLOTS ARE 1D
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
                    plane = data.return_axis_names()[0]
            number = self.__check_axd_1D(data.data.label, getattr(data, plane))
            self._PlottingUtils__update_params(
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
                ax_letter = 'A'
            if ax_letter not in self.axd:
                ax_letter = 'A'
            self._PlottingUtils__update_params(
                                                ax_letter=ax_letter,
                                                plane=plane,
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

    def plotGWDecomposition(self, qt, **kwargs):
        """
        Plot the GW decomposition. ONLY in 2D.
        """
        self.Close()
        self._PlotCreation__setup_axd(5, 1)
        self.Xscale('linear', 'A')
        self.Yscale('log', 'E')
        t, AE220, f_h, nuc_h, conv_h, out_h  = \
            self._Data__get_GW_decomposition_data(qt, **kwargs)
        if qt == 'h+eq':
            y_strain = r'$h_{+}^{eq}$ [cm]'
        elif qt == 'h+pol':
            y_strain = r'$h_{+}^{pol}$ [cm]'
        elif qt == 'hxeq':
            y_strain = r'$h_{\times}^{eq}$ [cm]'
        elif qt == 'hxpol':
            y_strain = r'$h_{\times}^{pol}$ [cm]'
        ylim = GW_limit(f_h)
        if self.sim_dim == 2:
            dec_label = r'$A^{E2}_{20}(r, t)$ [cm]'
        else:
            dec_label = y_strain.split(' ')
            dec_label[0] = dec_label[0][:-1]
            dec_label[0] += '(r)$'
            dec_label = (' ').join(dec_label)
        if 'D' not in kwargs:
            dec_label = r'$\mathcal{D}' + dec_label[1:]
            y_strain = r'$\mathcal{D}' + y_strain[1:]
            D = 1
        else:
            D = kwargs['D']
        _, convect_radius = \
            self._Data__get_1D_radii_data('innercore_radius_avg')
        _, nuc_radius = \
            self._Data__get_1D_radii_data('PNS_nucleus_radius_avg')
        self._PlottingUtils__update_params('A', t, f_h,
                                           None, False, None,
                                           1, None, None, self.sim_dim)
        self._PlottingUtils__update_params('B', t, nuc_h,
                                           None, False, None,
                                           1, None, None, self.sim_dim)
        self._PlottingUtils__update_params('C', t, conv_h,
                                           None, False, None,
                                           1, None, None, self.sim_dim)
        self._PlottingUtils__update_params('D', t, out_h,
                                           None, False, None,
                                           1, None, None, self.sim_dim)
        color = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for (axd_letter, label, c) in zip(['A', 'B', 'C', 'D'],
                                       [r'$h_f$', r'$h_{nuc}$', r'$h_{conv}$',
                                        r'$h_{out}$'],
                                       color[:4]):
            self._PlottingUtils__plot1D(axd_letter)
            self.axd[axd_letter].get_lines()[0].set_color(c)
            self.update_legend([label], axd_letter)
            self.ylim(ylim, axd_letter)
        X, Y = np.meshgrid(t, u.convert_to_km(self.cell.radius(self.ghost)))
        self._PlottingUtils__update_params('E', (X, Y),
                                            AE220,
                                            'R', False,
                                            (-3/D, 3/D), 2, 
                                            'seismic',
                                            dec_label,
                                            self.sim_dim)
        self.axd['E'].plot(t, u.convert_to_km(convect_radius), ls='dashed',
                           color='black', lw=0.75)
        self.axd['E'].plot(t, u.convert_to_km(nuc_radius),
                           color='black', lw=0.75)
        self.fig.text(0.04, 0.685, y_strain, va="center", rotation='vertical')
        self.labels('t-t$_b$ [s]', 'R [km]', 'E')
        self._PlottingUtils__plot2D('E')
        self.xlim((-0.005, t.max()), 'A')
        
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
    
    def plotHHT(self, **kwargs):
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
        Zxx, f, t = self._Data__get_data_from_name('HH_spectrum', **kwargs)
        f /= 1e3
        ## 2D plot of spectrogram
        strain = kwargs['strain'][1] + ',' + kwargs['strain'][2:]
        label = r'$\frac{\mathrm{dE_{GW_{' + strain +  r'}}}}{\mathrm{df}}$ [B$\cdot$HZ$^{-1}$]'
        self._PlottingUtils__update_params(plot, (t, f),
                                           Zxx, cbars[plot], False, 
                                           (0, Zxx.max() * 0.45),
                                           -3, 'magma',
                                            label,
                                            self.sim_dim)
        
        self.labels('t-t$_b$ [s]', '$f$ [kHz]', plot)
        self._PlottingUtils__plot2Dmesh(plot)
        self.ylim((0, 2), plot)
        self.Xscale('linear', plot)
        self.Yscale('linear', plot)
        self.xlim((-0.005, t.max()), plot)

        ## Overplot the instantaneous frequency
        if 'inst_freq' in kwargs:
            if not 'alpha' in kwargs:
                kwargs['alpha'] = 0.5
            if not 'color' in kwargs:
                kwargs['color'] = 'gainsboro'
            data = self._Data__get_data_from_name('instantaneous_frequency',
                                                  **kwargs)
            self._PlottingUtils__update_params(plot, data[0],
                                               list(data[1] * 1e-3),
                                               None, None,
                                               None, 1, None, None,
                                               self.sim_dim, lw=.5,
                                               **kwargs)
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

    def add_2Dfield(self, file, axd_letter, comp, plane):
        """
        Adds a 2D field to the figure.
        """
        if axd_letter not in self.axd:
            raise ValueError('The axis letter is not in the figure.')
        if 2 not in self.plot_dim[axd_letter]:
            raise ValueError('The axis letter is not a 2D plot.')
        self.ghost.update_ghost_cells(t_l=3, t_r=3, p_l=3, p_r=3)
        if comp == 'velocity':
            self.__addVelocity_field(file, axd_letter, plane)
        elif comp == 'Bfield':
            self.__addBfield(file, axd_letter, plane)
    
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

    def __addVelocity_field(self, file, axd_letter, plane):
        """
        Plots the 2D velocity field in the plane specified by the user.
        """
        ## Get the plane indices
        plane, index_theta, index_phi = get_plane_indices(self, plane)
        vr = self._Data__get_data_from_name('radial_velocity', file)
        vr = self._Data__plane_cut(vr, index_theta, index_phi)
        if self.sim_dim == 2:
            va = self._Data__get_data_from_name('theta_velocity', file)
            va = self._Data__plane_cut(va, index_theta, index_phi)
            angle = self.cell.theta(self.ghost)
            vx  = (vr * np.sin(angle)[:, None] +  va * np.cos(angle)[:, None])
            vy = (vr * np.cos(angle)[:, None] -  va * np.sin(angle)[:, None])
        else:
            theta = self.cell.theta(self.ghost)
            phi = self.cell.phi(self.ghost)
            vtheta = self._Data__get_data_from_name('theta_velocity', file)
            vtheta = self._Data__plane_cut(vtheta, index_theta, index_phi)
            vphi = self._Data__get_data_from_name('phi_velocity', file)
            vphi = self._Data__plane_cut(vphi, index_theta, index_phi)
            
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
        self._PlottingUtils__plot2Dfield(axd_letter)

    def __addBfield(self, file, axd_letter, plane):
        """
        Plots the 2D magnetic field in the plane specified by the user.
        """
        streamlines = self.loaded_data.stream_function(file, plane)
        self._PlottingUtils__update_fields_params(axd_letter,
                                                  streamlines, 'B')
        self._PlottingUtils__plot2Dfield(axd_letter)
             
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
