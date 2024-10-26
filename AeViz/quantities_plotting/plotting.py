import os
import numpy as np
from AeViz.load_utils.data_load_utils import Data
from AeViz.plot_utils.plotting_utils import PlottingUtils
import matplotlib.pyplot as plt
from AeViz.grid.grid import grid
from AeViz.units import u
from AeViz.quantities_plotting.plotting_helpers import (recognize_quantity,
                                                        setup_cbars,
                                                        setup_cbars_profile,
                                                        setup_cbars_spectrogram,
                                                        setup_cbars_HHT,
                                                        normalize_indices,
                                                        get_data_to_plot,
                                                        get_plane_indices,
                                                        show_figure,
                                                        get_qt_for_label,
                                                        plot_panel,
                                                        plot_profile_panel)
from AeViz.plot_utils.utils import plot_labels, xaxis_labels, GW_limit
import cv2

class Plotting(PlottingUtils, Data):
    def __init__(self):
        PlottingUtils.__init__(self)
        Data.__init__(self)
    
    def plot1D(self, file, qt, xaxis, index1, index2, **kwargs):
        """
        Plots a line for the quantity in the xaxis. This can be either a
        radial average, angular or a single radius.
        """
        axd_letters = ['A', 'B', 'C', 'D']
        legend = None

        ## GET THE DATA
        if 'radius' not in qt and 'spheres' not in qt:
            if 'GW' in qt:
                post_data = self._Data__get_data_from_name(
                    "_".join(qt.split('_')[:-1]), file, **kwargs)
            else:
                post_data = self._Data__get_data_from_name(qt, file, **kwargs)
            if 'nu_integrated'in qt and 'all' in qt:
                legend = [r'$\nu_e$', r'$\overline{\nu}_e$', r'$\nu_x$']
            if 'nu_integrated_lum' in qt:
                post_data[:, 1:] *= 1e-53
            if 'PNS_angular_mom_all' == qt:
                legend = ['L$_x$', 'L$_y$', 'L$_z$', r'L$_\mathrm{tot}$']
            if 'kick_velocity' in qt and 'all' in qt:
                legend = ['tot',  'hydro', r'$\nu$']
        else:
            post_data = list(self._Data__get_1D_radii_data(qt, **kwargs))
            if 'all' in qt:
                legend = ['max', 'min', 'avg']
            qt = "_".join(qt.split('_')[:-1])
            for i in range(1, len(post_data)):
                post_data[i] = u.convert_to_km(post_data[i])
        
        ##POSTPROCESS THE DATA
        if xaxis != 'time':
            index1, index2 = normalize_indices(index1, index2)
            data = get_data_to_plot(index1, index2, post_data, xaxis,
                                    (self.cell.dr_integration(self.ghost), 
                                    self.cell.dtheta_integration(self.ghost),
                                    self.cell.dphi(self.ghost)))
        elif qt == 'PNS_angular_mom_all' or 'radius' in qt or \
            'kick_velocity_all' in qt:
            data = post_data[1:]
        elif 'GW' in qt:
            if self.sim_dim == 2:
                data = post_data[:, 1:]
            else:
                if qt.endswith('h+eq'):
                    data = post_data[:,1]
                elif qt.endswith('h+pol'):
                    data = post_data[:,2]
                elif qt.endswith('hxeq'):
                    data = post_data[:,3]
                elif qt.endswith('hxpol'):
                    data = post_data[:,4]
        elif type(post_data) == list or type(post_data) == tuple:
            data = post_data[1]
        elif 'radius' not in qt and 'spheres' not in qt and \
            'nu_integrated' not in qt:
            data = post_data[:, 1]
        else:
            data = post_data[:, 1:]
        
        ## CHECK THE XAXIS IS CORRECT AND GET THE GRID
        if xaxis == 'radius':
            grid = u.convert_to_km(self.cell.radius(self.ghost))
            xlabel = 'R [km]'
            scale = 'log'
        elif xaxis == 'theta':
            if self.sim_dim == 1:
                raise ValueError('Cannot plot theta in 1D.')
            grid = self.cell.theta(self.ghost)
            xlabel = r'$\theta$ [rad]'
            scale = 'linear'
        elif xaxis == 'phi':
            if self.sim_dim == 1 or self.sim_dim == 2:
                raise ValueError('Cannot plot phi in 1D or 2D.')
            grid = self.cell.phi(self.ghost)
            xlabel = '$\phi$ [rad]'
            scale = 'linear'
        elif xaxis == 'time':
            if 'radius' in qt or 'spheres' in qt :
                grid = post_data[0]
            elif type(post_data) == list or type(post_data) == tuple:
                grid = post_data[0]
            else:
                grid = post_data[:, 0]
            xlabel = 't [s]'
            scale = 'linear'
        else:
            raise ValueError('xaxis must be radius, theta, phi or time.')
        

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
            ylabel_qt = get_qt_for_label(qt, **kwargs)
            number = self.__check_axd_1D(ylabel_qt, xaxis)

            self._PlottingUtils__update_params(axd_letters[number], grid, data,
                                            None, plot_labels[ylabel_qt]['log'],
                                            None, 1, None, None, self.sim_dim,
                                            **kwargs)
            self._PlottingUtils__plot1D(axd_letters[number])
            
            ## SET THE LIMITS
            self.xlim((-0.005, grid.max()), axd_letters[number])
            if 'GW' not in qt:
                self.ylim(plot_labels[ylabel_qt]['lim'], axd_letters[number])
            else:
                self.ylim(plot_labels[ylabel_qt]['lim'](data), axd_letters[number])
            ## SET THE LABELS
            self.labels(xlabel, plot_labels[ylabel_qt]['label'],
                        axd_letters[number])
            ## SET THE SCALES
            self.Xscale(scale, axd_letters[number])
            self.Yscale(plot_labels[ylabel_qt]['log'], axd_letters[number])
            ## SET THE LEGEND        
            self.update_legend(legend, axd_letters[number])
        else:
            if 'plot' in kwargs:
                ax_letter = kwargs['plot']
            else:
                ax_letter = 'A'
            if ax_letter not in self.axd:
                ax_letter = 'A'
            self._PlottingUtils__update_params(ax_letter, grid, data,
                                            None, None, None, 1, None,
                                            None, self.sim_dim,
                                            **kwargs)
            self._PlottingUtils__plot1D(ax_letter)
 
    def plot2D(self, file, plane, qt1=None, qt2=None, qt3=None, qt4=None,
               **kwargs):
        """
        Plot a contourf plot of the quantity in the plane.
        """
        ## Set the ghost cells
        self.ghost.update_ghost_cells(t_l=3, t_r=3, p_l=3, p_r=3)
        ## Set the grid
        gr = grid(self.sim_dim, u.convert_to_km(self.cell.radius(self.ghost)),
                  self.cell.theta(self.ghost), self.cell.phi(self.ghost))
        X, Y = gr.cartesian_grid_2D(plane)
        ## Get the number of quantities to plot
        qt1, qt2, qt3, qt4 = recognize_quantity(qt1, qt2, qt3, qt4, True)
        number_of_quantities = sum(x is not None for x in [qt1, qt2, qt3, qt4])
        ## Get the indices associated to the plane
        plane, index_theta, index_phi = get_plane_indices(self, plane)
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
            plot_panel(self, 'A', file, qt1, (X, Y), (index_theta, index_phi),
                       cbars, plot_labels, plane, **kwargs)
            ## Default limits
            self.xlim((0, 100), "A")
        if qt2 is not None:
            plot_panel(self, 'B', file, qt2, (X, Y), (index_theta, index_phi),
                       cbars, plot_labels, plane, **kwargs)
        if qt3 is not None and qt4 is None:
            plot_panel(self, 'C', file, qt3, (X, Y), (index_theta, index_phi),
                       cbars, plot_labels, plane, **kwargs)
            plot_panel(self, 'D', file, qt3, (X, Y), (index_theta, index_phi),
                       cbars, plot_labels, plane, **kwargs)
        elif qt3 is not None:
            plot_panel(self, 'C', file, qt3, (X, Y), (index_theta, index_phi),
                       cbars, plot_labels, plane, **kwargs)
        if qt4 is not None:
            plot_panel(self, 'D', file, qt4, (X, Y), (index_theta, index_phi),
                       cbars, plot_labels, plane, **kwargs)
        self.ghost.restore_default()
        self.xlim(self.xlims["A"], "A")
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
            plot_profile_panel(self, 'A', qt1, cbars, plot_labels, **kwargs)
        if qt2 is not None:
            plot_profile_panel(self, 'B', qt2, cbars, plot_labels, **kwargs)
        if qt3 is not None:
            plot_profile_panel(self, 'C', qt3, cbars, plot_labels, **kwargs)
        if qt4 is not None:
            plot_profile_panel(self, 'D', qt4, cbars, plot_labels, **kwargs)
        if redo:
            for ax_letter in self.axd:
                if ax_letter.islower():
                    continue
                self._PlottingUtils__update_cbar_position(ax_letter,
                                                          cbars[ax_letter]) 
            self._PlottingUtils__redo_plot()
        show_figure()

    def plotGWDecomposition(self, qt):
        """
        Plot the GW decomposition. ONLY in 2D.
        """
        self.Close()
        self._PlotCreation__setup_axd(5, 1)
        self.Xscale('linear', 'A')
        self.Yscale('log', 'E')
        t, AE220, f_h, nuc_h, conv_h, out_h  = \
            self._Data__get_GW_decomposition_data(qt)
        if qt == 'h+eq':
            y_strain = r'$h_{+}^{eq}$ [cm]'
        elif qt == 'h+pol':
            y_strain = r'$h_{+}^{pol}$ [cm]'
        elif qt == 'hx+eq':
            y_strain = r'$h_{\times}^{eq}$ [cm]'
        elif qt == 'hx+pol':
            y_strain = r'$h_{\times}^{pol}$ [cm]'
        ylim = GW_limit(f_h)
        if self.sim_dim == 2:
            dec_label = r'$A^{E2}_{20}(r, t)$ [cm]'
        else:
            dec_label = None
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
                                            (-3, 3), 2, 
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

    def plotGWspectrogram(self, qt):
        """
        Plots the GW spectrogram and the GW signal.
        """
        redo = False
        if self.axd is not None:
            number_spect = sum([-2 in self.plot_dim[ax_letter] 
                                 for ax_letter in self.axd if ax_letter 
                                 in self.plot_dim])
            number_GW = sum([((1 in self.plot_dim[ax_letter]) and 
                             (-2 not in self.plot_dim[ax_letter]))
                             for ax_letter in self.axd if ax_letter in 
                             self.plot_dim])
            if number_spect != number_GW:
                self.Close()
                number = 1
            else:
                number = number_GW + 1
                redo = True
        else:
            number = 1
        cbars, plots = setup_cbars_spectrogram(number)
        self._PlotCreation__setup_axd(number, 5)
        post_data_GWs = self._Data__get_data_from_name(
                    'GW_Amplitudes')
        post_data_spect = self._Data__get_data_from_name(
                    'GW_spectrogram', 3.086e+22)
        
        if self.sim_dim == 2:
            dataGWs = post_data_GWs[:, 1:]
            t, f, Zxx = post_data_spect[0], post_data_spect[1], \
                post_data_spect[2]
        else:
            if qt == 'h+eq':
                dataGWs = post_data_GWs[:,1]
                t, f, Zxx = post_data_spect[0][:, 0], \
                    post_data_spect[1][..., 0], post_data_spect[2][..., 0]
            elif qt == 'h+pol':
                dataGWs = post_data_GWs[:,2]
                t, f, Zxx = post_data_spect[0][:, 1], \
                    post_data_spect[1][..., 1], post_data_spect[2][..., 1]
            elif qt == 'hxeq':
                dataGWs = post_data_GWs[:,3]
                t, f, Zxx = post_data_spect[0][:,2], \
                    post_data_spect[1][..., 2], post_data_spect[2][..., 2]
            elif qt == 'hxpol':
                dataGWs = post_data_GWs[:,4]
                t, f, Zxx = post_data_spect[0][:, 3], \
                    post_data_spect[1][..., 3], post_data_spect[2][..., 3]
        f /= 1e3
        ## 1D plot of GWs
        self._PlottingUtils__update_params(plots[0], post_data_GWs[:,0],
                                           dataGWs, None, False, None,
                                           1, None, None, self.sim_dim)
        
        self.labels(None, plot_labels['GW_Amplitudes_'+qt]['label'], plots[0])
        self.ylim(plot_labels['GW_Amplitudes_'+qt]['lim'](dataGWs), plots[0])
        self.xlim((-0.005, post_data_GWs[:,0].max()), plots[0])
        self.Xscale('linear', plots[0])
        self.Yscale('linear', plots[0])
        self._PlottingUtils__plot1D(plots[0])
        ## 2D plot of spectrogram
        self._PlottingUtils__update_params(plots[1], (t, f),
                                           Zxx, cbars[plots[1]], True, 
                                           (1e-24, 1e-20),
                                           -2, 'magma',
                r'$\frac{\mathrm{dE_{GW}}}{\mathrm{df}}$ [B$\cdot$HZ$^{-1}$]',
                                           self.sim_dim)
        self.labels('t-t$_b$ [s]', '$f$ [kHz]', plots[1])
        self._PlottingUtils__plot2Dmesh(plots[1])
        self.ylim((0, 2), plots[1])
        self.Xscale('linear', plots[1])
        self.Yscale('linear', plots[1])
        self.xlim((-0.005, post_data_GWs[:,0].max()), plots[1])
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
        t, Y, data = self._Data__get_data_from_name('data_for_barcode', **kwargs)
        if kwargs['msum']:
            self.labels('t-t$_b$ [s]', 'l', plot)
        else:
            self.labels('t-t$_b$ [s]', 'lm', plot)
        limits = max(np.abs(data.max()), np.abs(data.min())) * 0.8
        if kwargs['zero_norm']:
            cbar_label = r'$\tilde{\rho}/\tilde{\rho}_{00}$'
        self._PlottingUtils__update_params(plot, (t, Y),
                                           data, cbars[plot], False, 
                                           (-limits, limits),
                                           -3, 'seismic',
                                           cbar_label,
                                           self.sim_dim)
        self.Xscale('linear', plot)
        self.Yscale('linear', plot)
        self._PlottingUtils__plot2Dmesh(plot, **kwargs) 
        show_figure()

    def add_2Dfield(self, file, axd_letter, comp, plane):
        """
        Adds a 2D field to the figure.
        """
        if axd_letter not in self.axd:
            raise ValueError('The axis letter is not in the figure.')
        if self.plot_dim[axd_letter] != 2:
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
            self.plot_dim['C'], self.grid['C'], self.data['C'] = self.plot_dim['B'], \
                self.grid['B'], self.data['B']
            self.cbar_lv['C'], self.cbar_position['C'], self.cbar_log['C'] = \
                self.cbar_lv['B'], self.cbar_position['B'], self.cbar_log['B']
            self.cmap_color['C'], self.cbar_label['C'] = self.cmap_color['B'], \
                self.cbar_label['B']
            self.xlims['C'], self.ylims['C'] = self.xlims['B'], self.ylims['B']
            self.logX['C'], self.logY['C'] = self.logX['B'], self.logY['B']
            self.xlabels['C'], self.ylabels['C'] = self.xlabels['B'], \
                self.ylabels['B']
            
            self.plot_dim['B'], self.grid['B'], self.data['B'] = self.plot_dim['A'], \
                self.grid['A'], self.data['A']
            self.cbar_lv['B'], self.cbar_position['B'], self.cbar_log['B'] = \
                self.cbar_lv['A'], self.cbar_position['A'], self.cbar_log['A']
            self.cmap_color['B'], self.cbar_label['B'] = self.cmap_color['A'],\
                self.cbar_label['A']
            self.xlims['B'], self.ylims['B'] = self.xlims['A'], self.ylims['A']
            self.logX['B'], self.logY['B'] = self.logX['A'], self.logY['A']
            self.xlabels['B'], self.ylabels['B'] = self.xlabels['A'],\
                self.ylabels['A']

            self.plot_dim['A'], self.grid['A'], self.data['A'] = None, None, None
            self.cbar_lv['A'], self.cbar_position['A'], self.cbar_log['A'] = \
                None, None, None
            self.cmap_color['A'], self.cbar_label['A'] = None, None
            self.xlims['A'], self.ylims['A'] = None, None
            self.logX['A'], self.logY['A'] = None, None
            self.xlabels['A'], self.ylabels['A'] = None, None
            self._PlottingUtils__redo_plot()
            return number - 1
        elif self.number == 3 and self.form_factor == 2:
            if (plot_labels[qt]['label'] != self.axd['A'].get_ylabel()) or \
                        (xaxis_labels[xaxis] != self.axd['A'].get_xlabel()):
                self.Close()
                show_figure()
                self._PlotCreation__setup_axd(number, 1)
        else:
            number = 1
            for axd_letter in self.axd:
                if (plot_labels[qt]['label'] == \
                    self.axd[axd_letter].get_ylabel()) and \
                        (xaxis_labels[xaxis] == \
                            self.axd[axd_letter].get_xlabel()):
                    return number - 1
                number += 1
            if number > 4:
                raise ValueError('No more axes available.')
            self.number = number
            self._PlottingUtils__redo_plot()
        return number - 1
