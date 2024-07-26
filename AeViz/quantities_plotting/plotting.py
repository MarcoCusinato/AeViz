from AeViz.quantities_plotting import TERMINAL
import os
import numpy as np
from AeViz.load_utils.data_load_utils import Data
from AeViz.plot_utils.plotting_utils import PlottingUtils
import matplotlib.pyplot as plt
from AeViz.grid.grid import grid
from AeViz.units.units import units
from AeViz.utils.math_utils import function_average
from AeViz.plot_utils.utils import plot_labels, xaxis_labels, GW_limit
import cv2

u = units()

def recognize_quantity(qt1, qt2, qt3, qt4, pars):
    if qt1 == 'NUE':
        qt1, qt2, qt3, qt4 = 'NUEX', 'NUEY', 'NUEZ', 'NUEE'
    if qt1 == 'NUA':
        qt1, qt2, qt3, qt4 = 'NUAX', 'NUAY', 'NUAZ', 'NUAE'
    if qt1 == 'NUX':
        qt1, qt2, qt3, qt4 = 'NUXX', 'NUXY', 'NUXZ', 'NUXE'
    if qt1 == 'B':
        qt1, qt2, qt3, qt4 = 'BX', 'BY', 'BZ', None
    if pars:
        if qt1 == 'VEL':
            qt1, qt2, qt3, qt4 = 'VELX', 'VELY', 'VELZ', None
        if qt1 == 'V':
            qt1, qt2, qt3, qt4 = 'VX', 'VY', 'VZ', None
        if qt1 == 'X':
            qt1, qt2, qt3, qt4 = 'X_n', 'X_p', 'X_alpha', 'X_h'
        if qt1 == 'CPOT':
            qt1, qt2, qt3, qt4 = 'CPOT_e', 'CPOT_n', 'CPOT_p', 'CPOT_nu'
    return qt1, qt2, qt3, qt4

def setup_cbars(qt1, qt2, qt3, qt4):
    if qt4 is not None or qt3 is not None:
        number, form_factor = 4, 2
        cbars = {'A': 'T', 'B': 'B', 'C': 'T', 'D': 'B'}
    elif qt2 is not None:
        number, form_factor = 2, 2
        cbars = {'A': 'L', 'B': 'R'}
    elif qt1 is not None:
        number, form_factor = 1, 2
        cbars = {'A': 'R'}
    else:
        raise ValueError('No quantity given.')
    return number, form_factor, cbars

def setup_cbars_profile(qt1, qt2, qt3, qt4):
    if qt4 is not None:
        number, form_factor = 4, 4
        cbars = {'A': 'L', 'B': 'R', 'C': 'L', 'D': 'R'}
    elif qt3 is not None:
        number, form_factor = 3, 4
        cbars = {'A': 'R', 'B': 'L', 'C': 'R'}
    elif qt2 is not None:
        number, form_factor = 2, 4
        cbars = {'A': 'R', 'B': 'R'}
    elif qt1 is not None:
        number, form_factor = 1, 4
        cbars = {'A': 'R'}
    else:
        raise ValueError('No quantity given.')
    return number, form_factor, cbars

def setup_cbars_spectrogram(number):
    if number == 1:
        cbars = {'B': 'R'}
        plots = ["A", "B"]
    elif number == 2:
        cbars = {'B': 'L', 'D': 'R'}
        plots = ["C", "D"]
    elif number == 3:
        cbars = {'B': 'L', 'D': 'R', 'F': 'L'}
        plots = ["E", "F"]
    elif number == 4:
        cbars = {'B': 'L', 'D': 'R', 'F': 'L', 'H': 'R'}
        plots = ["G", "H"]
    return cbars, plots

def normalize_indices(index1, index2):
    if type(index1) == range:
        index1 = list(index1)
    if type(index2) == range:
        index2 = list(index2)
    if type(index1) == list and type(index2) == list:
        index2 = index2[0]
    return index1, index2

def get_data_to_plot(index1, index2, post_data, xaxis, dV):
    if post_data.ndim == 1:
        data = post_data
    elif post_data.ndim == 2:
        if xaxis == 'radius':
            if type(index1) == list:
                data = []
                for i in index1:
                    data.append(post_data[i, :])
            elif index1 is None:
                data = function_average(post_data, 2, 'Omega', dV[1][:, None] \
                    * dV[2])
            else:
                data = post_data[index1, :]
        elif xaxis == 'theta':
            
            if type(index1) == list:
                data = []
                for i in index1:
                    data.append(post_data[:, i])
            elif index1 is None:
                data = function_average(post_data, 2, 'theta', dV[0][None, :] \
                    * dV[2])
            else:
                data = post_data[:, index1]
    elif post_data.ndim == 3:
        if xaxis == 'radius':
            if type(index1) == list:
                data = []
                for i in index1:
                    data.append(post_data[i, index2, :])
            elif type(index2) == list:
                data = []
                for i in index2:
                    data.append(post_data[index1, i, :])
            elif index1 is None and index2 is None:
                data = function_average(post_data, 3, 'Omega', 
                                        dV[0][:, None, None] * \
                                            dV[1][None, :, None])
            else:
                data = post_data[index1, index2, :]
        elif xaxis == 'theta':
            if type(index1) == list:
                data = []
                for i in index1:
                    data.append(post_data[i, :, index2])
            elif type(index2) == list:
                data = []
                for i in index2:
                    data.append(post_data[index1, :, i])
            elif index1 is None and index2 is None:
                data = function_average(post_data, 3, 'theta',
                                        dV[0][:, None, None] \
                                            * dV[2][None, None, :])
            else:
                data = post_data[index1, :, index2]
        elif xaxis == 'phi':
            if type(index1) == list:
                data = []
                for i in index1:
                    data.append(post_data[:, i, index2])
            elif type(index2) == list:
                data = []
                for i in index2:
                    data.append(post_data[:, index1, i])
            elif index1 is None and index2 is None:
                data = function_average(post_data, 3, 'phi',
                                        dV[1][None, :, None] \
                                            * dV[2][None, None, None])
            else:
                data = post_data[:, index1, index2]
    return data
        
def show_figure():
    """
    Show the figure if the module is imported from the terminal.
    """
    if TERMINAL:
        plt.ion()
        plt.show()
        

class Plotting(PlottingUtils, Data):
    def __init__(self):
        PlottingUtils.__init__(self)
        Data.__init__(self)
    
    def plot1D(self, file, qt, xaxis, index1, index2):

        axd_letters = ['A', 'B', 'C', 'D']
        legend = None
        
        if 'radius' not in qt and 'spheres' not in qt:
            if 'GW' in qt:
                post_data = self._Data__get_data_from_name(
                    "_".join(qt.split('_')[:-1]), file)
            else:
                post_data = self._Data__get_data_from_name(qt, file)
            if 'nu_integrated'in qt and 'all' in qt:
                legend = [r'$\nu_e$', r'$\overline{\nu}_e$', r'$\nu_x$']
            if 'nu_integrated_lum' in qt:
                post_data[:, 1:] *= 1e-53
            if 'PNS_angular_mom_all' == qt:
                legend = ['L$_x$', 'L$_y$', 'L$_z$', r'L$_\mathrm{tot}$']
            if 'kick_velocity' in qt and 'all' in qt:
                legend = ['tot',  'hydro', r'$\nu$']
        else:
            post_data = list(self._Data__get_1D_radii_data(qt))
            if 'all' in qt:
                legend = ['max', 'min', 'avg']
            qt = "_".join(qt.split('_')[:-1])
            for i in range(1, len(post_data)):
                post_data[i] = u.convert_to_km(post_data[i])
        number = self.__check_axd_1D(qt, xaxis)
        if xaxis == 'radius':
            grid = u.convert_to_km(self.cell.radius(self.ghost))
            self.labels('R [km]', plot_labels[qt]['label'],
                        axd_letters[number])
            self.Xscale('log', axd_letters[number])
        elif xaxis == 'theta':
            if self.sim_dim == 1:
                raise ValueError('Cannot plot theta in 1D.')
            grid = self.cell.theta(self.ghost)
            self.labels(r'$\theta$ [rad]', plot_labels[qt]['label'],
                        axd_letters[number])
            self.Xscale('linear', axd_letters[number])
        elif xaxis == 'phi':
            if self.sim_dim == 1 or self.sim_dim == 2:
                raise ValueError('Cannot plot phi in 1D or 2D.')
            grid = self.cell.phi(self.ghost)
            self.labels('$\phi$ [rad]', plot_labels[qt]['label'],
                        axd_letters[number])
            self.Xscale('linear', axd_letters[number])
        elif xaxis == 'time':
            if 'radius' in qt or 'spheres' in qt :
                grid = post_data[0]
            elif type(post_data) == list or type(post_data) == tuple:
                grid = post_data[0]
            else:
                grid = post_data[:, 0]
            
            self.labels('t [s]', plot_labels[qt]['label'], axd_letters[number])
            self.Xscale('linear', axd_letters[number])
        else:
            raise ValueError('xaxis must be radius, theta, phi or time.')

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
        
        self._PlottingUtils__update_params(axd_letters[number], grid, data,
                                           None, plot_labels[qt]['log'], None,
                                           1, None, None, self.sim_dim)
        if 'GW' not in qt:
            self.ylim(plot_labels[qt]['lim'], axd_letters[number])
        else:
            self.ylim(plot_labels[qt]['lim'](data), axd_letters[number])
        self._PlottingUtils__plot1D(axd_letters[number])
        self.update_legend(legend, axd_letters[number])
        
        self.xlim((-0.005, grid.max()), axd_letters[number])
        self.Yscale(plot_labels[qt]['log'], axd_letters[number])
    
    def __check_axd_1D(self, qt, xaxis):
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

    def plot2D(self, file, plane, qt1=None, qt2=None, qt3=None, qt4=None):
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
        if self.sim_dim == 2:
            plane = 'xz'
            index_phi = None
            index_theta = None
        elif plane == 'xy':
            index_phi = None
            index_theta = len(self.cell.theta(self.ghost)) // 2
        elif plane == 'xz':
            index_phi = self.ghost.ghost - self.ghost.p_l
            index_theta = None
        elif plane == 'yz':
            index_theta = None
            index_phi = (len(self.cell.phi(self.ghost)) - (2 * \
                self.ghost.ghost - self.ghost.p_l - self.ghost.p_r)) // 4 + \
                    self.ghost.ghost - self.ghost.p_l
        ## Create the plot
        redo = False
        if self.axd is not None:
            redo = True
            if 'C' in self.axd and 'D' in self.axd:
                if np.all(self.data['C'] == self.data['D']):
                    del self.axd['D']
            
            if (number_of_quantities == 4) or \
                (number_of_quantities == 3 and 'B' in self.axd) or \
                (number_of_quantities == 2 and 'C' in self.axd) or \
                (number_of_quantities == 1 and 'D' in self.axd) or \
                (np.any([self.plot_dim[axd_letter] != 2 for axd_letter 
                        in [ax_l for ax_l in self.axd if 
                            not ax_l.islower()]])):
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
            ## get the data
            data = self._Data__get_data_from_name(qt1, file)
            ## get the plane cut
            data = self._Data__plane_cut(data, index_theta, index_phi)
            self._PlottingUtils__update_params('A', (X, Y), 
                                            data,
                                            cbars['A'], plot_labels[qt1]['log'],
                                            plot_labels[qt1]['lim'], 2, 
                                            plot_labels[qt1]['cmap'], 
                                            plot_labels[qt1]['label'],
                                            self.sim_dim)  
            self.labels('X [km]', 'Z [km]', 'A')
            self._PlottingUtils__plot2D('A')
            self.Xscale('linear', 'A')
            self.Yscale('linear', 'A')
        if qt2 is not None:
            data = self._Data__get_data_from_name(qt2, file)
            data = self._Data__plane_cut(data, index_theta, index_phi)
            print('here')
            self._PlottingUtils__update_params('B', (X, Y),
                                        data,
                                        cbars['B'], plot_labels[qt2]['log'],
                                        plot_labels[qt2]['lim'], 2,
                                        plot_labels[qt2]['cmap'],
                                        plot_labels[qt2]['label'],
                                        self.sim_dim)
            self.labels('X [km]', 'Z [km]', 'B')
            self._PlottingUtils__plot2D('B')
            self.Xscale('linear', 'B')
            self.Yscale('linear', 'B')
        if qt3 is not None and qt4 is None:
            data = self._Data__get_data_from_name(qt3, file)
            data = self._Data__plane_cut(data, index_theta, index_phi)
            self._PlottingUtils__update_params('C', (X, Y),
                                        data,
                                        cbars['C'], plot_labels[qt3]['log'],
                                        plot_labels[qt3]['lim'], 2, 
                                        plot_labels[qt3]['cmap'],
                                        plot_labels[qt3]['label'],
                                        self.sim_dim)
            self.labels('X [km]', 'Z [km]', 'C')
            self._PlottingUtils__plot2D('C')
            self.Xscale('linear', 'C')
            self.Yscale('linear', 'C')
            self._PlottingUtils__update_params('D', (X, Y),
                                        data,
                                        cbars['D'], plot_labels[qt3]['log'],
                                        plot_labels[qt3]['lim'], 2,
                                        plot_labels[qt3]['cmap'],
                                        plot_labels[qt3]['label'],
                                        self.sim_dim)
            self.labels('X [km]', 'Z [km]', 'D')
            self._PlottingUtils__plot2D('D')
            self.Xscale('linear', 'D')
            self.Yscale('linear', 'D')
        elif qt3 is not None:
            data = self._Data__get_data_from_name(qt3, file)
            data = self._Data__plane_cut(data, index_theta, index_phi)
            self._PlottingUtils__update_params('C', (X, Y),
                                        data,
                                        cbars['C'], plot_labels[qt3]['log'],
                                        plot_labels[qt3]['lim'], 2, 
                                        plot_labels[qt3]['cmap'],
                                        plot_labels[qt3]['label'],
                                        self.sim_dim)
            self.labels('X [km]', 'Z [km]', 'C')
            self._PlottingUtils__plot2D('C')
            self.Xscale('linear', 'C')
            self.Yscale('linear', 'C')
        if qt4 is not None:
            data = self._Data__get_data_from_name(qt4, file)
            data = self._Data__plane_cut(data, index_theta, index_phi)
            self._PlottingUtils__update_params('D', (X, Y),
                                        data,
                                        cbars['D'], plot_labels[qt4]['log'],
                                        plot_labels[qt4]['lim'], 2,
                                        plot_labels[qt4]['cmap'],
                                        plot_labels[qt4]['label'],
                                        self.sim_dim)    
            self.labels('X [km]', 'Z [km]', 'D')
            self._PlottingUtils__plot2D('D')
            self.Xscale('linear', 'D')
            self.Yscale('linear', 'D')
        self.ghost.restore_default()
        self.xlim((0, 100), "A")
        if redo:
            for ax_letter in self.axd:
                if ax_letter.islower():
                    continue
                self._PlottingUtils__update_cbar_position(ax_letter,
                                                          cbars[ax_letter]) 
            self._PlottingUtils__redo_plot()
        
        show_figure()
    
    def plotProfile(self, qt1=None, qt2=None, qt3=None, qt4=None):
        number_of_quantities = sum(x is not None for x in [qt1, qt2, qt3, qt4])
        redo = False
        
        if self.axd is not None:
            redo = True            
            if (number_of_quantities == 4) or \
                (number_of_quantities == 3 and 'B' in self.axd) or \
                (number_of_quantities == 2 and 'C' in self.axd) or \
                (number_of_quantities == 1 and 'D' in self.axd) or \
                (np.any([self.plot_dim[axd_letter] != -1 for axd_letter 
                        in [ax_l for ax_l in self.axd if ax_l not in 
                            ['a', 'b', 'c', 'd']]])):
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
        X, Y = None, None
        if qt1 is not None:
            time, _, prof = self._Data__get_profile(qt1)
            if X is None or Y is None:
                X, Y = np.meshgrid(time,
                                u.convert_to_km(self.cell.radius(self.ghost)))
            self._PlottingUtils__update_params('A', (X, Y), prof,
                                            cbars['A'], plot_labels[qt1]['log'],
                                            plot_labels[qt1]['lim'], -1, 
                                            plot_labels[qt1]['cmap'], 
                                            plot_labels[qt1]['label'],
                                            self.sim_dim)
            self.labels('t-t$_b$ [s]', 'R [km]', 'A')
            self._PlottingUtils__plot2D('A')
            self.xlim((-0.005, time.max()), 'A')
            self.Yscale('log', 'A')
            self.Xscale('linear', 'A')
        if qt2 is not None:
            time, _, prof = self._Data__get_profile(qt2)
            if X is None or Y is None:
                X, Y = np.meshgrid(time,
                                u.convert_to_km(self.cell.radius(self.ghost)))
            self._PlottingUtils__update_params('B', (X, Y), prof,
                                            cbars['B'], plot_labels[qt2]['log'],
                                            plot_labels[qt2]['lim'], -1,
                                            plot_labels[qt2]['cmap'],
                                            plot_labels[qt2]['label'],
                                            self.sim_dim)
            self.labels('t-t$_b$ [s]', 'R [km]', 'B')
            self._PlottingUtils__plot2D('B')
            self.xlim((-0.005, time.max()), 'B')
            self.Yscale('log', 'B')
            self.Xscale('linear', 'B')
        if qt3 is not None:
            time, _, prof = self._Data__get_profile(qt3)
            if X is None or Y is None:
                X, Y = np.meshgrid(time,
                                u.convert_to_km(self.cell.radius(self.ghost)))
            self._PlottingUtils__update_params('C', (X, Y), prof,
                                            cbars['C'], plot_labels[qt3]['log'],
                                            plot_labels[qt3]['lim'], -1,
                                            plot_labels[qt3]['cmap'],
                                            plot_labels[qt3]['label'],
                                            self.sim_dim)
            self.labels('t-t$_b$ [s]', 'R [km]', 'C')
            self._PlottingUtils__plot2D('C')
            self.xlim((-0.005, time.max()), 'C')
            self.Yscale('log', 'C')
            self.Xscale('linear', 'C')
        if qt4 is not None:
            time, _, prof = self._Data__get_profile(qt4)
            if X is None or Y is None:
                X, Y = np.meshgrid(time,
                                u.convert_to_km(self.cell.radius(self.ghost)))
            self._PlottingUtils__update_params('D', (X, Y), prof,
                                            cbars['D'], plot_labels[qt4]['log'],
                                            plot_labels[qt4]['lim'], -1,
                                            plot_labels[qt4]['cmap'],
                                            plot_labels[qt4]['label'],
                                            self.sim_dim)
            self.labels('t-t$_b$ [s]', 'R [km]', 'D')
            self._PlottingUtils__plot2D('D')
            self.xlim((-0.005, time.max()), 'D')
            self.Yscale('log', 'D')
            self.Xscale('linear', 'D')
        if redo:
            for ax_letter in self.axd:
                if ax_letter.islower():
                    continue
                self._PlottingUtils__update_cbar_position(ax_letter,
                                                          cbars[ax_letter]) 
            self._PlottingUtils__redo_plot()
        show_figure()

    def plotGWDecomposition(self, qt):
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
        redo = False
        if self.axd is not None:
            number_spect = sum([self.plot_dim[ax_letter] == -2 
                                 for ax_letter in self.axd if ax_letter 
                                 in self.plot_dim])
            number_GW = sum([self.plot_dim[ax_letter] == 1
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
                if ax_letter.islower() or self.plot_dim[ax_letter] == 1:
                    continue
                self._PlottingUtils__update_cbar_position(ax_letter,
                                                          cbars[ax_letter]) 
            self._PlottingUtils__redo_plot()
        show_figure()
        
    def add_2Dfield(self, file, axd_letter, comp, plane, index1):
        if axd_letter not in self.axd:
            raise ValueError('The axis letter is not in the figure.')
        if self.plot_dim[axd_letter] != 2:
            raise ValueError('The axis letter is not a 2D plot.')
        self.ghost.update_ghost_cells(t_l=3, t_r=3, p_l=3, p_r=3)
        if comp == 'velocity':
            self.__addVelocity_field(file, axd_letter, plane, index1)
        elif comp == 'Bfield':
            self.__addBfield(file, axd_letter, plane, index1)
    
    def __addVelocity_field(self, file, axd_letter, plane, index1):
        if plane == 'xy':
            vr = self._Data__get_data_from_name('radial_velocity', file)\
                [:, index1, :]
            va = self._Data__get_data_from_name('phi_velocity', file)\
                [:, index1, :]
            angle = self.cell.phi(self.ghost)
        elif plane == 'xz':
            if index1 is None or self.sim_dim == 2:
                vr = self._Data__get_data_from_name('radial_velocity', file)
                va = self._Data__get_data_from_name('theta_velocity', file)
            else:
                vr = self._Data__get_data_from_name('radial_velocity', file)\
                    [index1, :, :]
                va = self._Data__get_data_from_name('theta_velocity', file)\
                    [index1, :, :]
            angle = self.cell.theta(self.ghost)
        vx  = (vr * np.sin(angle)[:, None] +  va * np.cos(angle)[:, None])
        vy = (vr * np.cos(angle)[:, None] -  va * np.sin(angle)[:, None])
        self._PlottingUtils__update_fields_params(axd_letter,
                                                  (vx, vy), 'v')
        self._PlottingUtils__plot2Dfield(axd_letter)

    def __addBfield(self, file, axd_letter, plane, index1):
        streamlines = self.loaded_data.stream_function(file, plane)
        self._PlottingUtils__update_fields_params(axd_letter,
                                                  streamlines, 'B')
        self._PlottingUtils__plot2Dfield(axd_letter)
        
    def make_movie(self, qt1=None, qt2=None, qt3=None, qt4=None, top=None,
              plane='xz', start_time=None, end_time=None,
              vfield=False, Bfield=False, top_time=False, lims=None):
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
        self._PlotCreation__close_figure()
        self._PlottingUtils__reset_params()

        
