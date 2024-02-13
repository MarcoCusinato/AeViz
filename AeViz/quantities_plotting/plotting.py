from AeViz import TERMINAL
from AeViz.load_utils.data_load_utils import Data
from AeViz.plot_utils.plotting_utils import PlottingUtils
import matplotlib.pyplot as plt
from AeViz.grid.grid import grid
from AeViz.units.units import units
from AeViz.utils.math_utils import function_average
from AeViz.plot_utils.utils import plot_labels, xaxis_labels

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
        plt.show()
        plt.ion()

class Plotting(PlottingUtils, Data):
    def __init__(self):
        PlottingUtils.__init__(self)
        Data.__init__(self)
    
    def  plot1D(self, qt, xaxis, index1, index2):

        axd_letters = ['A', 'B', 'C', 'D']
        number = self.__check_axd_1D(qt, xaxis)
        post_data = self._Data__get_data_from_name(qt)
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
        else:
            raise ValueError('xaxis must be radius, theta or phi.')
        
        index1, index2 = normalize_indices(index1, index2)
        
        data = get_data_to_plot(index1, index2, post_data, xaxis,
                                (self.cell.dr_integration(self.ghost), 
                                    self.cell.dtheta_integration(self.ghost),
                                    self.cell.dphi(self.ghost)))
        
        self._PlottingUtils__update_params(axd_letters[number], grid, data,
                                           None, plot_labels[qt]['log'], None,
                                           1, None, None) 
        
        self._PlottingUtils__plot1D(axd_letters[number])
        
        self.ylim(plot_labels[qt]['lim'], axd_letters[number])
        self.xlim((grid.min(), grid.max()), axd_letters[number])
        self.Yscale(plot_labels[qt]['log'], axd_letters[number])
    
    def __check_axd_1D(self, qt, xaxis):
        number = 1 
        if not self.fig_is_open():
            show_figure()
            self._PlotCreation__setup_axd(number, 1)
        elif self.number == 2 and self.form_factor == 2:
            self.number = 3
            self.dim['C'], self.grid['C'], self.data['C'] = self.dim['B'], \
                self.grid['B'], self.data['B']
            self.cbar_lv['C'], self.cbar_position['C'], self.cbar_log['C'] = \
                self.cbar_lv['B'], self.cbar_position['B'], self.cbar_log['B']
            self.cmap_color['C'], self.cbar_label['C'] = self.cmap_color['B'], \
                self.cbar_label['B']
            self.xlims['C'], self.ylims['C'] = self.xlims['B'], self.ylims['B']
            self.logX['C'], self.logY['C'] = self.logX['B'], self.logY['B']
            self.xlabels['C'], self.ylabels['C'] = self.xlabels['B'], \
                self.ylabels['B']
            
            self.dim['B'], self.grid['B'], self.data['B'] = self.dim['A'], \
                self.grid['A'], self.data['A']
            self.cbar_lv['B'], self.cbar_position['B'], self.cbar_log['B'] = \
                self.cbar_lv['A'], self.cbar_position['A'], self.cbar_log['A']
            self.cmap_color['B'], self.cbar_label['B'] = self.cmap_color['A'],\
                self.cbar_label['A']
            self.xlims['B'], self.ylims['B'] = self.xlims['A'], self.ylims['A']
            self.logX['B'], self.logY['B'] = self.logX['A'], self.logY['A']
            self.xlabels['B'], self.ylabels['B'] = self.xlabels['A'],\
                self.ylabels['A']

            self.dim['A'], self.grid['A'], self.data['A'] = None, None, None
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

    def plot2DwithPar(self, qt1=None, qt2=None, qt3=None, qt4=None):
        show_figure()
        self.ghost.update_ghost_cells(t_l=3, t_r=3, p_l=3, p_r=3)
        gr = grid(self.sim_dim, u.convert_to_km(self.cell.radius(self.ghost)),
                  self.cell.theta(self.ghost), self.cell.phi(self.ghost))
        X, Y = gr.cartesian_grid()

        qt1, qt2, qt3, qt4 = recognize_quantity(qt1, qt2, qt3, qt4, True)

        number, form_factor, cbars = setup_cbars(qt1, qt2, qt3, qt4)
        self._PlotCreation__setup_axd(number, form_factor)
        
        self._PlottingUtils__update_params('A', (X, Y), 
                                           self._Data__get_data_from_name(qt1),
                                           cbars['A'], plot_labels[qt1]['log'],
                                           plot_labels[qt1]['lim'], 2, 
                                           plot_labels[qt1]['cmap'], 
                                           plot_labels[qt1]['label'])  
        self._PlottingUtils__plot2D('A')
        self.labels('X [km]', 'Z [km]', 'A')
        if qt2 is not None:
            self._PlottingUtils__update_params('B', (X, Y),
                                        self._Data__get_data_from_name(qt2),
                                        cbars['B'], plot_labels[qt2]['log'],
                                        plot_labels[qt2]['lim'], 2,
                                        plot_labels[qt2]['cmap'],
                                        plot_labels[qt2]['label'])
            self._PlottingUtils__plot2D('B')
            self.labels('X [km]', 'Z [km]', 'B')
        if qt3 is not None and qt4 is None:
            self._PlottingUtils__update_params('C', (X, Y),
                                        self._Data__get_data_from_name(qt3),
                                        cbars['C'], plot_labels[qt3]['log'],
                                        plot_labels[qt3]['lim'], 2, 
                                        plot_labels[qt3]['cmap'],
                                        plot_labels[qt3]['label'])
            self._PlottingUtils__plot2D('C')
            self.labels('X [km]', 'Z [km]', 'C')
            self._PlottingUtils__update_params('D', (X, Y),
                                        self._Data__get_data_from_name(qt3),
                                        cbars['D'], plot_labels[qt3]['log'],
                                        plot_labels[qt3]['lim'], 2,
                                        plot_labels[qt3]['cmap'],
                                        plot_labels[qt3]['label'])
            self._PlottingUtils__plot2D('D')
            self.labels('X [km]', 'Z [km]', 'D')
        if qt4 is not None:
            self._PlottingUtils__update_params('C', (X, Y),
                                        self._Data__get_data_from_name(qt3),
                                        cbars['C'], plot_labels[qt3]['log'],
                                        plot_labels[qt3]['lim'], 2,
                                        plot_labels[qt3]['cmap'],
                                        plot_labels[qt3]['label'])
            self._PlottingUtils__plot2D('C')
            self._PlottingUtils__update_params('D', (X, Y),
                                        self._Data__get_data_from_name(qt4),
                                        cbars['D'], plot_labels[qt4]['log'],
                                        plot_labels[qt4]['lim'], 2,
                                        plot_labels[qt4]['cmap'],
                                        plot_labels[qt4]['label'])
            self._PlottingUtils__plot2D('D')
            self.labels('X [km]', 'Z [km]', 'D')
        self.ghost.restore_default()
        self.xlim((0, 100), "A")

        for ax_letter in self.axd:
            if ax_letter in ['a', 'b', 'c', 'd']:
                continue
            self.Xscale('linear', ax_letter)
            self.Yscale('linear', ax_letter)


    def plot2DnoPar(self, qt1=None, type1='hydro', qt2=None, type2='hydro', 
                    qt3=None, type3='hydro', qt4=None, type4='hydro'):
        show_figure()

        gr = grid(self.sim_dim, self.radius, self.theta, self.phi)
        X, Y = gr.cartesian_grid()

        qt1, qt2, qt3, qt4 = recognize_quantity(qt1, qt2, qt3, qt4, False)

        number, form_factor, cbars = setup_cbars(qt1, qt2, qt3, qt4)
        self._PlotCreation__setup_axd(number, form_factor)
        
    def Close(self):
        self._PlotCreation__close_figure()

        
