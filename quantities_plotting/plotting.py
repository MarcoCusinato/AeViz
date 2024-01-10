from load_utils.data_load_utils import Data
from plot_utils.plotting_utils import PlottingUtils
from plot_utils.plot_creation import PlotCreation
import matplotlib.pyplot as plt
from scidata.grid.grid import grid


plot_labels = {'RHO': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'ENE': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'VX': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'VY': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'VZ': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'YE': {'log': False, 'lim': (0.0, 0.5), 'cmap': 'viridis', 'label': r'Y$_e$'},
        #'YZ': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'ERR': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'LRTZ': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'DENS': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'EINT': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'ENTH': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'PELE': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'TELE': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'NELE': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'PION': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'TION': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'NION': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'VELX': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'VELY': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'VELZ': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'T': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'ENTR': {'log': False, 'lim': (1.5, 25), 'cmap': 'gist_rainbow_r', 'label': r'S [k$_b$/bry]'},
        'GAMM': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'HEAT': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'DELP': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'SMOMX': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'SMOMY': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'SMOMZ': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'PGAS': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'VSOUND': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'X_n': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'X_p': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'X_alpha': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'X_h': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'Abar': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'Zbar': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'CPOT_e': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'CPOT_n': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'CPOT_p': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'CPOT_nu': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'BHEX': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'}}

class Plotting(PlottingUtils, Data):
    def __init__(self):
        super().__init__()

    def plot2D(self, qt1=None, qt2=None, qt3=None, qt4=None):
        plt.show()
        plt.ion()

        gr = grid(self.sim_dim, self.radius, self.theta, self.phi)
        X, Y = gr.cartesian_grid()
        if qt4 is not None:
            self._PlotCreation__setup_axd(4, 2)
            cbars = {'A': 'T', 'B': 'B', 'C': 'T', 'D': 'B'}
        elif qt3 is not None:
            self._PlotCreation__setup_axd(4, 2)
            cbars = {'A': 'T', 'B': 'B', 'C': 'T', 'D': 'B'}
        elif qt2 is not None:
            self._PlotCreation__setup_axd(2, 2)
            cbars = {'A': 'L', 'B': 'R'}
        elif qt1 is not None:
            self._PlotCreation__setup_axd(1, 2)
            cbars = {'A': 'R'}
        else:
            raise ValueError('No quantity given.')
        
        self._PlottingUtils__update_params('A', (X, Y), self._Data__get_data_from_name(qt1), cbars['A'], 
                                           plot_labels[qt1]['log'], plot_labels[qt1]['lim'], 2, plot_labels[qt1]['cmap'], plot_labels[qt1]['label'])  
        self._PlottingUtils__plot2D('A')
        self.labels('X [km]', 'Z [km]', 'A')
        if qt2 is not None:
            self._PlottingUtils__update_params('B', (X, Y), self._Data__get_data_from_name(qt2), cbars['B'], 
                                               plot_labels[qt2]['log'], plot_labels[qt2]['lim'], 2, plot_labels[qt2]['cmap'], plot_labels[qt2]['label'])
            self._PlottingUtils__plot2D('B')
            self.labels('X [km]', 'Z [km]', 'B')
        if qt3 is not None and qt4 is None:
            self._PlottingUtils__update_params('C', (X, Y), self._Data__get_data_from_name(qt3), cbars['C'], 
                                                  plot_labels[qt3]['log'], plot_labels[qt3]['lim'], 2, plot_labels[qt3]['cmap'], plot_labels[qt3]['label'])
            self._PlottingUtils__plot2D('C')
            self.labels('X [km]', 'Z [km]', 'C')
            self._PlottingUtils__update_params('D', (X, Y), self._Data__get_data_from_name(qt3), cbars['D'], 
                                                  plot_labels[qt3]['log'], plot_labels[qt3]['lim'], 2, plot_labels[qt3]['cmap'], plot_labels[qt3]['label'])
            self._PlottingUtils__plot2D('D')
            self.labels('X [km]', 'Z [km]', 'D')
        if qt4 is not None:
            self._PlottingUtils__update_params('D', (X, Y), self._Data__get_data_from_name(qt4), cbars['D'],
                                                  plot_labels[qt4]['log'], plot_labels[qt4]['lim'], 2, plot_labels[qt4]['cmap'], plot_labels[qt4]['label'])
            self._PlottingUtils__plot2D('D')
            self.labels('X [km]', 'Z [km]', 'D')
        self.xlim((0, 100), "A")

        
        
