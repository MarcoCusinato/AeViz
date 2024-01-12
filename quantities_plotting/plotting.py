from load_utils.data_load_utils import Data
from plot_utils.plotting_utils import PlottingUtils
from plot_utils.plot_creation import PlotCreation
import matplotlib.pyplot as plt
from scidata.grid.grid import grid
from scidata.cell.ghost import ghost


plot_labels = {'RHO': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g$\cdot$cm$^{-3}$]'},
        'ENE': {'log': True, 'lim': (1e24, 1e35), 'cmap': 'nipy_spectral', 'label': r'$E$ [erg]'},
        'VX': {'log': True, 'lim': (-1e10, 5e10), 'cmap': 'Spectral_r', 'label': r'$v_r$ [cm$\cdot$s$^{-1}$]'},
        'VY': {'log': True, 'lim': (-1e10, 1e10), 'cmap': 'Spectral_r', 'label': r'$v_\theta$ [cm$\cdot$s$^{-1}$]'},
        'VZ': {'log': True, 'lim': (1e7, 1e10), 'cmap': 'cividis', 'label': r'$v_\phi$ [cm$\cdot$s$^{-1}$]'},
        'YE': {'log': False, 'lim': (0.0, 0.5), 'cmap': 'gist_rainbow', 'label': r'Y$_e$'},
        'YN': {'log': False, 'lim': (0, 1.0), 'cmap': 'gist_rainbow', 'label': r'Y$_N$'},
        'ERR': {'log': False, 'lim': (0, 1), 'cmap': 'seismic', 'label': r'Error'},
        'LRTZ': {'log': False, 'lim': (1, 1.1), 'cmap': 'gist_rainbow', 'label': r'$\gamma$'},
        'DENS': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'EINT': {'log': True, 'lim': (1e24, 1e35), 'cmap': 'nipy_spectral', 'label': r'$E_{\mathrm{int}}$ [erg]'},
        'ENTH': {'log': True, 'lim': (1e25, 1e36), 'cmap': 'gist_stern', 'label': r'$H$ [erg'},
        'PELE': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'TELE': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'NELE': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'PION': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'TION': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'NION': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'VELX': {'log': True, 'lim': (-1e10, 5e10), 'cmap': 'Spectral_r', 'label': r'$v_r$ [cm$\cdot$s$^{-1}$]'},
        'VELY':  {'log': True, 'lim': (-1e10, 1e10), 'cmap': 'Spectral_r', 'label': r'$v_\theta$ [cm$\cdot$s$^{-1}$]'},
        'VELZ':  {'log': True, 'lim': (1e7, 1e10), 'cmap': 'cividis', 'label': r'$v_\phi$ [cm$\cdot$s$^{-1}$]'},
        'T': {'log': False, 'lim': (0, 40), 'cmap': 'inferno', 'label': r'T [MeV]'},
        'ENTR': {'log': False, 'lim': (1.5, 25), 'cmap': 'gist_rainbow_r', 'label': r'S [k$_b$/bry]'},
        'GAMM': {'log': False, 'lim': (0.5, 3.5), 'cmap': 'cividis', 'label': r'$\Gamma$'},
        'HEAT': {'log': True, 'lim': (-1e31, 1e32), 'cmap': 'Spectral_r', 'label': r'$Q_\nu$ [erg]'},
        'DELP': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r''},
        'JX': {'log': True, 'lim': (-1e24, 1e24), 'cmap': 'coolwarm', 'label': r'$j_r$ [g$\cdot$cm$^{2}\cdot$s$^{-2}$]'},
        'JY': {'log': True, 'lim': (-1e20, 1e20), 'cmap': 'Spectral_r', 'label': r'$j_\theta$ [g$\cdot$cm$^{2}\cdot$s$^{-2}$]'},
        'JZ': {'log': True, 'lim': (-1e23, 1e23), 'cmap': 'seismic', 'label': r'$j_\phi$ [g$\cdot$cm$^{2}\cdot$s$^{-2}]$'},
        'PGAS': {'log': True, 'lim': (1e25, 1e34), 'cmap': 'gist_rainbow_r', 'label': r'P$_\mathrm{gas}$ [dyn cm$^{-2}$]'},
        'VSOUND': {'log': True, 'lim': (1e8, 1e10), 'cmap': 'nipy_spectral', 'label': r'v$_\mathrm{sound}$ [cm$\cdot$s$^{-1}$]'},
        'X_n': {'log': False, 'lim': (0, 1), 'cmap': 'cividis', 'label': r'X$_n$'},
        'X_p': {'log': False, 'lim': (0, 1), 'cmap': 'viridis', 'label': r'X$_p$'},
        'X_alpha': {'log': False, 'lim': (0, 1), 'cmap': 'plasma', 'label': r'X$_\alpha$'},
        'X_h': {'log': False, 'lim': (0, 1), 'cmap': 'magma', 'label': r'X$_h$'},
        'Abar': {'log': True, 'lim': (4, 80), 'cmap': 'gist_rainbow_r', 'label': r'$\bar{A}$'},
        'Zbar': {'log': False, 'lim': (1, 34), 'cmap': 'nipy_spectral', 'label': r'$\bar{Z}$'},
        'CPOT_e': {'log': True, 'lim': (0.1, 300), 'cmap': 'coolwarm', 'label': r'$\mu_e$ [erg$\cdot$g$^{-1}$]'},
        'CPOT_n': {'log': True, 'lim': (-2e2, 3e3), 'cmap': 'bwr', 'label': r'$\mu_n$ [erg$\cdot$g$^{-1}$]'},
        'CPOT_p': {'log': True, 'lim': (-2e2, 3e3), 'cmap': 'seismic', 'label': r'$\mu_p$ [erg$\cdot$g$^{-1}$]'},
        'CPOT_nu': {'log': True, 'lim': (-2e2, 3e3), 'cmap': 'Spectral_r', 'label': r'$\mu_\nu$ [erg$\cdot$g$^{-1}$]'},
        'BHEX': {'log': True, 'lim': (1e4, 1e15), 'cmap': 'viridis', 'label': r'$\rho$ [g cm$^{-3}$]'},
        'NUEX': {'log': True, 'lim': (-1e39, 1e39), 'cmap': 'Spectral_r', 'label': r'F$_{\nu_e,r}$ [$\#_\nu\cdot$ s$^{-1}$]'},
        'NUEY': {'log': True, 'lim': (-1e38, 1e38), 'cmap': 'Spectral_r', 'label': r'F$_{\nu_e,\theta}$ [$\#_\nu\cdot$ s$^{-1}$]'},
        'NUEZ': {'log': True, 'lim': (-1e38, 1e38), 'cmap': 'Spectral_r', 'label': r'F$_{\nu_e,\phi}$ [$\#_\nu\cdot$ s$^{-1}$]'},
        'NUEE': {'log': True, 'lim': (1e25, 1e34), 'cmap': 'viridis', 'label': r'E$_{\nu_e}$ [erg$\cdot$ s$^{-1}$]'},
        'NUAX': {'log': True, 'lim': (-1e39, 1e39), 'cmap': 'Spectral_r', 'label': r'F$_{\overline{\nu}_e,r}$ [$\#_\nu\cdot$ s$^{-1}$]'},
        'NUAY': {'log': True, 'lim': (-1e38, 1e38), 'cmap': 'Spectral_r', 'label': r'F$_{\overline{\nu}_e,\theta}$ [$\#_\nu\cdot$ s$^{-1}$]'},
        'NUAZ': {'log': True, 'lim': (-1e38, 1e38), 'cmap': 'Spectral_r', 'label': r'F$_{\overline{\nu}_e,\phi}$ [$\#_\nu\cdot$ s$^{-1}$]'},
        'NUAE': {'log': True, 'lim': (1e25, 1e34), 'cmap': 'viridis', 'label': r'E$_{\overline{\nu}_e}$ [erg$\cdot$ s$^{-1}$]'},
        'NUXX': {'log': True, 'lim': (-1e39, 1e39), 'cmap': 'Spectral_r', 'label': r'F$_{\nu_x,r}$ [$\#_\nu\cdot$ s$^{-1}$]'},
        'NUXY': {'log': True, 'lim': (-1e38, 1e38), 'cmap': 'Spectral_r', 'label': r'F$_{\nu_x,\theta}$ [$\#_\nu\cdot$ s$^{-1}$]'},
        'NUXZ': {'log': True, 'lim': (-1e38, 1e38), 'cmap': 'Spectral_r', 'label': r'F$_{\nu_x,\phi}$ [$\#_\nu\cdot$ s$^{-1}$]'},
        'NUXE': {'log': True, 'lim': (1e25, 1e34), 'cmap': 'viridis', 'label': r'E$_{\nu_x}$ [erg$\cdot$ s$^{-1}$]'},
        'BX': {'log': True, 'lim': (-1e15, 1e15), 'cmap': 'coolwarm', 'label': r'B$_r$ [G]'},
        'BY': {'log': True, 'lim': (-1e15, 1e15), 'cmap': 'coolwarm', 'label': r'B$_\theta$ [G]'},
        'BZ': {'log': True, 'lim': (-1e15, 1e15), 'cmap': 'coolwarm', 'label': r'B$_\phi$ [G]'},
}

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

def get_data_to_plot(index1, index2, post_data, xaxis):
    if post_data.ndim == 1:
        data = post_data
    elif post_data.ndim == 2:
        if xaxis == 'radius':
            if type(index1) == list:
                data = []
                for i in index1:
                    data.append(post_data[i, :])
            elif index1 is None:
                data = post_data.mean(axis=0)
            else:
                data = post_data[index1, :]
        elif xaxis == 'theta':
            
            if type(index1) == list:
                data = []
                for i in index1:
                    data.append(post_data[:, i])
            elif index1 is None:
                data = post_data.mean(axis=1)
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
                data = post_data.mean(axis=(0, 1))
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
                data = post_data.mean(axis=(0, 2))
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
                data = post_data.mean(axis=(1, 2))
            else:
                data = post_data[:, index1, index2]
    return data
        


class Plotting(PlottingUtils, Data):
    def __init__(self):
        super().__init__()
        

    def  plot1D(self, qt, xaxis, index1, index2):
        axd_letters = ['A', 'B', 'C', 'D']
        number = 1 
        if not self.fig_is_open():
            plt.show()
            plt.ion()
            self._PlotCreation__setup_axd(number, 1)
        else:
            while True:
                try:
                    if plot_labels[qt]['label'] != self.axd[axd_letters[number-1]].get_ylabel():
                        number += 1
                    else:
                        break
                except:
                    if axd_letters[number-1] not in self.axd:
                        break
                if number > 4:
                    raise ValueError('No more axes available.')
            print(number)
            self.number = number
            self._PlottingUtils__redo_plot()
        
        gh = ghost(self.gh_cells)
        post_data = gh.remove_ghost_cells(self._Data__get_data_from_name(qt), self.sim_dim)
        if xaxis == 'radius':
            grid = gh.remove_ghost_cells(self.radius, self.sim_dim, 'radius')
            self.labels('R [km]', None, axd_letters[number-1])
            self.Xscale('log', axd_letters[number-1])
        elif xaxis == 'theta':
            if self.sim_dim == 1:
                raise ValueError('Cannot plot theta in 1D.')
            grid = gh.remove_ghost_cells(self.theta, self.sim_dim, 'theta')
            self.labels('$\theta$ [rad]', None, axd_letters[number-1])
        elif xaxis == 'phi':
            if self.sim_dim == 1 or self.sim_dim == 2:
                raise ValueError('Cannot plot phi in 1D or 2D.')
            grid = gh.remove_ghost_cells(self.phi, self.sim_dim, 'phi')
            self.labels('$\phi$ [rad]', None, axd_letters[number-1])
        else:
            raise ValueError('xaxis must be radius, theta or phi.')
        
        index1, index2 = normalize_indices(index1, index2)

        data = get_data_to_plot(index1, index2, post_data, xaxis)

        
        self._PlottingUtils__update_params(axd_letters[number-1], grid, data, None, 
                                           plot_labels[qt]['log'], plot_labels[qt]['lim'], 1, None, plot_labels[qt]['label']) 
        self._PlottingUtils__plot1D(axd_letters[number-1])
        

    def plot2DwithPar(self, qt1=None, qt2=None, qt3=None, qt4=None):
        plt.show()
        plt.ion()

        gr = grid(self.sim_dim, self.radius, self.theta, self.phi)
        X, Y = gr.cartesian_grid()

        qt1, qt2, qt3, qt4 = recognize_quantity(qt1, qt2, qt3, qt4, True)

        number, form_factor, cbars = setup_cbars(qt1, qt2, qt3, qt4)
        self._PlotCreation__setup_axd(number, form_factor)
        
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
            self._PlottingUtils__update_params('C', (X, Y), self._Data__get_data_from_name(qt3), cbars['C'], 
                                                  plot_labels[qt3]['log'], plot_labels[qt3]['lim'], 2, plot_labels[qt3]['cmap'], plot_labels[qt3]['label'])
            self._PlottingUtils__plot2D('C')
            self._PlottingUtils__update_params('D', (X, Y), self._Data__get_data_from_name(qt4), cbars['D'],
                                                  plot_labels[qt4]['log'], plot_labels[qt4]['lim'], 2, plot_labels[qt4]['cmap'], plot_labels[qt4]['label'])
            self._PlottingUtils__plot2D('D')
            self.labels('X [km]', 'Z [km]', 'D')
        self.xlim((0, 100), "A")

    def plot2DnoPar(self, qt1=None, type1='hydro', qt2=None, type2='hydro', 
                    qt3=None, type3='hydro', qt4=None, type4='hydro'):
        plt.show()
        plt.ion()

        gr = grid(self.sim_dim, self.radius, self.theta, self.phi)
        X, Y = gr.cartesian_grid()

        qt1, qt2, qt3, qt4 = recognize_quantity(qt1, qt2, qt3, qt4, False)

        number, form_factor, cbars = setup_cbars(qt1, qt2, qt3, qt4)
        self._PlotCreation__setup_axd(number, form_factor)
        
    def Close(self):
        self._PlotCreation__close_figure()

        
