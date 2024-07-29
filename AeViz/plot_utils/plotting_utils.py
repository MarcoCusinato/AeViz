from AeViz.plot_utils.plot_creation import PlotCreation
from matplotlib import ticker
import numpy as np
from matplotlib.colors import LogNorm, SymLogNorm, Normalize


def cbar_loaction(loc):
    location = {'T': 'top', 'B': 'bottom', 'L': 'left', 'R': 'right'}
    return location[loc]

def set2Dlims(ax, xlim, ylim, number, form_factor, sim_dim):
    
    if number == 5 and form_factor == 1:
        if xlim == None:
            ax["E"].set_ylim(ylim)
        if ylim == None:
            ax["E"].set_xlim(xlim)
    else:
        if sim_dim == 2:
            set2dlims2Dsim(ax, xlim, ylim, number, form_factor)
        else:
            set2Dlims3Dsim(ax, xlim, ylim, number, form_factor)
                
def set2dlims2Dsim(ax, xlim, ylim, number, form_factor):
    if xlim == None:
        ylim = list(ylim)
        ylim.sort()
        xlim = [0]
        if number == 4 and form_factor == 2:
            ylim[0] = 0
        else:
            ylim[0] = -ylim[1]
        xlim.append(ylim[1])
        
    if ylim == None:
        xlim = list(xlim)
        xlim.sort()
        xlim[0] = 0
        if number == 4 and form_factor == 2:
            ylim = xlim
        else:
            ylim = [-xlim[1], xlim[1]]

    if number == 1 and form_factor == 2:
        ax["A"].set(xlim=(xlim[0], xlim[1]), ylim=(ylim[0], ylim[1]), aspect=1)
    elif number == 2 and form_factor == 2:
        ax["A"].set(xlim=(xlim[1], xlim[0]), ylim=(ylim[0], ylim[1]), aspect=1)
        ax["B"].set(xlim=(xlim[0], xlim[1]), ylim=(ylim[0], ylim[1]), aspect=1)
    elif number == 3 and form_factor == 2:
        ax["B"].set(xlim=(xlim[1], xlim[0]), ylim=(ylim[0], ylim[1]), aspect=1)
        ax["C"].set(xlim=(xlim[0], xlim[1]), ylim=(ylim[0], ylim[1]), aspect=1)
    elif number == 4 and form_factor == 2:
        ax["A"].set(xlim=(xlim[1], xlim[0]), ylim=(ylim[0], ylim[1]), aspect=1)
        ax["B"].set(xlim=(xlim[1], xlim[0]), ylim=(-ylim[1], ylim[0]),
                    aspect=1)
        ax["C"].set(xlim=(xlim[0], xlim[1]), ylim=(ylim[0], ylim[1]), aspect=1)
        ax["D"].set(xlim=(xlim[0], xlim[1]), ylim=(-ylim[1], ylim[0]),
                    aspect=1)

def set2Dlims3Dsim(ax, xlim, ylim, number, form_factor):
    if xlim == None:
        ylim = list(ylim)
        ylim.sort()
        ylim[0] = -ylim[1]
        xlim = ylim
    if ylim == None:
        xlim = list(xlim)
        xlim.sort()
        xlim[0] = -xlim[1]
        ylim = xlim
            
    if number == 1 and form_factor == 2:
        ax["A"].set(xlim=(xlim[0], xlim[1]), ylim=(ylim[0], ylim[1]), aspect=1)
    elif number == 2 and form_factor == 2:
        ax["A"].set(xlim=(xlim[0], 0), ylim=(ylim[0], ylim[1]), aspect=1)
        ax["B"].set(xlim=(0, xlim[1]), ylim=(ylim[0], ylim[1]), aspect=1)
    elif number == 3 and form_factor == 2:
        ax["B"].set(xlim=(xlim[0], 0), ylim=(ylim[0], ylim[1]), aspect=1)
        ax["C"].set(xlim=(0, xlim[1]), ylim=(ylim[0], ylim[1]), aspect=1)
    elif number == 4 and form_factor == 2:
        ax["A"].set(xlim=(xlim[0], 0), ylim=(0, ylim[1]), aspect=1)
        ax["B"].set(xlim=(xlim[0], 0), ylim=(ylim[0], 0), aspect=1)
        ax["C"].set(xlim=(0, xlim[1]), ylim=(0, ylim[1]), aspect=1)
        ax["D"].set(xlim=(0, xlim[1]), ylim=(ylim[0], 0), aspect=1)

class PlottingUtils(PlotCreation):
    def __init__(self):
        self.__reset_params()
        PlotCreation.__init__(self)     

    def __update_params(self, ax_letter, grid, data, cbar_position, 
                        cbar_log, cbar_levels, dim, cmap, cbar_label,
                        sim_dim):
        self.grid[ax_letter] = grid
        self.cbar_position[ax_letter] = cbar_position
        self.cbar_log[ax_letter] = cbar_log
        self.cbar_lv[ax_letter] = cbar_levels
        self.data[ax_letter] = data
        self.plot_dim[ax_letter] = dim
        self.cmap_color[ax_letter] = cmap
        self.cbar_label[ax_letter] = cbar_label
        self.sim_dimension[ax_letter] = sim_dim
         
    def __update_cbar_position(self, ax_letter, cbar_position):
        self.cbar_position[ax_letter] = cbar_position
    
    def __update_fields_params(self, ax_letter, field, field_type):
        self.field[ax_letter] = field
        self.field_type[ax_letter] = field_type
    
    def __reset_params(self):
        self.plot_dim = {}
        self.grid = {}
        self.data = {}
        self.cbar_lv = {}
        self.cbar_position = {}
        self.cbar_log = {}
        self.cmap_color = {}
        self.cbar_label = {}
        self.xlims = {}
        self.ylims = {}
        self.logX = {}
        self.logY = {}
        self.xlabels = {}
        self.ylabels = {}
        self.legend = {}
        self.field = {}
        self.field_type = {}
        self.sim_dimension = {}

    def __plot2Dmesh(self, ax_letter):
        if self.cbar_log[ax_letter]:
            if self.cbar_lv[ax_letter][0] < 0 and \
                self.cbar_lv[ax_letter][1] > 0:
                lntresh = 10 ** (np.round(
                    min(np.log10(-self.cbar_lv[ax_letter][0]),
                    np.log10(self.cbar_lv[ax_letter][1]))) - 6)
                norm = SymLogNorm(lntresh, vmin=self.cbar_lv[ax_letter][0],
                                  vmax=self.cbar_lv[ax_letter][1])
            elif self.cbar_lv[ax_letter][0] < 0 and \
                self.cbar_lv[ax_letter][1] <= 0:
                if self.cbar_lv[ax_letter][1] == 0:
                    self.cbar_lv[ax_letter] = list(self.cbar_lv[ax_letter])
                    self.cbar_lv[ax_letter][1] = -10 ** (
                        np.log10(-self.cbar_lv[ax_letter][1]) - 10)
                    self.cbar_lv[ax_letter] = tuple(self.cbar_lv[ax_letter])
                norm = SymLogNorm(-self.cbar_lv[ax_letter][1],
                                  vmin=self.cbar_lv[ax_letter][0],
                                  vmax=self.cbar_lv[ax_letter][1])
            else:
                norm = LogNorm(vmin=self.cbar_lv[ax_letter][0],
                                vmax=self.cbar_lv[ax_letter][1])
            fmt = lambda x, pos: '{:.0e}'.format(x)
        else:
            norm = Normalize(vmin=self.cbar_lv[ax_letter][0],
                             vmax=self.cbar_lv[ax_letter][1])
            fmt = lambda x, pos: '{:.1f}'.format(x)
        pcm = self.axd[ax_letter].pcolormesh(self.grid[ax_letter][0],
                                            self.grid[ax_letter][1],
                                            self.data[ax_letter], norm=norm,
                                            cmap=self.cmap_color[ax_letter],
                                            shading='gouraud')
        cbar = self.fig.colorbar(pcm, cax=self.axd[ax_letter.lower()],
                                 format=ticker.FuncFormatter(fmt),
                                 location=cbar_loaction(
                                     self.cbar_position[ax_letter]))
        cbar.set_label(self.cbar_label[ax_letter])
        
    def __plot2Dfield(self, ax_letter):
        if self.field_type[ax_letter] == 'v':
            skip = (slice(None, None, 5), slice(None, None, 5))
            self.axd[ax_letter].quiver(self.grid[ax_letter][0][skip],
                                      self.grid[ax_letter][1][skip],
                                      self.field[ax_letter][0][skip],
                                      self.field[ax_letter][1][skip],
                                      linewidths=0.01,
                                      color='black',
                                      angles='xy'
                                      )
        elif self.field_type[ax_letter] == 'B':
            self.axd[ax_letter].contour(self.grid[ax_letter][0],
                                        self.grid[ax_letter][1],
                                        self.field[ax_letter], 45,
                                        colors = 'black',
                                        linewidths=0.2)
        
    def __plot2D(self, ax_letter):
        if self.cbar_log[ax_letter]:
            if self.cbar_lv[ax_letter][0] < 0 and \
                self.cbar_lv[ax_letter][1] > 0:
                lntresh = 10 ** (np.round(
                    min(np.log10(-self.cbar_lv[ax_letter][0]),
                    np.log10(self.cbar_lv[ax_letter][1]))) - 6)
                cbar_levels = np.concatenate((
                    -np.logspace(np.log10(-self.cbar_lv[ax_letter][0]),
                                 np.log10(lntresh), 45, endpoint=False),
                    np.linspace(-lntresh, lntresh, 10, endpoint=False),
                    np.logspace(np.log10(lntresh),
                                np.log10(self.cbar_lv[ax_letter][1]), 45)))
                norm = SymLogNorm(lntresh, vmin=self.cbar_lv[ax_letter][0],
                                  vmax=self.cbar_lv[ax_letter][1])
            elif self.cbar_lv[ax_letter][0] < 0 and \
                self.cbar_lv[ax_letter][1] <= 0:
                if self.cbar_lv[ax_letter][1] == 0:
                    self.cbar_lv[ax_letter] = list(self.cbar_lv[ax_letter])
                    self.cbar_lv[ax_letter][1] = -10 ** (
                        np.log10(-self.cbar_lv[ax_letter][1]) - 10)
                    self.cbar_lv[ax_letter] = tuple(self.cbar_lv[ax_letter])
                cbar_levels = -np.logspace(np.log10(-self.cbar_lv[ax_letter][0]),
                                          np.log10(-self.cbar_lv[ax_letter][1]),
                                          100)
                norm = SymLogNorm(-self.cbar_lv[ax_letter][1],
                                  vmin=self.cbar_lv[ax_letter][0],
                                  vmax=self.cbar_lv[ax_letter][1])
            else:
                cbar_levels = np.logspace(np.log10(self.cbar_lv[ax_letter][0]),
                                          np.log10(self.cbar_lv[ax_letter][1]),
                                          100)
                norm = LogNorm(vmin=self.cbar_lv[ax_letter][0],
                                vmax=self.cbar_lv[ax_letter][1])
            fmt = lambda x, pos: '{:.0e}'.format(x)
        else:
            cbar_levels = np.linspace(self.cbar_lv[ax_letter][0],
                                      self.cbar_lv[ax_letter][1], 100)
            norm = Normalize(vmin=self.cbar_lv[ax_letter][0],
                             vmax=self.cbar_lv[ax_letter][1])
            if self.cbar_label[ax_letter] == r'Y$_e$':
                fmt = lambda x, pos: '{:.2f}'.format(x)
            else: 
                fmt = lambda x, pos: '{:.1f}'.format(x)
        pcm = self.axd[ax_letter].contourf(self.grid[ax_letter][0],
                                           self.grid[ax_letter][1],
                                           self.data[ax_letter], norm=norm,
                                           levels=cbar_levels,
                                           antialiased=True,
                                           cmap=self.cmap_color[ax_letter],
                                           extend='both')
        cbar = self.fig.colorbar(pcm, cax=self.axd[ax_letter.lower()],
                                 format=ticker.FuncFormatter(fmt), 
                                location=cbar_loaction(
                                    self.cbar_position[ax_letter]))
        cbar.set_label(self.cbar_label[ax_letter])
        if self.cbar_position[ax_letter] in ['L', 'R'] and self.plot_dim == 2:
            self.axd[ax_letter].yaxis.labelpad = -10
        if self.cbar_position[ax_letter] in ['T', 'B']:
            for lb in cbar.ax.xaxis.get_ticklabels()[::2]:
                lb.set_visible(False)

    def __plot1D(self, ax_letter):
        if type(self.data[ax_letter]) == list:
            for data in self.data[ax_letter]:
                self.axd[ax_letter].plot(self.grid[ax_letter], data)
        else:
            self.axd[ax_letter].plot(self.grid[ax_letter],
                                     self.data[ax_letter])

    def __redo_plot(self):
        self._PlotCreation__close_figure()
        self._PlotCreation__setup_axd(self.number, self.form_factor)
        for ax_letter in self.axd:
            if ax_letter not in self.plot_dim.keys():
                continue
            if self.data[ax_letter] is None:
                continue
            if self.plot_dim[ax_letter] == 2:
                self.__plot2D(ax_letter)
                set2Dlims(self.axd, self.xlims[ax_letter], None, self.number,
                          self.form_factor, self.sim_dimension[ax_letter])
            elif self.plot_dim[ax_letter] == -1:
                self.__plot2D(ax_letter)
                self.xlim(self.xlims[ax_letter], ax_letter)
                self.Yscale(self.logY[ax_letter], ax_letter)
                self.Xscale(self.logX[ax_letter], ax_letter)
            elif self.plot_dim[ax_letter] == -2:
                self.__plot2Dmesh(ax_letter)
                self.xlim(self.xlims[ax_letter], ax_letter)
                self.ylim(self.ylims[ax_letter], ax_letter)
                self.Xscale(self.logX[ax_letter], ax_letter)
                self.Yscale(self.logY[ax_letter], ax_letter)
            else:
                self.ylim(self.ylims[ax_letter], ax_letter)
                self.xlim(self.xlims[ax_letter], ax_letter)
                self.Xscale(self.logX[ax_letter], ax_letter)
                self.Yscale(self.logY[ax_letter], ax_letter)
                self.__plot1D(ax_letter)
            if ax_letter in self.legend:
                self.update_legend(self.legend[ax_letter], ax_letter)
            if ax_letter in self.field:
                self.__plot2Dfield(ax_letter)
            self.labels(self.xlabels[ax_letter], self.ylabels[ax_letter],
                        ax_letter)
        self._PlotCreation__setup_aspect()

    def xlim(self, xlim, axd_letter="A"):
        if self.plot_dim[axd_letter] == 2:
            set2Dlims(self.axd, xlim, None, self.number, self.form_factor,
                      self.sim_dimension[axd_letter])
            self.__save_lims()
        else:
            self.axd[axd_letter].set_xlim(xlim)
        self.__save_xlims(axd_letter)
        self._PlotCreation__setup_aspect()
    
    def ylim(self, ylim, axd_letter="A"):
        if self.plot_dim[axd_letter] == 2:
            set2Dlims(self.axd, None, ylim, self.number, self.form_factor,
                      self.sim_dimension[axd_letter])
            self.__save_lims()
        else:
            self.axd[axd_letter].set_ylim(ylim)
        self.__save_ylims(axd_letter)
        self._PlotCreation__setup_aspect()
    
    def Xscale(self, scale, ax_letter="A"):
        self.__save_lims(ax_letter)
        if scale == 'log':
            if self.xlims[ax_letter][0] < 0:
                lntresh = 10 ** (np.round(
                    min(np.log10(-self.ylims[ax_letter][0]),
                    np.log10(self.ylims[ax_letter][1]))) - 6)
                self.axd[ax_letter].set_xscale('symlog', linthresh=lntresh)
            else:
                self.axd[ax_letter].set_xscale('log')
        else:
            self.axd[ax_letter].set_xscale('linear')
        self.logX[ax_letter] = self.axd[ax_letter].get_xscale()
    
    def Yscale(self, scale, ax_letter="A"):
        self.__save_lims(ax_letter)
        if scale in ['log', 'symlog'] or scale == True:
            if self.ylims[ax_letter][0] < 0:
                lntresh = 10 ** (np.round(
                    min(np.log10(-self.ylims[ax_letter][0]),
                    np.log10(self.ylims[ax_letter][1]))) - 6)
                self.axd[ax_letter].set_yscale('symlog', linthresh=lntresh)
            else:
                self.axd[ax_letter].set_yscale('log')
        else:
            self.axd[ax_letter].set_yscale('linear')
        self.logY[ax_letter] = self.axd[ax_letter].get_yscale()
 
    def cmap(self, cmap, axd_letter="A"):
        self.cmap_color[axd_letter] = cmap
        self.__redo_plot()

    def cbar_levels(self, cbar_levels, axd_letter="A"):
        self.cbar_lv[axd_letter] = cbar_levels
        self.__redo_plot()
    
    def labels(self, xlabel, ylabel, axd_letter="A"):
        if xlabel is not None:
            self.axd[axd_letter].set_xlabel(xlabel)
        if ylabel is not None:
            self.axd[axd_letter].set_ylabel(ylabel)
        self.__save_labels(axd_letter)

    def __save_xlims(self, ax_letter=None):
        if ax_letter is None:
            for ax_letter in self.axd:
                self.xlims[ax_letter] = self.axd[ax_letter].get_xlim()
        else:
            self.xlims[ax_letter] = self.axd[ax_letter].get_xlim()
    
    def __save_ylims(self, ax_letter=None):
        if ax_letter is None:
            for ax_letter in self.axd:
                self.ylims[ax_letter] = self.axd[ax_letter].get_ylim()
        else:
            self.ylims[ax_letter] = self.axd[ax_letter].get_ylim()
    
    def __save_lims(self, ax_letter=None):
        if ax_letter is None:
            for ax_letter in self.axd:
                self.xlims[ax_letter] = self.axd[ax_letter].get_xlim()
                self.ylims[ax_letter] = self.axd[ax_letter].get_ylim()
        else:
            self.xlims[ax_letter] = self.axd[ax_letter].get_xlim()
            self.ylims[ax_letter] = self.axd[ax_letter].get_ylim()
 
    def __save_labels(self, ax_letter=None):
        if ax_letter is None:
            for ax_letter in self.axd:
                self.xlabels[ax_letter] = self.axd[ax_letter].get_xlabel()
                self.ylabels[ax_letter] = self.axd[ax_letter].get_ylabel()
        else:
            self.xlabels[ax_letter] = self.axd[ax_letter].get_xlabel()
            self.ylabels[ax_letter] = self.axd[ax_letter].get_ylabel()

    def __save_scale(self, ax_letter=None):
        if ax_letter is None:
            for ax_letter in self.axd:
                self.logX[ax_letter] = self.axd[ax_letter].get_xscale()
                self.logY[ax_letter] = self.axd[ax_letter].get_yscale()
        else:
            self.logX[ax_letter] = self.axd[ax_letter].get_xscale()
            self.logY[ax_letter] = self.axd[ax_letter].get_yscale()
    
    def update_legend(self, legend, axd_letter="A"):
        if legend is None:
            pass
        elif len(self.axd[axd_letter].lines) != len(legend):
            pass
        else:
            self.legend[axd_letter] = legend
            for ll in range(len(self.legend[axd_letter])):
                self.axd[axd_letter].lines[ll].set_label(
                    self.legend[axd_letter][ll])
            self.axd[axd_letter].legend(loc='upper right')
    
    def __save_params(self):
        self.__save_lims()
        self.__save_scale()
        self.__save_labels()

