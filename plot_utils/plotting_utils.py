from plot_creation import PlotCreation
from load_utils.data_load_utils import Data
from matplotlib import ticker
import numpy as np

def cbar_loaction(loc):
    location = {'T': 'top', 'B': 'bottom', 'L': 'left', 'R': 'right'}
    return location[loc]

def set2Dlims(ax, xlim, ylim, number, form_factor):
    if xlim==None:
        ylim = list(ylim)
        ylim.sort()
        xlim = [0]
        if number == 4 and form_factor == 1:
            ylim[0] = 0
        else:
            ylim[0] = -ylim[1]
        xlim.append(ylim[1])
        
    if ylim==None:
        xlim = list(xlim)
        xlim.sort()
        xlim[0] = 0
        if number == 4 and form_factor == 1:
            ylim = xlim
        else:
            ylim = [-xlim[1], xlim[1]]

    if number == 1 and form_factor == 2:
        ax["A"].set(xlim=(xlim[0], xlim[1]), ylim=(ylim[0], ylim[1]), aspect=1)
    elif number == 2 and form_factor == 2:
        ax["A"].set(xlim=(xlim[1], xlim[0]), ylim=(ylim[0], ylim[1]), aspect=1)
        ax["B"].set(xlim=(xlim[0], xlim[1]), ylim=(ylim[0], ylim[1]), aspect=1)
    elif number == 3 and form_factor == 1:
        ax["B"].set(xlim=(xlim[1], xlim[0]), ylim=(ylim[0], ylim[1]), aspect=1)
        ax["C"].set(xlim=(xlim[0], xlim[1]), ylim=(ylim[0], ylim[1]), aspect=1)
    elif number == 4 and form_factor == 1:
        ax["A"].set(xlim=(xlim[1], xlim[0]), ylim=(ylim[0], ylim[1]), aspect=1)
        ax["B"].set(xlim=(xlim[0], xlim[1]), ylim=(ylim[0], ylim[1]), aspect=1)
        ax["C"].set(xlim=(xlim[1], xlim[0]), ylim=(-ylim[1], ylim[0]), aspect=1)
        ax["D"].set(xlim=(xlim[0], xlim[1]), ylim=(-ylim[1], ylim[0]), aspect=1)


class PlottingUtils(PlotCreation):
    def __init__(self):
        self.__reset_params()

    def __update_params(self, ax_letter, grid, data, cbar_position, 
                        cbar_log, cbar_levels, dim, cmap, cbar_label):
        self.grid[ax_letter] = grid
        self.cbar_position[ax_letter] = cbar_position
        self.cbar_log[ax_letter] = cbar_log
        self.cbar_lv[ax_letter] = cbar_levels
        self.data[ax_letter] = data
        self.dim[ax_letter] = dim
        self.cmap_color[ax_letter] = cmap
        self.cbar_label = cbar_label
        if self.dim[ax_letter] == 2:
            self.data[ax_letter] = self.data[ax_letter].T
    
    def __reset_params(self):
        self.dim = {}
        self.grid = {}
        self.data = {}
        self.cbar_lv = {}
        self.cbar_position = {}
        self.cbar_log = {}
        self.cmap_color = {}
        self.cbar_label = {}
        self.xlims = {}
        self.ylims = {}
        self.xlabels = {}
        self.ylabels = {}

    def __plot2D(self, ax_letter):
        if self.cbar_log[ax_letter]:
            cbar_levels = np.logspace(self.cbar_lv[ax_letter][0], self.cbar_lv[ax_letter][1], 100)
            locator = ticker.LogLocator()
            fmt = lambda x, pos: '{:.0e}'.format(x)
        else:
            cbar_levels = np.linspace(self.cbar_lv[ax_letter][0], self.cbar_lv[ax_letter][1], 100)
            locator = ticker.LinearLocator()
            fmt = lambda x, pos: '{:.1f}'.format(x)
        
        pcm = self.axd[ax_letter].contourf(self.grid[ax_letter][0], self.grid[ax_letter][1], self.data[ax_letter],
                                            levels=cbar_levels, locator=locator, antialiased=True, cmap=self.cmap_color[ax_letter],
                                            extend='both')
       
        cbar = self.fig.colorbar(pcm, cax=self.axd[ax_letter.lower()], format=ticker.FuncFormatter(fmt), 
                                   location=cbar_loaction[self.cbar_position[ax_letter]])
        cbar.set_label(self.cbar_label[ax_letter])

    def __plot1D(self, ax_letter):
        if type(self.data[ax_letter]) == list:
            for data in self.data[ax_letter]:
                self.axd[ax_letter].plot(self.grid[ax_letter], data)
        else:
            self.axd[ax_letter].plot(self.grid[ax_letter], self.data[ax_letter])

    def __redo_plot(self):
        self.__close_figure()
        self.__setup_axd(self.number, self.form_factor)
        for ax_letter in self.axd:
            self.xlim(self.xlims[ax_letter], ax_letter)
            self.ylim(self.ylims[ax_letter], ax_letter)
            self.labels(self.xlabels[ax_letter], self.ylabels[ax_letter], ax_letter)
            if self.dim[ax_letter] == 2:
                self.__plot2D(ax_letter)
            else:
                self.__plot1D(ax_letter)
        self.__setup_aspect()

    def xlim(self, xlim, axd_letter="A"):
        if self.dim[axd_letter] == 2:
            set2Dlims(self.axd, xlim, None, self.number, self.form_factor)
        else:
            self.axd[axd_letter].set_xlim(xlim)
        self.__save_lims()
    
    def ylim(self, ylim, axd_letter="A"):
        if self.dim[axd_letter] == 2:
            set2Dlims(self.axd, None, ylim, self.number, self.form_factor)
        else:
            self.axd[axd_letter].set_ylim(ylim)
        self.__save_lims()
    
    def cmap(self, cmap, axd_letter="A"):
        self.cmap_color[axd_letter] = cmap
        self.__redo_plot()

    def cbar_levels(self, cbar_levels, axd_letter="A"):
        self.cbar_lv[axd_letter] = cbar_levels
        self.__redo_plot()
    
    def __save_lims(self):
        for ax_letter in self.axd:
            self.xlims[ax_letter] = self.axd[ax_letter].get_xlim()
            self.ylims[ax_letter] = self.axd[ax_letter].get_ylim()
    
    def __save_labels(self):
        for ax_letter in self.axd:
            self.xlabels[ax_letter] = self.axd[ax_letter].get_xlabel()
            self.ylabels[ax_letter] = self.axd[ax_letter].get_ylabel()


