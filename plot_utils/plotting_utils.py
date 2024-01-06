from plot_creation import PlotCreation
from load_utils.data_load_utils import Data
from matplotlib import ticker
import numpy as np

def cbar_loaction(loc):
    location = {'T': 'top', 'B': 'bottom', 'L': 'left', 'R': 'right'}
    return location[loc]

class PlottingUtils(PlotCreation):
    def __init__(self):
        self.__reset_params()

    def __update_params(self, ax_letter, grid, data, cbar_position, 
                        cbar_log, cbar_levels, dim, cmap, cbar_label):
        self.grid[ax_letter] = grid
        self.cbar_position[ax_letter] = cbar_position
        self.cbar_log[ax_letter] = cbar_log
        self.cbar_levels[ax_letter] = cbar_levels
        self.data[ax_letter] = data
        self.dim[ax_letter] = dim
        self.cmap[ax_letter] = cmap
        self.cbar_label = cbar_label
        if self.dim[ax_letter] == 2:
            self.data[ax_letter] = self.data[ax_letter].T
    
    def __reset_params(self):
        self.dim = {}
        self.grid = {}
        self.data = {}
        self.cbar_levels = {}
        self.cbar_position = {}
        self.cbar_log = {}
        self.cmap = {}
        self.cbar_label = {}

    def __plot2D(self, ax_letter):
        if self.cbar_log[ax_letter]:
            cbar_levels = np.logspace(self.cbar_levels[ax_letter][0], self.cbar_levels[ax_letter][1], 100)
            locator = ticker.LogLocator()
            fmt = lambda x, pos: '{:.0e}'.format(x)
        else:
            cbar_levels = np.linspace(self.cbar_levels[ax_letter][0], self.cbar_levels[ax_letter][1], 100)
            locator = ticker.LinearLocator()
            fmt = lambda x, pos: '{:.1f}'.format(x)
        
        pcm = self.axd[ax_letter].contourf(self.grid[ax_letter][0], self.grid[ax_letter][1], self.data[ax_letter],
                                            levels=cbar_levels, locator=locator, antialiased=True, cmap=self.cmap[ax_letter],
                                            extend='both')
       
        cbar = self.fig.colorbar(pcm, cax=self.axd[ax_letter.lower()], format=ticker.FuncFormatter(fmt), 
                                   location=cbar_loaction[self.cbar_position[ax_letter]])
        cbar.set_label(self.cbar_label[ax_letter])


