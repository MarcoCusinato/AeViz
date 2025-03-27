from AeViz.plot_utils.plot_creation import PlotCreation
from matplotlib import ticker
import numpy as np
from matplotlib.colors import LogNorm, SymLogNorm, Normalize
from AeViz.plot_utils.limits_utils import set2Dlims
from AeViz.plot_utils.figure_utils import cbar_loaction


class PlottingUtils(PlotCreation):
    """
    Class with utilities for plotting the data and manipulating the
    essential plot parametes, eg limits, scales, labels, etc.
    """
    def __init__(self):
        """
        No parameters are needed to initialize the class.
        Once we initialize it we call PlotCreation class to have the
        figure and axes ready. Also we initialize the dictionaries
        containing the plot parameters.
        """
        self.__reset_params()
        PlotCreation.__init__(self)     

    def xlim(self, xlim, axd_letter="A"):
        """
        Sets the x limits of the axes corresponding to the axd_letter
        plot. In case of 2D plots, it sets the limits of all the axes.
        The limits are saved in the xlims dictionary.
        """
        if 2 in self.plot_dim[axd_letter]:
            for idx in range(len(self.plot_dim[axd_letter])):
                if self.plot_dim[axd_letter][idx] == 2:
                    sdim = self.sim_dimension[axd_letter][idx]
                    grid = self.grid[axd_letter][idx][0]
                    if np.isclose(np.abs(grid.max()), np.abs(grid.min())):
                        plane = "xy"
                    else:
                        plane = "xz"
                    break
            set2Dlims(self.axd, xlim, None, self.number, self.form_factor,
                      sdim, plane)
            self.__save_lims()
        else:
            self.axd[axd_letter].set_xlim(xlim)
        self.__save_xlims(axd_letter)
        self._PlotCreation__setup_aspect()
    
    def ylim(self, ylim, axd_letter="A"):
        """
        Sets the y limits of the axes corresponding to the axd_letter
        plot. In case of 2D plots, it sets the limits of all the axes.
        The limits are saved in the ylims dictionary.
        """
        if 2 in self.plot_dim[axd_letter]:
            for idx in range(len(self.plot_dim[axd_letter])):
                if self.plot_dim[axd_letter][idx] == 2:
                    sdim = self.sim_dimension[axd_letter][idx]
                    grid = self.grid[axd_letter][idx]
                    if np.isclose(np.abs(grid.max()), np.abs(grid.min())):
                        plane = "xy"
                    else:
                        plane = "xz"
                    break
            set2Dlims(self.axd, None, ylim, self.number, self.form_factor,
                      sdim, plane)
            self.__save_lims()
        else:
            self.axd[axd_letter].set_ylim(ylim)
        self.__save_ylims(axd_letter)
        self._PlotCreation__setup_aspect()
    
    def Xscale(self, scale, ax_letter="A"):
        """
        Sets the x axis scale of the plot at the cooresponding letter.
        If the scale is 'log' and the limits are negative, we use a 
        custom symlog scale. 
        """
        self.__save_lims(ax_letter)
        if scale == 'log' or scale == True:
            if self.xlims[ax_letter][0] < 0:
                lntresh = 10 ** (np.round(
                    min(np.log10(-self.xlims[ax_letter][0].value),
                    np.log10(self.xlims[ax_letter][1].value))) - 6)
                self.axd[ax_letter].set_xscale('symlog', linthresh=lntresh)
            else:
                self.axd[ax_letter].set_xscale('log')
        else:
            self.axd[ax_letter].set_xscale('linear')
        self.logX[ax_letter] = self.axd[ax_letter].get_xscale()
    
    def Yscale(self, scale, ax_letter="A"):
        """
        Sets the y axis scale of the plot at the cooresponding letter.
        If the scale is 'log' and the limits are negative, we use a 
        custom symlog scale. 
        """
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
        """
        Change the colormap of the plot at the corresponding letter, if
        it is a 2D plot. The plot then needs to be redone.
        """
        self.cmap_color[axd_letter] = cmap
        self.__redo_plot()

    def cbar_levels(self, cbar_levels, axd_letter="A"):
        """
        Change the colorbar levels of the plot at the corresponding
        letter, if it is a 2D plot. The plot then needs to be redone.
        """
        self.cbar_lv[axd_letter] = cbar_levels
        self.__redo_plot()
    
    def labels(self, xlabel, ylabel, axd_letter="A"):
        """
        Let's the user change the default labels of the plot at the
        corresponding letter. For whatever reason...
        """
        if xlabel is not None:
            self.axd[axd_letter].set_xlabel(xlabel)
        if ylabel is not None:
            self.axd[axd_letter].set_ylabel(ylabel)
        self.__save_labels(axd_letter)
    
    def update_legend(self, legend, axd_letter="A"):
        """
        Plot the legend of the plot at the corresponding letter. The
        number of legend entries must be the same as the number of
        lines in the plot.
        """
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

    def __update_params(self, ax_letter, plane, data, cbar_position, dim,
                        sim_dim, **kwargs):
        """
        Most important method of the class. It updates the dictionaries
        containing the plot parameters and data information.
        It is meat to be called BEFORE plotting.
        """
        if plane in  ['xy', 'yx', 'xz', 'zx', 'yz', 'zy']:
            if plane in ['yx', 'zx', 'zy']:
                plane = plane[::-1]
            plane = plane.upper()
            grid = (getattr(data, plane[0]), getattr(data, plane[1]))
        else:
            grid = getattr(data, plane)
        if ax_letter not in self.plot_dim:
            self.grid[ax_letter] = [grid]
            self.data[ax_letter] = [data.data]
            self.plot_dim[ax_letter] = [dim]
            self.sim_dimension[ax_letter] = [sim_dim]
            self.cbar_position[ax_letter] = cbar_position
            self.cbar_log[ax_letter] = data.data.log
            self.cbar_lv[ax_letter] = None
            self.cmap_color[ax_letter] = data.data.cmap
            self.cbar_label[ax_letter] = data.data.label
            if 'alpha' in kwargs:
                self.alpha[ax_letter] = [kwargs['alpha']]
            else:
                self.alpha[ax_letter] = [1]
            if 'color' in kwargs or 'c' in kwargs:
                self.line_color[ax_letter] = [kwargs.get('color',
                                              kwargs.get('c'))]
            else:
                self.line_color[ax_letter] = [None]
            if 'lw' in kwargs:
                self.lw[ax_letter] = [kwargs['lw']]
            else:
                self.lw[ax_letter] = [None]
        else:
            self.grid[ax_letter].append(grid)
            self.data[ax_letter].append([data.data])
            self.plot_dim[ax_letter].append(dim)
            self.sim_dimension[ax_letter].append(sim_dim)
            if 'alpha' in kwargs:
                self.alpha[ax_letter].append(kwargs['alpha'])
            else:
                self.alpha[ax_letter].append(1)
            if 'color' in kwargs or 'c' in kwargs:
                self.line_color[ax_letter].append(kwargs.get('color',
                                                         kwargs.get('c')))
            else:
                self.line_color[ax_letter].append(None)
            if 'lw' in kwargs:
                self.lw[ax_letter].append(kwargs['lw'])
            else:
                self.lw[ax_letter].append(None)
         
    def __update_cbar_position(self, ax_letter, cbar_position):
        """
        Changes the position of the colorbar of the plot.
        """
        self.cbar_position[ax_letter] = cbar_position
    
    def __update_fields_params(self, ax_letter, field, field_type):
        """
        If we want to plot field lines or arrows, we need to have the
        saved in the corresponding dictionary.
        """
        self.field[ax_letter] = field
        self.field_type[ax_letter] = field_type
    
    def __clear_param_key(self, ax_letter):
        """
        Clears the dictionaries containing the plot parameters.
        """
        keys = [self.plot_dim,
                self.grid,
                self.data,
                self.sim_dimension,
                self.cbar_position,
                self.cbar_log,
                self.cbar_lv,
                self.cmap_color,
                self.cbar_label,
                self.xlims,
                self.ylims,
                self.logX,
                self.logY,
                self.xlabels,
                self.ylabels,
                self.legend,
                self.alpha,
                self.line_color,
                self.lw,
                self.field,
                self.field_type]
        for key in keys:
            if ax_letter in key:
                key.pop(ax_letter)

    def __reset_params(self):
        """
        Initializes the dictionaries containing the plot parameters. Or
        clears them.
        """
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
        self.alpha = {}
        self.line_color = {}
        self.lw = {}
        self.field = {}
        self.field_type = {}
        self.sim_dimension = {}
    
    def __normalize_format_cbar(self, ax_letter):
        """
        Sets up the normalization and the tick format of the colorbar.
        Also returns a 1D array with the custom cbar levels.
        """
        self.cbar_lv[ax_letter] = list(self.cbar_lv[ax_letter])
        if self.cbar_lv[ax_letter][0] is None:
            self.cbar_lv[ax_letter][0] = np.nanmin(self.data[ax_letter])
        if self.cbar_lv[ax_letter][1] is None:
            self.cbar_lv[ax_letter][1] = np.nanmax(self.data[ax_letter])
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
        self.cbar_lv[ax_letter] = tuple(self.cbar_lv[ax_letter])
        return norm, fmt, cbar_levels

    def __plot2Dmesh(self, ax_letter, **kwargs):
        """
        Fills up the selected axes with a 2D mesh plot.
        """
        ## -2 is spectrogram, -3 is HHT spectrogram
        if (self.plot_dim[ax_letter].count(-2) + \
            self.plot_dim[ax_letter].count(-3)) != 1:
            raise ValueError('Too many mesh data')
        try:
            indx = self.plot_dim[ax_letter].index(-2)
        except:
            indx = self.plot_dim[ax_letter].index(-3)
        norm, fmt, _ = self.__normalize_format_cbar(ax_letter)
        if 'lmin' in kwargs:
            sh = 'auto'
        else:
            sh = 'gouraud'
        pcm = self.axd[ax_letter].pcolormesh(self.grid[ax_letter][indx][0],
                                            self.grid[ax_letter][indx][1],
                                            self.data[ax_letter][indx], norm=norm,
                                            cmap=self.cmap_color[ax_letter],
                                            shading=sh)
        cbar = self.fig.colorbar(pcm, cax=self.axd[ax_letter.lower()],
                                 format=ticker.FuncFormatter(fmt),
                                 location=cbar_loaction(
                                     self.cbar_position[ax_letter]))
        cbar.set_label(self.cbar_label[ax_letter])
        
    def __plot2Dfield(self, ax_letter):
        """
        Add a 2D field to the plot. It can be a velocity field or a
        magnetic field.
        """
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
        """
        Adds a contourf plot to the selected axes.
        """
        ## Check whether we have more than one 2D data per plot
        ## 2 is slice -1 is profile
        if (self.plot_dim[ax_letter].count(2) + \
            self.plot_dim[ax_letter].count(-1)) != 1:
            raise ValueError('Too many 2D data')
        norm, fmt, cbar_levels = self.__normalize_format_cbar(ax_letter)
        ## Find what to plot
        try:
            indx = self.plot_dim[ax_letter].index(2)
        except:
            indx = self.plot_dim[ax_letter].index(-1)
        
        pcm = self.axd[ax_letter].contourf(self.grid[ax_letter][indx][0],
                                            self.grid[ax_letter][indx][1],
                                            self.data[ax_letter][indx],
                                            norm=norm,
                                            levels=cbar_levels,
                                            antialiased=True,
                                            cmap=self.cmap_color[ax_letter],
                                            extend='both')

        cbar = self.fig.colorbar(pcm, cax=self.axd[ax_letter.lower()],
                                 format=ticker.FuncFormatter(fmt), 
                                location=cbar_loaction(
                                    self.cbar_position[ax_letter]))
        cbar.set_label(self.cbar_label[ax_letter])
        ## Moved the label to avoid overlapping with the cbar
        if self.cbar_position[ax_letter] in ['L', 'R'] and self.plot_dim == 2:
            self.axd[ax_letter].yaxis.labelpad = -10
        ## Hides every other tick label to avoid overlapping
        if self.cbar_position[ax_letter] in ['T', 'B']:
            for lb in cbar.ax.xaxis.get_ticklabels()[::2]:
                lb.set_visible(False)

    def __plot1D(self, ax_letter, redo=False):
        """
        Adds a 1D plot to the selected axes. If we pass a list of data,
        it will plot all of them.
        """
        if not redo:
            for indx in reversed(range(len(self.plot_dim[ax_letter]))):
                if self.plot_dim[ax_letter][indx] != 1:
                    continue
                kw = {'alpha': self.alpha[ax_letter][indx]}
                if self.lw[ax_letter][indx] is not None:
                    kw['lw'] = self.lw[ax_letter][indx]
                if self.line_color[ax_letter][indx] is not None:
                    kw['color'] = self.line_color[ax_letter][indx]
                if type(self.data[ax_letter][indx]) == list:
                    for data in self.data[ax_letter][indx]:
                        self.axd[ax_letter].plot(self.grid[ax_letter][indx],
                                                data, **kw)
                else:
                    self.axd[ax_letter].plot(self.grid[ax_letter][indx],
                                            self.data[ax_letter][indx],
                                            **kw)
                break
        else:
            for indx in range(len(self.plot_dim[ax_letter])):
                if self.plot_dim[ax_letter][indx] != 1:
                    continue
                kw = {'alpha': self.alpha[ax_letter][indx]}
                if self.lw[ax_letter][indx] is not None:
                    kw['lw'] = self.lw[ax_letter][indx]
                if self.line_color[ax_letter][indx] is not None:
                    kw['color'] = self.line_color[ax_letter][indx]
                if type(self.data[ax_letter][indx]) == list:
                    for data in self.data[ax_letter][indx]:
                        self.axd[ax_letter].plot(self.grid[ax_letter][indx],
                                                data, **kw)
                else:
                    self.axd[ax_letter].plot(self.grid[ax_letter][indx],
                                            self.data[ax_letter][indx],
                                            **kw)

    def __redo_plot(self):
        """
        We replot everything in the figure. What is done depends on the
        type of plot we are dealing with.
        """
        self._PlotCreation__close_figure()
        self._PlotCreation__setup_axd(self.number, self.form_factor)
        for ax_letter in self.axd:
            if ax_letter not in self.plot_dim.keys():
                continue
            if self.data[ax_letter] is None:
                continue
            dm = self.plot_dim[ax_letter]
            all_1D = all([d == 1 for d in dm])
            if 2 in dm:
                self.__plot2D(ax_letter)
                for indx in range(len(dm)):
                    if dm[indx] == 2:
                        sdim = self.sim_dimension[ax_letter][indx]
                        grid = self.grid[ax_letter][indx][0]
                        if np.isclose(np.abs(grid.max()), np.abs(grid.min())):
                            plane = "xy"
                        else:
                            plane = "xz"
                        break
                set2Dlims(self.axd, self.xlims[ax_letter], None, self.number,
                          self.form_factor, sdim, plane)
            elif -1 in dm:
                self.__plot2D(ax_letter)
                self.xlim(self.xlims[ax_letter], ax_letter)
                self.Yscale(self.logY[ax_letter], ax_letter)
                self.Xscale(self.logX[ax_letter], ax_letter)
            elif (-2 in dm) or (-3 in dm):
                self.__plot2Dmesh(ax_letter)
                self.xlim(self.xlims[ax_letter], ax_letter)
                self.ylim(self.ylims[ax_letter], ax_letter)
                self.Xscale(self.logX[ax_letter], ax_letter)
                self.Yscale(self.logY[ax_letter], ax_letter)
            if 1 in dm:
                self.__plot1D(ax_letter, redo=True)
                if all_1D:
                    self.ylim(self.ylims[ax_letter], ax_letter)
                    self.xlim(self.xlims[ax_letter], ax_letter)
                    self.Xscale(self.logX[ax_letter], ax_letter)
                    self.Yscale(self.logY[ax_letter], ax_letter)
            if ax_letter in self.legend:
                self.update_legend(self.legend[ax_letter], ax_letter)
            if ax_letter in self.field:
                self.__plot2Dfield(ax_letter)
            #self.labels(self.xlabels[ax_letter], self.ylabels[ax_letter],
            #            ax_letter)
        self._PlotCreation__setup_aspect()

    def __save_xlims(self, ax_letter=None):
        """
        Saves the x limits into the dictionary.
        """
        if ax_letter is None:
            for ax_letter in self.axd:
                xlims = self.axd[ax_letter].get_xlim()
                unit = self.axd[ax_letter].xaxis.get_units()[1]
                self.xlims[ax_letter] = (xlims[0] * unit, xlims[1] * unit)
        else:
            xlims = self.axd[ax_letter].get_xlim()
            unit = self.axd[ax_letter].xaxis.get_units()[1]
            self.xlims[ax_letter] = (xlims[0] * unit, xlims[1] * unit)
    
    def __save_ylims(self, ax_letter=None):
        """
        Saves the y limits into the dictionary.
        """
        if ax_letter is None:
            for ax_letter in self.axd:
                ylims = self.axd[ax_letter].get_ylim()
                unit = self.axd[ax_letter].yaxis.get_units()[1]
                self.ylims[ax_letter] = (ylims[0] * unit, ylims[1] * unit)
        else:
            ylims = self.axd[ax_letter].get_ylim()
            unit = self.axd[ax_letter].yaxis.get_units()[1]
            self.ylims[ax_letter] = (ylims[0] * unit, ylims[1] * unit)
    
    def __save_lims(self, ax_letter=None):
        """
        Saves the x and y limits into the dictionary.
        """
        self.__save_xlims(ax_letter)
        self.__save_ylims(ax_letter)
 
    def __save_labels(self, ax_letter=None):
        """
        Saves the x and y labels into the dictionary.
        """
        if ax_letter is None:
            for ax_letter in self.axd:
                xlabel = self.axd[ax_letter].xaxis.get_units()[0]
                ylabel = self.axd[ax_letter].yaxis.get_units()[0]
                self.xlabels[ax_letter] = xlabel
                self.ylabels[ax_letter] = ylabel
        else:
            xlabel = self.axd[ax_letter].xaxis.get_units()[0]
            ylabel = self.axd[ax_letter].yaxis.get_units()[0]
            self.xlabels[ax_letter] = xlabel
            self.ylabels[ax_letter] = ylabel

    def __save_scale(self, ax_letter=None):
        """
        Saves the x and y scales into the dictionary.
        """
        if ax_letter is None:
            for ax_letter in self.axd:
                self.logX[ax_letter] = self.axd[ax_letter].get_xscale()
                self.logY[ax_letter] = self.axd[ax_letter].get_yscale()
        else:
            self.logX[ax_letter] = self.axd[ax_letter].get_xscale()
            self.logY[ax_letter] = self.axd[ax_letter].get_yscale()
    
    def __save_params(self):
        """
        Saves the limits, scales and labels of the plots.
        """
        self.__save_lims()
        self.__save_scale()
        self.__save_labels()

