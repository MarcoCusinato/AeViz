from AeViz.plot_utils.plot_creation import PlotCreation
from matplotlib import ticker
import numpy as np
from matplotlib.colors import LogNorm, SymLogNorm, Normalize
from AeViz.plot_utils.limits_utils import set2Dlims
from AeViz.plot_utils.figure_utils import cbar_loaction
from AeViz.units import aerray
from AeViz.units import u


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
        
    def to(self, unit, plot='A', axis=None):
        """
        change the units of the selected plot either of the main data or axis.add()
        unit: the unit to change to
        plot:the letter of the plot
        axis: the axis to target, if None, would target the main data
        """
        if axis is None:
            for i in range(len(self.data[plot])):
                try:
                    old_units = self.data[plot][i].unit
                    self.data[plot][i] = self.data[plot][i].to(unit)
                    if type(self.cbar_lv[plot]) == tuple:
                        self.cbar_lv[plot] = list(self.cbar_lv[plot])
                    for lv in range(len(self.cbar_lv[plot])):
                        self.cbar_lv[plot][lv] = self.cbar_lv[plot][lv] * \
                            old_units.to(unit)
                except:
                    pass
        else:
            if any([a == 2 for a in self.plot_dim[plot]]):
                axis = 'xy'
            if axis == 'x':
                for i in range(len(self.data[plot])):
                    try:
                        if type(self.grid[plot][i]) == tuple:
                            self.grid[plot][i] = list(self.grid[plot][i])
                        if (self.grid[plot][i]) == list:
                            old_units = self.grid[plot][i][0].unit
                            self.grid[plot][i][0] = self.grid[plot][i][0].to(unit)
                        else:
                            old_units = self.grid[plot][i].unit
                            self.grid[plot][i] = self.grid[plot][i].to(unit)
                    except:
                        pass
            elif axis == 'y':
                for i in range(len(self.data[plot])):
                    try:
                        if self.plot_dim[plot][i] == 1:
                            old_units = self.data[plot][i].unit
                            self.data[plot][i] = self.data[plot][i].to(unit)
                        else:
                            if type(self.grid[plot][i]) == tuple:
                                self.grid[plot][i] = list(self.grid[plot][i])
                            old_units = self.grid[plot][i][1].unit
                            self.grid[plot][i][1] = self.grid[plot][i][1].to(unit)
                    except:
                        pass
            elif axis in ['xy', 'yx']:
                for i in range(len(self.data[plot])):
                    try:
                        if self.plot_dim[plot][i] == 1:
                            old_units = self.data[plot][i].unit
                            self.data[plot][i] = self.data[plot][i].to(unit)
                            self.grid[plot][i] = self.grid[plot][i].to(unit)
                        else:
                            if type(self.grid[plot][i]) == tuple:
                                self.grid[plot][i] = list(self.grid[plot][i])
                            old_units = self.grid[plot][i][0].unit
                            self.grid[plot][i][1] = self.grid[plot][i][1].to(unit)
                            self.grid[plot][i][0] = self.grid[plot][i][0].to(unit)                    
                    except:
                        pass 
        self.__redo_plot()

    def xlim(self, xlim, axd_letter="A"):
        """
        Sets the x limits of the axes corresponding to the axd_letter
        plot. In case of 2D plots, it sets the limits of all the axes.
        The limits are saved in the xlims dictionary.
        """
        if xlim is None:
            self.xlims[axd_letter] = xlim
            return
        xlim = list(xlim)
        for i, xl in enumerate(xlim):
            if not isinstance(xl, aerray):
                xlim[i] = xlim[i] * self.axd[axd_letter].xaxis.get_units()[1]
            else:
                xlim[i] = xlim[i].to(self.axd[axd_letter].xaxis.get_units()[1])
        xlim = tuple(xlim)
        if 2 in self.plot_dim[axd_letter]:
            idx = self.plot_dim[axd_letter].index(2)
            plane = self.plot_dim[axd_letter][idx]
            sdim = self.sim_dimension[axd_letter][idx]
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
        if ylim is None:
            self.ylims[axd_letter] = ylim
            return
        ylim = list(ylim)
        for i, yl in enumerate(ylim):
            if not isinstance(yl, aerray):
                ylim[i] = ylim[i] * self.axd[axd_letter].yaxis.get_units()[1]
            else:
                ylim[i] = ylim[i].to(self.axd[axd_letter].yaxis.get_units()[1])
        ylim = tuple(ylim)
        if 2 in self.plot_dim[axd_letter]:
            idx = self.plot_dim[axd_letter].index(2)
            plane = self.plot_dim[axd_letter][idx]
            sdim = self.sim_dimension[axd_letter][idx]
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
        if scale is None:
            self.logX[ax_letter] = scale
            return
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
        if scale is None:
            self.logY[ax_letter] = scale
            return
        self.__save_lims(ax_letter)
        if scale in ['log', 'symlog'] or scale == True:
            if self.ylims[ax_letter][0] < 0:
                lntresh = 10 ** (np.round(
                    min(np.log10(-self.ylims[ax_letter][0].value),
                    np.log10(self.ylims[ax_letter][1].value))) - 6)
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

    def title(self, title, axd_letter='A'):
        """
        Set the title for the axes corresonding to the axd_letter plot.
        It is then saved in the titles dictionary.
        """
        if title is not None:
            self.axd[axd_letter].set_title(title)
        self.__save_title(axd_letter)

    def text(self, xt, yt, text, axd_letter='A', kwargs_text={}):
        """
        Write text on the plot corresponding to axd_letter.
        """
        
        if text is not None:
            self.axd[axd_letter].text(xt, yt, text, **kwargs_text)
        self.__save_text(xt, yt, text, kwargs_text, axd_letter)
    
    def labels(self, xlabel, ylabel, axd_letter="A"):
        """
        Let the user change the default labels of the plot at the
        corresponding letter. For whatever reason...
        """
        if xlabel is not None:
            unit = self.axd[axd_letter].xaxis.get_units()[1]
            self.axd[axd_letter].set_xlabel(xlabel + f' [{unit:latex}]')
        if ylabel is not None:
            unit = self.axd[axd_letter].yaxis.get_units()[1]
            self.axd[axd_letter].set_ylabel(ylabel + f' [{unit:latex}]')
        self.__save_labels(axd_letter)
    
    def update_legend(self, legend, axd_letter="A"):
        """
        Plot the legend of the plot at the corresponding letter. The
        number of legend entries must be the same as the number of
        lines in the plot.
        """
        if axd_letter in self.legend:
            old_legend = legend.copy()
            legend = self.legend[axd_letter]
            if len(self.axd[axd_letter].lines) != len(legend):
                if type(old_legend) == list:
                    for ll in old_legend:
                        legend.append(ll)
                else:
                    legend.append(old_legend)
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

    def __update_params(self, file=None, ax_letter=None, plane=None, data=None,
                        cbar_position=None, dim=None, sim_dim=None,
                        **kwargs):
        """
        Most important method of the class. It updates the dictionaries
        containing the plot parameters and data information.
        It is meat to be called BEFORE plotting.
        """
        if type(plane) == tuple:
            if type(plane[0]) == str:
                grid = (getattr(data, plane[0]), getattr(data, plane[1]))
            else:
                grid = getattr(data, data.return_axis_names()[0])
        elif plane in  ['xy', 'yx', 'xz', 'zx', 'yz', 'zy',
                        'xz_phi_avg', 'zx_phi_avg', 'yz_phi_avg', 'zy_phi_avg']:
            if plane in ['xz_phi_avg', 'zx_phi_avg', 'yz_phi_avg', 'zy_phi_avg']:
                plane = 'xz_phi_avg'
            if plane in ['yx', 'zx', 'zy']:
                plane = plane[::-1]
            plane = plane.upper()
            grid = (getattr(data, plane[0]), getattr(data, plane[1]))
        else:
            grid = getattr(data, plane)
        if ax_letter not in self.plot_dim:
            self.file[ax_letter] = [file]
            self.plane[ax_letter] = [plane]
            self.grid[ax_letter] = [grid]
            self.data[ax_letter] = [data.data]
            self.plot_dim[ax_letter] = [dim]
            self.sim_dimension[ax_letter] = [sim_dim]
            self.cbar_position[ax_letter] = cbar_position
            self.cbar_log[ax_letter] = data.data.log
            self.cbar_lv[ax_letter] = data.data.limits
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
            if 'ls' in kwargs:
                self.ls[ax_letter] = [kwargs['ls']]
            else:
                self.ls[ax_letter] = [None]
        else:
            self.file[ax_letter].append(file)
            self.plane[ax_letter].append(plane)
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
            if 'ls' in kwargs:
                self.ls[ax_letter].append(kwargs['ls'])
            else:
                self.ls[ax_letter].append(None)
         
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
        if ax_letter not in self.field:
            self.field[ax_letter] = [field]
            self.field_type[ax_letter] = [field_type]
        else:
            self.field[ax_letter].append(field)
            self.field_type[ax_letter].append(field_type)
    
    def __copy_param_key(self, ax_letter_in, ax_letter_out):
        """
        Clears the dictionaries containing the plot parameters.
        """
        keys = [
                self.file,
                self.plane,
                self.plot_dim,
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
                self.titles,
                self.legend,
                self.alpha,
                self.line_color,
                self.lw,
                self.field,
                self.field_type]
        for key in keys:
            if ax_letter_out in key:
                key[ax_letter_in] = key[ax_letter_out]
    
    def __clear_param_key(self, ax_letter):
        """
        Clears the dictionaries containing the plot parameters.
        """
        keys = [
                self.file,
                self.plane,
                self.plot_dim,
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
                self.titles,
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
        self.file = {}
        self.plane = {}
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
        self.titles = {}
        self.texts = {}
        self.legend = {}
        self.alpha = {}
        self.line_color = {}
        self.lw = {}
        self.ls = {}
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
                                            self.data[ax_letter][indx].value,
                                            norm=norm,
                                            cmap=self.cmap_color[ax_letter],
                                            shading=sh)
        cbar = self.fig.colorbar(pcm, cax=self.axd[ax_letter.lower()],
                                 format=ticker.FuncFormatter(fmt),
                                 location=cbar_loaction(
                                     self.cbar_position[ax_letter]))
        if self.data[ax_letter][indx].unit == u.dimensionless_unscaled:
            ulab = ''
        else:
            ulab = f' [{self.data[ax_letter][indx].unit:latex}]'
        cbar.set_label(self.cbar_label[ax_letter])
        self.set_labels(ax_letter + ulab)
        
    def __plot2Dfield(self, ax_letter, grid_number):
        """
        Add a 2D field to the plot. It can be a velocity field or a
        magnetic field.
        """
        for i in range(len(self.field_type[ax_letter])):
            if self.field_type[ax_letter][i] == 'v':
                skip = (slice(None, None, 5), slice(None, None, 5))
                self.axd[ax_letter].quiver(self.grid[ax_letter][grid_number][0][skip],
                                        self.grid[ax_letter][grid_number][1][skip],
                                        self.field[ax_letter][i][0][skip].value,
                                        self.field[ax_letter][i][1][skip].value,
                                        linewidths=0.01,
                                        color='black',
                                        angles='xy'
                                        )
            elif self.field_type[ax_letter][i] == 'B':
                self.axd[ax_letter].contour(self.grid[ax_letter][grid_number][0],
                                            self.grid[ax_letter][grid_number][1],
                                            self.field[ax_letter][i], 45,
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
                                            self.data[ax_letter][indx].value,
                                            norm=norm,
                                            levels=cbar_levels,
                                            antialiased=True,
                                            cmap=self.cmap_color[ax_letter],
                                            extend='both')

        cbar = self.fig.colorbar(pcm, cax=self.axd[ax_letter.lower()],
                                 format=ticker.FuncFormatter(fmt), 
                                location=cbar_loaction(
                                    self.cbar_position[ax_letter]))
        if self.data[ax_letter][indx].unit == u.dimensionless_unscaled:
            ulab = ''
        else:
            ulab = f' [{self.data[ax_letter][indx].unit:latex}]'
        cbar.set_label(self.cbar_label[ax_letter] + ulab)
        ## Moved the label to avoid overlapping with the cbar
        if self.cbar_position[ax_letter] in ['L', 'R'] and self.plot_dim == 2:
            self.axd[ax_letter].yaxis.labelpad = -10
        ## Hides every other tick label to avoid overlapping
        if self.cbar_position[ax_letter] in ['T', 'B']:
            for lb in cbar.ax.xaxis.get_ticklabels()[::2]:
                lb.set_visible(False)
        if self.axd[ax_letter].xaxis.get_units() is None:
            self.axd[ax_letter].xaxis.set_units((self.grid[ax_letter][indx][0].label,
                                                self.grid[ax_letter][indx][0].unit))
        if self.axd[ax_letter].yaxis.get_units() is None:
            self.axd[ax_letter].yaxis.set_units((self.grid[ax_letter][indx][1].label,
                                                self.grid[ax_letter][indx][1].unit))
        self.set_labels(ax_letter)

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
                if self.ls[ax_letter][indx] is not None:
                    kw['ls'] = self.ls[ax_letter][indx]
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
                if self.ls[ax_letter][indx] is not None:
                    kw['ls'] = self.ls[ax_letter][indx]
                if type(self.data[ax_letter][indx]) == list:
                    for data in self.data[ax_letter][indx]:
                        self.axd[ax_letter].plot(self.grid[ax_letter][indx],
                                                data, **kw)
                else:
                    self.axd[ax_letter].plot(self.grid[ax_letter][indx],
                                            self.data[ax_letter][indx],
                                            **kw)
        self.set_labels(ax_letter)

    def __redo_plot(self):
        """
        We replot everything in the figure. What is done depends on the
        type of plot we are dealing with.
        """
        if any([ax.name == 'hammer' for ax in self.axd.values()]):
            self._PlotCreation__close_figure()
            self._PlotCreation__setup_axd(self.number, self.form_factor,
                                          prj='hammer')
        else:
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
                        if np.isclose((np.abs(grid.max())).value,
                                      (np.abs(grid.min())).value):
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
                number = self.plot_dim[ax_letter].index(2)
                self.__plot2Dfield(ax_letter, number)
            if ax_letter in self.titles:
                self.title(self.titles[ax_letter], ax_letter)
            if ax_letter in self.texts:
                if isinstance(self.texts[ax_letter], list):
                    # Loop on the texts
                    for i in range(len(self.texts[ax_letter])):
                        if 'transform' in self.texts[ax_letter][i][3]:
                            self.texts[ax_letter][i][3]['transform'] = \
                                self.axd[ax_letter].transAxes
                        self.text(self.texts[ax_letter][i][0], \
                                self.texts[ax_letter][i][1], \
                                self.texts[ax_letter][i][2], \
                                kwargs_text=self.texts[ax_letter][i][3])
                else:
                    pass
            #self.labels(self.xlabels[ax_letter], self.ylabels[ax_letter],
            #            ax_letter)
        self._PlotCreation__setup_aspect()

    def __save_xlims(self, ax_letter=None):
        """
        Saves the x limits into the dictionary.
        """
        if ax_letter is None:
            for ax_letter in self.axd:
                if ax_letter.islower():
                    continue
                xlims = self.axd[ax_letter].get_xlim()
                unit = self.axd[ax_letter].xaxis.get_units()
                if unit is None:
                    continue
                else:
                    unit = unit[1]
                self.xlims[ax_letter] = (xlims[0] * unit, xlims[1] * unit)
        else:
            xlims = self.axd[ax_letter].get_xlim()
            unit = self.axd[ax_letter].xaxis.get_units()
            if unit is not None:
                unit = unit[1]
                self.xlims[ax_letter] = (xlims[0] * unit, xlims[1] * unit)
    
    def __save_ylims(self, ax_letter=None):
        """
        Saves the y limits into the dictionary.
        """
        if ax_letter is None:
            for ax_letter in self.axd:
                if ax_letter.islower():
                    continue
                ylims = self.axd[ax_letter].get_ylim()
                unit = self.axd[ax_letter].yaxis.get_units()
                if unit is None:
                    continue
                else:
                    unit = unit[1]
                self.ylims[ax_letter] = (ylims[0] * unit, ylims[1] * unit)
        else:
            ylims = self.axd[ax_letter].get_ylim()
            unit = self.axd[ax_letter].yaxis.get_units()
            if unit is not None:
                unit = unit[1]
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

    def __save_title(self, ax_letter=None):
        """
        Save the title
        """
        if ax_letter is None:
            for ax_letter in self.axd:
                title = self.axd[ax_letter].get_title()
                self.titles[ax_letter] = title
        else:
            title = self.axd[ax_letter].get_title()
            self.titles[ax_letter] = title

    def __save_text(self, xt, yt, text, kwargs, ax_letter):
        """
        Save the custom texts
        """
        t = ((xt, yt, text, kwargs))
        if ax_letter in self.texts:
            if not t in self.texts[ax_letter]:
                self.texts[ax_letter].append(t)
        else:
            self.texts[ax_letter] = [t]
        
    
    def set_labels(self, ax_letter=None):
        if ax_letter is None:
            for letter in self.axd:
                try:
                    un = self.axd[letter].xaxis.get_units()
                    if un[1] == u.dimensionless_unscaled:
                        lab = un[0]
                    else:
                        lab = un[0] + f' [{un[1]:latex}]'
                    self.axd[letter].set_xlabel(lab)
                except:
                    pass
                try:
                    un = self.axd[letter].yaxis.get_units()
                    if un[1] == u.dimensionless_unscaled:
                        lab = un[0]
                    else:
                        lab = un[0] + f' [{un[1]:latex}]'
                    self.axd[letter].set_ylabel(lab)
                except:
                    pass
        else:
            try:
                un = self.axd[ax_letter].xaxis.get_units()
                if un[1] == u.dimensionless_unscaled:
                    lab = un[0]
                else:
                    lab = un[0] + f' [{un[1]:latex}]'
                self.axd[ax_letter].set_xlabel(lab)
            except:
                pass
            try:
                un = self.axd[ax_letter].yaxis.get_units()
                if un[1] == u.dimensionless_unscaled:
                    lab = un[0]
                else:
                    lab = un[0] + f' [{un[1]:latex}]'
                self.axd[ax_letter].set_ylabel(lab)
            except:
                pass
                

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

