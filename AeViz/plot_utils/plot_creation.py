import matplotlib.pyplot as plt
from AeViz.plot_utils.figure_utils import (return_positioning, return_fig_size,
                                           change_x_cbar_aspect,
                                           change_y_cbar_aspect,
                                           change_x_aspect,
                                           change_y_aspect)

class PlotCreation(object):
    """
    Parent class for creating plots. Most of its methods are internal
    subroutines, meant just to make ''easier'' the plotting of the data.

    """
    def __init__(self):
        """
        We do not need any parameters to initialize the class.
        When called it initializes empty figure, axes dictionary and
        plot number and form factor.
        """
        self.fig = None
        self.axd = None
        self.number = None
        self.form_factor = None

    def fig_is_open(self):
        """
        Check if the figure has been created
        """
        if self.fig is None:
            return False
        else:
            return True
    
    def axd_is_open(self):
        """
        Check if the axes dictionary has been created
        """
        if self.axd is None:
            return False
        else:
            return True
        
    def change_figsize(self, multiplies=1):
        """
        This is a reskin of the matplotlib function set_size_inches.
        However, it preserves the aspect ratio of the figure.
        """
        if self.fig_is_open():
            self.fig.set_size_inches(self.fig.get_size_inches() * multiplies)
    
    def __setup_figure(self):
        """
        If the figure is not open, we create it.
        """
        if not self.fig_is_open():
            self.fig = plt.figure()
    
    def __setup_axd(self, number=1, form_factor=1):
        """
        Most important method of the class.
        What it does:
        - sets up the number of sybplot and their positioning
        - shares the axis if necessary
        - sets up the label position
        - sets up the aspect ratio
        - sets up the figure size
        """
        self.number = number
        self.form_factor = form_factor
        if self.fig_is_open():
            self.__close_figure()
        if not self.fig_is_open():
            self.__setup_figure()
        if not self.axd_is_open():
            position, width_ratio, height_ratio, gs_kw = \
                return_positioning(self.number, self.form_factor)
            self.axd = self.fig.subplot_mosaic(position,
                                               gridspec_kw=gs_kw,
                                               width_ratios=width_ratio, 
                                               height_ratios=height_ratio)
            if self.number == 1:
                if self.form_factor == 5:
                    self.__share_axis(self.axd["A"], [self.axd["B"]], True,
                                      False)
            elif self.number  == 2:
                if self.form_factor == 2:
                    self.__share_axis(self.axd["A"], [self.axd["B"]], False,
                                      True)
                elif self.form_factor == 5:
                    self.__share_axis(self.axd["A"], [self.axd["B"], 
                                                    self.axd["C"],
                                                    self.axd["D"]],
                                        True, False)
            elif self.number == 3:
                if self.form_factor == 2:
                    self.__share_axis(self.axd["B"], [self.axd["C"]], False,
                                      True)
                elif self.form_factor == 5:
                    self.__share_axis(self.axd["A"], [self.axd["B"], 
                                                    self.axd["C"],
                                                    self.axd["D"],
                                                    self.axd["E"],
                                                    self.axd["F"]],
                                        True, False)
            elif self.number == 4:
                if self.form_factor == 2:
                    self.__share_axis(self.axd["A"], [self.axd["B"]], True,
                                      False)
                    self.__share_axis(self.axd["C"], [self.axd["D"]], True,
                                      False)
                    self.__share_axis(self.axd["A"], [self.axd["C"]], False,
                                      True)
                    self.__share_axis(self.axd["B"], [self.axd["D"]], False,
                                      True)
                elif self.form_factor == 5:
                    for ax_letter in self.axd:
                        self.__share_axis(self.axd["A"], [self.axd["B"], 
                                                    self.axd["C"],
                                                    self.axd["D"],
                                                    self.axd["E"],
                                                    self.axd["F"],
                                                    self.axd["G"],
                                                    self.axd["H"]],
                                        True, False)
            elif self.number == 5:
                if self.form_factor == 1:
                    self.__share_axis(self.axd["A"], [self.axd["B"], 
                                                      self.axd["C"],
                                                      self.axd["D"],
                                                      self.axd["E"]],
                                        True, False)
            self.__setup_label_position()
            self.fig.set_size_inches(return_fig_size(self.number,
                                                     self.form_factor))
    
    def __close_figure(self):
        """
        Closes the figure if it is open, and restores the figure and axes
        dictionary to empty values.
        """
        if self.fig_is_open():
            plt.close(self.fig)
            self.fig = None
            self.axd = None
    
    def __setup_aspect(self):
        if self.number == 1 and self.form_factor == 2:
            change_y_cbar_aspect(self.axd["A"], self.axd["a"])
        elif (self.number == 2 and self.form_factor == 2):
            change_y_cbar_aspect(self.axd["A"], self.axd["a"])
            change_y_cbar_aspect(self.axd["B"], self.axd["b"])
        elif self.number == 3 and self.form_factor == 2:
            change_y_cbar_aspect(self.axd["B"], self.axd["b"])
            change_y_cbar_aspect(self.axd["B"], self.axd["c"])
            change_x_aspect(self.axd["A"], self.axd["B"], self.axd["C"])
        elif self.number == 4 and self.form_factor == 2:
            change_x_cbar_aspect(self.axd["A"], self.axd["a"])
            change_x_cbar_aspect(self.axd["B"], self.axd["b"])
            change_x_cbar_aspect(self.axd["C"], self.axd["c"])
            change_x_cbar_aspect(self.axd["D"], self.axd["d"])
    
    def __share_axis(self, axd, list_axs_share, x, y):
        """
        Given an axes and a list of axes, we share the x and/or y axis.
        """
        if x:
            for ax in list_axs_share:
                ax.sharex(axd)
        if y:
            for ax in list_axs_share:
                ax.sharey(axd)
    
    def __setup_label_position(self):
        """
        Given the number of plots and form factor of the plot,
        we set the label position, and remove the ticks if necessary.
        """
        if self.number == 2:
            if self.form_factor == 2:
                self.axd["A"].tick_params(top=False, labeltop=False,
                        bottom=True, labelbottom=True,
                        left=True, labelleft=True,
                        right=False, labelright=False)
                self.axd["B"].tick_params(top=False, labeltop=False,
                        bottom=True, labelbottom=True,
                        left=False, labelleft=False,
                        right=True, labelright=True)
                self.axd["A"].yaxis.set_label_position('left')
                self.axd["B"].yaxis.set_label_position('right')
        elif self.number == 3 and self.form_factor == 2:
            self.axd["A"].tick_params(top=True, labeltop=True,
                        bottom=False, labelbottom=False,
                        left=True, labelleft=True,
                        right=False, labelright=False)
            self.axd["B"].tick_params(top=False, labeltop=False,
                        bottom=True, labelbottom=True,
                        left=True, labelleft=True,
                        right=False, labelright=False)
            self.axd["C"].tick_params(top=False, labeltop=False,
                        bottom=True, labelbottom=True,
                        left=False, labelleft=False,
                        right=True, labelright=True)
            self.axd["A"].xaxis.set_label_position('top')
            self.axd["B"].xaxis.set_label_position('bottom')
            self.axd["C"].xaxis.set_label_position('bottom')
            self.axd["B"].yaxis.set_label_position('left')
            self.axd["C"].yaxis.set_label_position('right')
        elif self.number == 4 and self.form_factor == 2:
            self.axd["A"].tick_params(top=True, labeltop=True,
                        bottom=False, labelbottom=False,
                        left=True, labelleft=True,
                        right=False, labelright=False)
            self.axd["B"].tick_params(top=False, labeltop=False,
                        bottom=True, labelbottom=True,
                        left=True, labelleft=True,
                        right=False, labelright=False)
            self.axd["C"].tick_params(top=True, labeltop=True,
                        bottom=False, labelbottom=False,
                        left=False, labelleft=False,
                        right=True, labelright=True)
            self.axd["D"].tick_params(top=False, labeltop=False,
                        bottom=True, labelbottom=True,
                        left=False, labelleft=False,
                        right=True, labelright=True)
            self.axd["A"].xaxis.set_label_position('top')
            self.axd["B"].xaxis.set_label_position('bottom')
            self.axd["C"].xaxis.set_label_position('top')
            self.axd["D"].xaxis.set_label_position('bottom')
            self.axd["A"].yaxis.set_label_position('left')
            self.axd["B"].yaxis.set_label_position('left')
            self.axd["C"].yaxis.set_label_position('right')
            self.axd["D"].yaxis.set_label_position('right')
        elif self.number == 5 and self.form_factor == 1:
                self.axd["A"].tick_params(top=False, labeltop=False,
                                          bottom=False, labelbottom=False,
                                          left=True, labelleft=True,
                                          right=False, labelright=False)
                self.axd["B"].tick_params(top=False, labeltop=False,
                                            bottom=False, labelbottom=False,
                                            left=True, labelleft=True,
                                            right=False, labelright=False)
                self.axd["C"].tick_params(top=False, labeltop=False,
                                            bottom=False, labelbottom=False,
                                            left=True, labelleft=True,
                                            right=False, labelright=False)
                self.axd["D"].tick_params(top=False, labeltop=False,
                                            bottom=True, labelbottom=True,
                                            left=True, labelleft=True,
                                            right=False, labelright=False)
    
    
