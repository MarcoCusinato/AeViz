import matplotlib.pyplot as plt

def return_positioning(number, form_factor):
    positioning = {1: {1: """A""",
                       2: """A.a""",
                       4: """A.a""",
                       5:"""A..
                            B.b"""},
                   2: {1: """AB""",
                       2: """a.AB.b""",
                       3: """A
                             B""",
                       4: """A.a
                             B.b""",
                       5: """..A.C..
                             b.B.D.d"""
                       },
                    3: {1:"""A
                             B
                             C""",
                        2:"""..AA..
                             b.BC.c""",
                        4: """..AAA.a
                              b.B.C.c""",
                        5: """..A.C..
                              b.B.D.d
                              .......
                              ..E....
                              f.F...."""},
                    4: {1:"""AB
                             CD""",
                        2:"""ac
                             ..
                             AC
                             BD
                             ..
                             bd""",
                        4: """a.A.B.b
                              c.C.D.d""",
                        5: """..A.C..
                              b.B.D.d
                              .......
                              ..E.G..
                              f.F.H.h"""},
                    5: {1:"""A..
                             B..
                             C..
                             D..
                             E.e"""}}
    
    width_ratio = {1: {1: [1],
                       2: [1.1, 0.4, 0.08],
                       4: [1.1, 0.1, 0.08],
                       5: [1.1, 0.1, 0.08]},
                   2: {1: [1, 1],
                       2: [0.08, 0.4, 1.1, 1.1, 0.4, 0.08],
                       3: [1],
                       4: [1.1, 0.1, 0.08],
                       5: [0.08, 0.2, 1, 0.2, 1, 0.1, 0.08]},
                   3: {1: [1],
                       2: [0.08, 0.4, 1.1, 1.1, 0.4, 0.08],
                       4:[0.08, 0.4, 1.1, 0.4, 1.1, 0.4, 0.08],
                       5: [0.08, 0.2, 1, 0.2, 1, 0.1, 0.08]},
                   4: {1: [1, 1],
                       2: [1, 1],
                       4: [0.08, 0.4, 1.1, 0.4, 1.1, 0.4, 0.08],
                       5: [0.08, 0.2, 1, 0.2, 1, 0.1, 0.08]},
                    5: {1: [1, 0.03, 0.04]}}
    
    height_ratio = {1: {1: [1],
                        2: [1],
                        4: [1],
                        5: [0.25, 0.75]},
                    2: {1: [1],
                        2: [1],
                        3: [1, 1],
                        4: [1, 1],
                        5: [0.25, 0.75]},
                    3: {1: [1, 1, 1],
                        2: [0.25, 0.75],
                        4: [1, 1],
                        5: [0.25, 0.75, 0.2, 0.25, 0.75]},
                    4: {1: [1, 1],
                        2: [0.08, 0.4, 1.1, 1.1, 0.4, 0.08],
                        4: [1, 1],
                        5: [0.25, 0.75, 0.2, 0.25, 0.75]},
                    5: {1: [0.5, 0.5, 0.5, 0.5, 2]}}
    if form_factor == 4:
        gs_kw = {"wspace": 0,
                 "hspace": 0.2}
    elif '.' in positioning[number][form_factor]:
        gs_kw = {"wspace": 0,
                 "hspace": 0}
    else:
        gs_kw = {"wspace": 0.2,
                 "hspace": 0.2}
    return positioning[number][form_factor], width_ratio[number][form_factor],\
            height_ratio[number][form_factor], gs_kw


def return_fig_size(number, form_factor):
    fig_size = {1: {1: (8, 8),
                    2: (8, 8),
                    4: (6, 4),
                    5: (6, 6)},
                2: {1: (16, 8),
                    2: (8, 8),
                    3: (8, 8),
                    4: (6, 8),
                    5: (12, 6)},
                3: {1: (9, 8.42),
                    2: (9, 8.42),
                    4: (12, 8),
                    5: (12, 12)},
                4: {1: (16, 16),
                    2: (8, 11.6),
                    4: (12, 8),
                    5: (12, 12)},
                5: {1: (9, 16)}}
    return fig_size[number][form_factor]

def change_y_cbar_aspect(axd, axd_bar):
    axd_position = axd.get_position()
    y0, y1 = axd_position.y0, axd_position.y1
    bar_position = axd_bar.get_position()
    bar_position.y0, bar_position.y1 = y0, y1
    axd_bar.set_position(bar_position)

def change_x_cbar_aspect(axd, axd_bar):
    axd_position = axd.get_position()
    x0, x1 = axd_position.x0, axd_position.x1
    bar_position = axd_bar.get_position()
    bar_position.x0, bar_position.x1 = x0, x1
    axd_bar.set_position(bar_position)

def change_x_aspect(axd, axd_left, axd_right):
    axd_position = axd.get_position()
    axd_left_position = axd_left.get_position()
    axd_right_position = axd_right.get_position()
    axd_position.x0, axd_position.x1 = axd_left_position.x0, \
        axd_right_position.x1
    axd.set_position(axd_position)

def change_y_aspect(axd, axd_bottom, axd_top):
    axd_position = axd.get_position()
    axd_bottom_position = axd_bottom.get_position()
    axd_top_position = axd_top.get_position()
    axd_position.y0, axd_position.y1 = axd_bottom_position.y0,\
        axd_top_position.y1
    axd.set_position(axd_position)

class PlotCreation(object):
    def __init__(self):
        self.fig = None
        self.axd = None
        self.number = None
        self.form_factor = None


    def __setup_figure(self):
        if not self.fig_is_open():
            self.fig = plt.figure()
    
    def fig_is_open(self):
        if self.fig is None:
            return False
        else:
            return True
    
    def axd_is_open(self):
        if self.axd is None:
            return False
        else:
            return True
    
    def __setup_axd(self, number=1, form_factor=1):
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
        
    
    def change_figsize(self, multiplies=1):
        if self.fig_is_open():
            self.fig.set_size_inches(self.fig.get_size_inches()*multiplies)
            
    def __close_figure(self):
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
        if x:
            for ax in list_axs_share:
                ax.sharex(axd)
        if y:
            for ax in list_axs_share:
                ax.sharey(axd)
    
    def __setup_label_position(self):
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
    
    
