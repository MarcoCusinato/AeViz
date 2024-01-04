import matplotlib.pyplot as plt

def return_positioning(number, form_factor):
    positioning = {1: {1: """A""",
                       2: """A.a"""},
                   2: {1: """AB""",
                       2: """a.AB.b""",
                       3: """A
                             B""",
                       4: """A.a
                             B.b""",
                       },
                    3: {1:"""..AA..
                           b.BC.c""",
                        2:"""AA
                           BC""",},
                    4: {1:"""AB
                           CD""",
                        2:"""ab
                             ..
                             AB
                             CD
                             ..
                             cd"""}}
    width_ratio = {1: {1: [1],
                       2: [1.1, 0.4, 0.08]},
                   2: {1: [1, 1],
                       2: [0.08, 0.4, 1.1, 1.1, 0.4, 0.08],
                       3: [1],
                       4: [1.1, 0.4, 0.08]},
                   3: {1:  [0.08, 0.4, 1.1, 1.1, 0.4, 0.08],
                       2: [1, 1]},
                   4: {1: [1, 1],
                       2: [1, 1]}}
    
    height_ratio = {1: {1: [1],
                        2: [1]},
                    2: {1: [1],
                        2: [1],
                        3: [1, 1],
                        4: [1, 1]},
                    3: {1: [0.25, 0.75],
                        2: [0.25, 0.75]},
                    4: {1: [1, 1],
                        2: [0.08, 0.4, 1.1, 1.1, 0.4, 0.08]}}
    if '.' in positioning[number][form_factor]:
        gs_kw = {"wspace": 0,
                 "hspace": 0}
    else:
        gs_kw = {"wspace": 0.2,
                 "hspace": 0.2}
    return positioning[number][form_factor], width_ratio[number][form_factor], height_ratio[number][form_factor], gs_kw
    

class PlotCreation:
    def __init__(self):
        self.fig = None
        self.axd = None
        self.fig_is_open = self.fig_is_open()
        self.axd_is_open = self.axd_is_open()

    def setup_figure(self):
        if not self.fig_is_open:
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
    
    def setup_axd(self, number=1, form_factor=1):
        if not self.axd_is_open:
            position, width_ratio, height_ratio, gs_kw = return_positioning(number, form_factor)
            self.axd = self.fig.subplot_mosaic(position,
                                               gridspec_kw=gs_kw,
                                               width_ratios=width_ratio, 
                                               height_ratios=height_ratio)