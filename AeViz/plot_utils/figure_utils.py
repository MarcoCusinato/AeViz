"""
This is a utility module. It contains helper functions to assist in
figure creation and manipulation.
"""

def cbar_loaction(loc):
    """
    Helper function, returns the location of the colorbar.
    """
    location = {'T': 'top', 'B': 'bottom', 'L': 'left', 'R': 'right'}
    return location[loc]

def return_positioning(number, form_factor):
    """
    We decide the axes positioning in the figure, the with and height
    of each subplot and the grid spacing. Everything is contained in
    dictionaries. First entry is the number of axes, second entry is
    the form factor of the figure. Only exception is the 5th figure,
    which is used to plot GWs spectra and strains.
    """
    positioning = {1: {1: """A""",
                       2: """A.a""",
                       4: """A.a""",
                       5:"""A..
                            B.b""",
                       6: """A.a"""
                                },
                   2: {1: """AB""",
                       2: """a.AB.b""",
                       3: """A
                             B""",
                       4: """A.a
                             B.b""",
                       5: """..A.C..
                             b.B.D.d""",
                       6: """a.A.B.b"""
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
                              f.F....""",
                        6: """a.A.B.b
                              .......
                              c.C...."""},
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
                              f.F.H.h""",
                        6: """a.A.B.b
                              .......
                              c.C.D.d"""},
                    5: {1:"""A..
                             B..
                             C..
                             D..
                             E.e"""}}
    
    width_ratio = {1: {1: [1],
                       2: [1.1, 0.4, 0.08],
                       4: [1.1, 0.1, 0.08],
                       5: [1.1, 0.1, 0.08],
                       6: [1.1, 0.4, 0.08]},
                   2: {1: [1, 1],
                       2: [0.08, 0.4, 1.1, 1.1, 0.4, 0.08],
                       3: [1],
                       4: [1.1, 0.1, 0.08],
                       5: [0.08, 0.2, 1, 0.2, 1, 0.1, 0.08],
                       6: [0.08, 0.4, 1.1, 0.4, 1.1, 0.4, 0.08]
                       },
                   3: {1: [1],
                       2: [0.08, 0.4, 1.1, 1.1, 0.4, 0.08],
                       4:[0.08, 0.4, 1.1, 0.4, 1.1, 0.4, 0.08],
                       5: [0.08, 0.2, 1, 0.2, 1, 0.1, 0.08],
                       6: [0.08, 0.4, 1.1, 0.4, 1.1, 0.4, 0.08]},
                   4: {1: [1, 1],
                       2: [1, 1],
                       4: [0.08, 0.4, 1.1, 0.4, 1.1, 0.4, 0.08],
                       5: [0.08, 0.2, 1, 0.2, 1, 0.1, 0.08],
                       6: [0.08, 0.4, 1.1, 0.4, 1.1, 0.4, 0.08]},
                    5: {1: [1, 0.03, 0.04]}}
    
    height_ratio = {1: {1: [1],
                        2: [1],
                        4: [1],
                        5: [0.25, 0.75],
                        6: [1]},
                    2: {1: [1],
                        2: [1],
                        3: [1, 1],
                        4: [1, 1],
                        5: [0.25, 0.75],
                        6: [1]},
                    3: {1: [1, 1, 1],
                        2: [0.25, 0.75],
                        4: [1, 1],
                        5: [0.25, 0.75, 0.2, 0.25, 0.75],
                        6: [1, 0.2, 1]},
                    4: {1: [1, 1],
                        2: [0.08, 0.4, 1.1, 1.1, 0.4, 0.08],
                        4: [1, 1],
                        5: [0.25, 0.75, 0.2, 0.25, 0.75],
                        6: [1, 0.2, 1]},
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
    """
    We decide the figure size based on the number of axes and the
    form factor of the figure.
    """
    fig_size = {1: {1: (8, 8),
                    2: (8, 8),
                    4: (6, 4),
                    5: (6, 6),
                    6: (8, 8)},
                2: {1: (16, 8),
                    2: (8, 8),
                    3: (8, 8),
                    4: (6, 8),
                    5: (12, 6),
                    6: (8, 8)},
                3: {1: (9, 8.42),
                    2: (9, 8.42),
                    4: (12, 8),
                    5: (12, 12),
                    6: (9, 8.42)},
                4: {1: (16, 16),
                    2: (8, 11.6),
                    4: (12, 8),
                    5: (12, 12),
                    6: (8, 11.6)},
                5: {1: (9, 16)}}
    
    return fig_size[number][form_factor]

def change_y_cbar_aspect(axd, axd_bar):
    """
    When setting the colorobars as a subplot, sometimes it appears
    taller than the axes. This function fixes that.
    """
    axd_position = axd.get_position()
    y0, y1 = axd_position.y0, axd_position.y1
    bar_position = axd_bar.get_position()
    bar_position.y0, bar_position.y1 = y0, y1
    axd_bar.set_position(bar_position)

def change_x_cbar_aspect(axd, axd_bar):
    """
    When setting the colorobars as a subplot, sometimes it appears
    wider than the axes. This function fixes that.
    """
    axd_position = axd.get_position()
    x0, x1 = axd_position.x0, axd_position.x1
    bar_position = axd_bar.get_position()
    bar_position.x0, bar_position.x1 = x0, x1
    axd_bar.set_position(bar_position)

def change_x_aspect(axd, axd_left, axd_right):
    """
    Ensures the x-dimensional aspect ratio is the same for all the axes.
    """
    axd_position = axd.get_position()
    axd_left_position = axd_left.get_position()
    axd_right_position = axd_right.get_position()
    axd_position.x0, axd_position.x1 = axd_left_position.x0, \
        axd_right_position.x1
    axd.set_position(axd_position)

def change_y_aspect(axd, axd_bottom, axd_top):
    """
    Ensures the y-dimensional aspect ratio is the same for all the axes.
    """
    axd_position = axd.get_position()
    axd_bottom_position = axd_bottom.get_position()
    axd_top_position = axd_top.get_position()
    axd_position.y0, axd_position.y1 = axd_bottom_position.y0,\
        axd_top_position.y1
    axd.set_position(axd_position)
