def plot_qt(AeViz, file, qt, projection, index1, index2, plane, **kwargs):
    """
    Plots the quantity or quantities qt from a specific file or
    simulation in a specific projection.
    """
    if type(qt) == str:
        if projection == '1D':
            AeViz.plot1D(file, qt, plane, index1, index2, **kwargs)
        elif projection == '2D':
            if plane == 'time':
                AeViz.plotProfile(qt, **kwargs)
            else:
                AeViz.plot2D(file, plane, qt, **kwargs)
    elif type(qt) == list:
        if projection == '1D':
            for q in qt:
                AeViz.plot1D(file, q, plane, index1, index2, **kwargs)
        elif projection == '2D':
            if plane == 'time':
                if len(qt) == 4:
                    AeViz.plotProfile(qt[0], qt[1], qt[2], qt[3], **kwargs)
                elif len(qt) == 3:
                    AeViz.plotProfile(qt[0], qt[1], qt[2], **kwargs)
                elif len(qt) == 2:
                    AeViz.plotProfile(qt[0], qt[1], **kwargs)
                else:
                    AeViz.plotProfile(qt[0], **kwargs)
            else:
                if len(qt) == 4:
                    AeViz.plot2D(file, plane, qt[0], qt[1], qt[2],
                                 qt[3], **kwargs)
                elif len(qt) == 3:
                    AeViz.plot2D(file, plane, qt[0], qt[1], qt[2], **kwargs)
                elif len(qt) == 2:
                    AeViz.plot2D(file, plane, qt[0], qt[1], **kwargs)
                else:
                    AeViz.plot2D(file, plane, qt[0], **kwargs)