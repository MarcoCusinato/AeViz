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

def AeViz_plot_panel(AeViz, qt, file, projection, plane, **kwargs):
    """
    This plot a single panel/line on a figure.
    """
    kwargs.setdefault('spectrogram', False)
    kwargs.setdefault('histogram', False)
    if projection == '1D':
        if kwargs['histogram']:
            kwargs.setdefault('yquantity', 'mass')
            yquantity = kwargs['yquantity']
            kwargs.pop('yquantity')
            AeViz.plot1Dhistogram(file, qt, yquantity, **kwargs)
        elif plane == 'time':
            if hasattr(AeViz.loaded_data, 'global_' + qt):
                qt = 'global_' + qt
            if kwargs['spectrogram']:
                AeViz.plot1DSpectrogram(qt, **kwargs)
            else:
                return AeViz.plot1D(file, qt, plane, **kwargs)
        else:
            return AeViz.plot1D(file, qt, plane, **kwargs)
    elif projection == '2D':
        if file is not None:
            if type(plane) == tuple:
                plane = plane[0]
            if kwargs['histogram']:
                kwargs.setdefault('yquantity', qt)
                kwargs.setdefault('cquantity', 'mass')
                yquantity = kwargs['yquantity']
                cquantity = kwargs['cquantity']
                kwargs.pop('yquantity')
                kwargs.pop('cquantity')
                AeViz.plot2Dhistogram(file, qt, yquantity, cquantity, **kwargs)
            elif type(plane) != str:
                AeViz.plotHammer(file, plane, qt, **kwargs)
            else:
                AeViz.plot2D(file, plane, qt, **kwargs)
        elif plane == 'time':
            if kwargs['spectrogram']:
                if hasattr(AeViz.loaded_data, 'global_' + qt):
                    return AeViz.plotSpectrogram('global_' + qt, **kwargs)
                else:
                    return AeViz.plotSpectrogram(qt, **kwargs)
            else:
                AeViz.plotProfile(qt, **kwargs)
        else:
            AeViz.plot2D(file, plane, qt, **kwargs)

def AeViz_plot_radius_panel(AeViz, qt, projection, rad, **kwargs):
    kwargs.setdefault('spectrogram', False)
    kwargs.setdefault('rad', rad)
    if rad == 'full':
        kwargs.setdefault('plot', 'A')
    if rad == 'full':
        if kwargs['spectrogram']:
            raise TypeError('Cannot plot all radii in spectrogram mode.')
        if 2 not in AeViz.plot_dim[kwargs['plot']]:
            raise TypeError('Cannot plot full radius in a non 2 dimensional plot.')
        AeViz.plot1D_2Dradius(qt, **kwargs)
    else: 
        if kwargs['spectrogram'] and projection == '2D':
            AeViz.plotSpectrogram(qt, **kwargs)
        elif kwargs['spectrogram'] and projection == '1D':
            AeViz.plot1DSpectrogram(qt, **kwargs)
        else:
            AeViz.plot1D(None, qt, 'time', **kwargs)
