import matplotlib.pyplot as plt
from AeViz.utils.math_utils import function_average
from AeViz.quantities_plotting import TERMINAL
from AeViz.units import u
import numpy as np

def recognize_quantity(qt1, qt2, qt3, qt4, pars):
    """
    Distinguish between the different quantities that can be plotted.
    In particular, we deisentangle quantities that have different
    components.
    """
    if qt1 == 'NUE':
        qt1, qt2, qt3, qt4 = 'NUEX', 'NUEY', 'NUEZ', 'NUEE'
    if qt1 == 'NUA':
        qt1, qt2, qt3, qt4 = 'NUAX', 'NUAY', 'NUAZ', 'NUAE'
    if qt1 == 'NUX':
        qt1, qt2, qt3, qt4 = 'NUXX', 'NUXY', 'NUXZ', 'NUXE'
    if qt1 == 'B':
        qt1, qt2, qt3, qt4 = 'BX', 'BY', 'BZ', None
    if pars:
        if qt1 == 'VEL':
            qt1, qt2, qt3, qt4 = 'VELX', 'VELY', 'VELZ', None
        if qt1 == 'V':
            qt1, qt2, qt3, qt4 = 'VX', 'VY', 'VZ', None
        if qt1 == 'X':
            qt1, qt2, qt3, qt4 = 'X_n', 'X_p', 'X_alpha', 'X_h'
        if qt1 == 'CPOT':
            qt1, qt2, qt3, qt4 = 'CPOT_e', 'CPOT_n', 'CPOT_p', 'CPOT_nu'
    return qt1, qt2, qt3, qt4

def setup_cbars(qt1, qt2, qt3, qt4):
    """
    Returns the colorbars configuration for the contourf plots, given
    the number of quantities to plot.
    """
    if qt4 is not None or qt3 is not None:
        number, form_factor = 4, 2
        cbars = {'A': 'T', 'B': 'B', 'C': 'T', 'D': 'B'}
    elif qt2 is not None:
        number, form_factor = 2, 2
        cbars = {'A': 'L', 'B': 'R'}
    elif qt1 is not None:
        number, form_factor = 1, 2
        cbars = {'A': 'R'}
    else:
        raise ValueError('No quantity given.')
    return number, form_factor, cbars

def setup_cbars_profile(qt1, qt2, qt3, qt4):
    """
    Returns the colorbars configuration for the profile plots, given
    the number of quantities to plot.
    """
    if qt4 is not None:
        number, form_factor = 4, 4
        cbars = {'A': 'L', 'B': 'R', 'C': 'L', 'D': 'R'}
    elif qt3 is not None:
        number, form_factor = 3, 4
        cbars = {'A': 'R', 'B': 'L', 'C': 'R'}
    elif qt2 is not None:
        number, form_factor = 2, 4
        cbars = {'A': 'R', 'B': 'R'}
    elif qt1 is not None:
        number, form_factor = 1, 4
        cbars = {'A': 'R'}
    else:
        raise ValueError('No quantity given.')
    return number, form_factor, cbars

def setup_cbars_spectrogram(number):
    """
    Returns the colorbars configuration for the Fourier spectra plots,
    given the number of quantities to plot.
    """
    if number == 1:
        cbars = {'B': 'R'}
        plots = ["A", "B"]
    elif number == 2:
        cbars = {'B': 'L', 'D': 'R'}
        plots = ["C", "D"]
    elif number == 3:
        cbars = {'B': 'L', 'D': 'R', 'F': 'L'}
        plots = ["E", "F"]
    elif number == 4:
        cbars = {'B': 'L', 'D': 'R', 'F': 'L', 'H': 'R'}
        plots = ["G", "H"]
    return cbars, plots

def setup_cbars_HHT(number):
    """
    Returns the colorbars configuration for the Hilbert spectra plots,
    given the number of quantities to plot.
    """
    if number == 1:
        plot = "A"
        cbars = {'A': 'R'}
    elif number == 2:
        plot = "B"
        cbars = {'A': 'L', 'B': 'R'}
    elif number == 3:
        plot = "C"
        cbars = {'A': 'L', 'B': 'R', 'C': 'L'}
    elif number == 4:
        plot = "D"
        cbars = {'A': 'L', 'B': 'R', 'C': 'L', 'D': 'R'}
    return plot, cbars

def normalize_indices(index1, index2):
    """
    Gets the indices in the right format, which is a list or a single
    value.
    """
    if type(index1) == range:
        index1 = list(index1)
    if type(index2) == range:
        index2 = list(index2)
    if type(index1) == list and type(index2) == list:
        index2 = index2[0]
    return index1, index2

def get_data_to_plot(index1, index2, post_data, xaxis, dV):
    if post_data.ndim == 1:
        data = post_data
    elif post_data.ndim == 2:
        if xaxis == 'radius':
            if type(index1) == list:
                data = []
                for i in index1:
                    data.append(post_data[i, :])
            elif index1 is None:
                data = function_average(post_data, 2, 'Omega', dV[1][:, None] \
                    * dV[2])
            else:
                data = post_data[index1, :]
        elif xaxis == 'theta':
            
            if type(index1) == list:
                data = []
                for i in index1:
                    data.append(post_data[:, i])
            elif index1 is None:
                data = function_average(post_data, 2, 'theta', dV[0][None, :] \
                    * dV[2])
            else:
                data = post_data[:, index1]
    elif post_data.ndim == 3:
        if xaxis == 'radius':
            if type(index1) == list:
                data = []
                for i in index1:
                    data.append(post_data[i, index2, :])
            elif type(index2) == list:
                data = []
                for i in index2:
                    data.append(post_data[index1, i, :])
            elif index1 is None and index2 is None:
                data = function_average(post_data, 3, 'Omega', 
                                        dV[2][:, None, None] * \
                                            dV[1][None, :, None])
            else:
                data = post_data[index1, index2, :]
        elif xaxis == 'theta':
            if type(index1) == list:
                data = []
                for i in index1:
                    data.append(post_data[i, :, index2])
            elif type(index2) == list:
                data = []
                for i in index2:
                    data.append(post_data[index1, :, i])
            elif index1 is None and index2 is None:
                data = function_average(post_data, 3, 'theta',
                                        dV[2][:, None, None] \
                                            * dV[0][None, None, :])
            else:
                data = post_data[index1, :, index2]
        elif xaxis == 'phi':
            if type(index1) == list:
                data = []
                for i in index1:
                    data.append(post_data[:, i, index2])
            elif type(index2) == list:
                data = []
                for i in index2:
                    data.append(post_data[:, index1, i])
            elif index1 is None and index2 is None:
                data = function_average(post_data, 3, 'phi',
                                        dV[1][None, :, None] \
                                            * dV[0][None, None, :])
            else:
                data = post_data[:, index1, index2]
    return data

def show_figure():
    """
    Show the figure if the module is imported from the terminal.
    """
    if TERMINAL:
        plt.ion()
        plt.show()
    else:
        plt.ioff()

def get_qt_for_label(qt, **kwargs):
    """
    Get the label for the quantity to plot.
    """
    if 'der' not in kwargs:
        return qt
    if kwargs['der'] in ['r', 'radius']:
        return 'dr_' + qt
    elif kwargs['der'] in ['theta', 'th']:
        return 'dtheta_' + qt
    elif kwargs['der'] in ['phi', 'ph']:
        return 'dphi_' + qt
    elif kwargs['der'] in ['t', 'time']:
        return 'dt_' + qt
    else:
        return qt

def plot_panel(plotting_object, letter, file, quantity, cbars, plane,
               **kwargs):
    """
    Single panel plotting for contourf plots.
    Inputs:
    - plotting_object: Plotting class object.
    - letter: Panel letter.
    - file: File to read the data from.
    - quantity: Quantity to plot.
    - grid: tuple X, Y grid to plot on.
    - indices: tuple of the two indices to plot.
    - cbars: dictionary with the colorbars configuration.
    - labels: dictionary with the labels configuration.
    - plane: plane to plot.
    - **kwargs: additional arguments.
    """
    ## GET THE DATA
    if not 'plane' in kwargs:
        kwargs['plane'] = plane
    data = plotting_object._Data__get_data_from_name(quantity, file, **kwargs)
    plotting_object._PlottingUtils__update_params(
                                            file=file,
                                            ax_letter=letter,
                                            plane=plane,
                                            data=data,
                                            cbar_position=cbars[letter],
                                            dim=2,
                                            sim_dim=plotting_object.sim_dim)
    ## MAKE THE PLOT
    plotting_object._PlottingUtils__plot2D(letter)
    plot_number = plotting_object.plot_dim[letter].index(2)
    plotting_object.xlim(tuple(
        plotting_object.grid[letter][plotting_object.plot_dim['A'].index(2)][0].limits),
                         letter)
    plotting_object.Xscale(plotting_object.grid[letter][plot_number][0].log,
                           letter)
    plotting_object.Yscale(plotting_object.grid[letter][plot_number][1].log,
                           letter)
    plotting_object._PlottingUtils__save_labels(letter)

def plot_profile_panel(plotting_object, letter, quantity, cbars, **kwargs):
    kwargs["mesh"] = True
    data = plotting_object._Data__get_profile(quantity, **kwargs)
    """
    if 'rho_spherical_harmonics' in quantity:
        lab = labels[quantity]['label'].replace('XX', str(kwargs['l']))
        if kwargs['m'] is not None:
            lab = lab.replace('YY', str(kwargs['m']))
        else:
            lab = lab.replace('YY', '')
    else:
        lab = labels[quantity]['label']
    """
    plotting_object._PlottingUtils__update_params(
                                                  ax_letter=letter,
                                                  plane=('time', 'radius'),
                                                  data=data,
                                                  cbar_position=cbars[letter],
                                                  dim=-1,
                                                  sim_dim=plotting_object.sim_dim)
    plotting_object._PlottingUtils__plot2D(letter)
    plot_number = plotting_object.plot_dim[letter].index(-1)
    plotting_object.xlim(tuple(
        plotting_object.grid[letter][plot_number][0].limits), letter)
    plotting_object.ylim(tuple(
        plotting_object.grid[letter][plot_number][1].limits), letter)
    plotting_object.Xscale(plotting_object.grid[letter][plot_number][0].log,
                           letter)
    plotting_object.Yscale(plotting_object.grid[letter][plot_number][1].log,
                           letter)

def plot_hammer_panel(plotting_object, letter, file, quantity, cbars, plane,
                      **kwargs):
    if not 'plane' in kwargs:
        kwargs['plane'] = plane
    data = plotting_object._Data__get_data_from_name(quantity, file, **kwargs)
    plotting_object._PlottingUtils__update_params(
                                                  ax_letter=letter,
                                                  plane=('phi', 'theta'),
                                                  data=data,
                                                  cbar_position=cbars[letter],
                                                  dim=-2,
                                                  sim_dim=plotting_object.sim_dim)
    plotting_object._PlottingUtils__plot2Dmesh(letter)
    plotting_object.xlim(None, letter)
    plotting_object.ylim(None, letter)
    plotting_object.Xscale(None, letter)
    plotting_object.Yscale(None, letter)
    
def remove_labelling(data, no_nu):
    label = data.data.label
    new_label = label
    remove_symbols = [r'\mathrm{max}', r'\mathrm{min}', r'\mathrm{avg}',
                      r'\mathrm{tot}',
                       ',max', ',min', ',avg', ',tot', ', max', ', min', ', avg',
                       ', tot', ',x', ',y', ',z', ', x', ', y', ', z', ',r',
                       r',\theta', r',\phi', r', \theta', r', \phi', ', r',
                       ',pol', ',tor', ', pol', ', tor']
    for symb in remove_symbols:
        new_label = new_label.replace(symb, '')

    data.data.set(label=new_label)
    if no_nu:
        return data, label
    
    replace_nu = [
         r'\mathrm{\nu_e}', r'\mathrm{\nu_x}', r'\mathrm{\overline{\nu}_e}',
         r'\nu_\mathrm{e}', r'\nu_\mathrm{x}', r'\overline{\nu}_\mathrm{e}',
                 ]
    for nu  in replace_nu:
        new_label = new_label.replace(nu, r'\nu')
    data.data.set(label=new_label)
    return data, label  