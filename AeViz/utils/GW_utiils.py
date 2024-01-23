import numpy as np
from units.units import units
from utils.math_utils import IDL_derivative
from scipy.signal import stft

u = units()

def GW_strain(sim_dim, column_change, data, index):
    assert sim_dim in [1, 2, 3], "Simulation MUST be 1, 2 or 3D."
    if sim_dim == 1:
        return GW_strain_1D(data)
    elif sim_dim == 2:
        return correct_zero(2, GW_strain_2D(data), index)
    else:
        GWs = GW_strain_3D(data)
        if column_change is not None:
            GWs[:column_change, 1] = GW_strain_2D(data)[:column_change, 1]
            GWs[:column_change,2:] = 0
        return correct_zero(3, GWs, index)

def GW_strain_1D(data):
    print("No GW for you :'(")
    return None

def GW_strain_2D(data):
    """
        GWs amplitudes calculated as the first partial time derivative
        of NE_220 and not with the second partial time derivative of ME_220.
        Moreover we only consider matter contribution to the amplitude, being
        the main contribution to it.
    """
    const = -0.125 *  np.sqrt(15/np.pi)
    return np.stack((data[:, 2], const * IDL_derivative(data[:,2], data[:,5])),
              axis = 1)

def GW_strain_3D(data):
    const = 1 #8 * np.sqrt(np.pi / 15)
    hD_pl_p = 2. * ( data[:,9] - data[:,13] )
    hD_pl_e = 2. * ( data[:,17] - data[:,13] )
    hD_cr_p = 2. * ( data[:,10] + data[:,12] )
    hD_cr_e = 2. * ( - data[:,14] - data[:,16] )

    hD_pl_p = const * IDL_derivative( data[:,2], hD_pl_p )
    hD_cr_p = const * IDL_derivative( data[:,2], hD_cr_p )
    hD_pl_e = const * IDL_derivative( data[:,2], hD_pl_e )
    hD_cr_e = const * IDL_derivative( data[:,2], hD_cr_e )

    return np.stack((data[:,2], hD_pl_e, hD_pl_p, hD_cr_e, hD_cr_p), axis = -1)

def correct_zero(sim_dim, GWs, index):
    if sim_dim == 1:
        pass
    else:
        if index is None:
            return GWs
        for i in range(1, GWs.shape[1]):
            GWs[:, i] -= GWs[:index, i].mean()
        return GWs
def GW_spectrogram(sim_dim, GWs, window_size):
    assert sim_dim in [1, 2, 3], "Simulation MUST be 1, 2 or 3D."
    if sim_dim == 1:
        return GW_spectrogram_1D(GWs)
    elif sim_dim == 2:
        return GW_spectrogram_2D(GWs, window_size)
    else:
        return GW_spectrogram_3D(GWs, window_size)

def GW_spectrogram_1D(GWs):
    print("And also no spectrogram for you :'(\nごめんなさい")
    return None

def GW_spectrogram_2D(GWs, window_size):
    """
        
    """
    fs = (GWs[1, 0] - GWs[0, 0]) ** (-1)
    frequency, time, Zxx = stft(GWs[:,1], fs = fs, nperseg=window_size)
    time = time / time[-1] * GWs[-1, 0]
    return time, frequency, np.abs(Zxx)

def GW_spectrogram_3D(GWs, window_size):
    fs = (GWs[1, 0] - GWs[0, 0]) ** (-1)
    freq_pl_e, time_pl_e, Zxx_pl_e = stft(GWs[:,1], fs = fs, nperseg=window_size)
    time_pl_e = time_pl_e / time_pl_e[-1] * GWs[-1, 0]
    freq_pl_p, time_pl_p, Zxx_pl_p = stft(GWs[:,2], fs = fs, nperseg=window_size)
    time_pl_p = time_pl_p / time_pl_p[-1] * GWs[-1, 0]
    freq_cr_e, time_cr_e, Zxx_cr_e = stft(GWs[:,3], fs = fs, nperseg=window_size)
    time_cr_e = time_cr_e / time_cr_e[-1] * GWs[-1, 0]
    freq_cr_p, time_cr_p, Zxx_cr_p = stft(GWs[:,4], fs = fs, nperseg=window_size)
    time_cr_p = time_cr_p / time_cr_p[-1] * GWs[-1, 0]
    
    return np.stack((time_pl_e, time_pl_p, time_cr_e, time_cr_p), axis=-1), \
           np.stack((freq_pl_e, freq_pl_p, freq_cr_e, freq_cr_p), axis=-1), \
           np.stack((np.abs(Zxx_pl_e), np.abs(Zxx_pl_p), np.abs(Zxx_cr_e),
                     np.abs(Zxx_cr_p)), axis=-1)

def calculate_AE220(simulation):
    """
    Calculates the AE220 from density and velocities for a
    2D simulation.
    ONLY 2D
    Returns
        time: arry of time step
        AE220: len(radius), len(time) array
    """
    radius = simulation.cell.radius(simulation.ghost)
    dV = -simulation.cell.dVolume_integration(simulation.ghost)
    costheta = np.cos(simulation.cell.theta(simulation.ghost))
    file_list = simulation.file_list_hdf()
    inner_radius, time, gh_inner = simulation.get_innercore_radius(innercore_radius=True, ret_time=True,
                                                    ghost_cells=True)
    convective_radius, gh_conv = simulation.get_convective_radius(convective_radius=True,
                                                            ghost_cells=True)
    NE220_full = np.zeros((len(radius), len(time)))
    NE220_inner = np.zeros((len(radius), len(time)))
    NE220_convection = np.zeros((len(radius), len(time)))
    NE220_outer = np.zeros((len(radius), len(time)))
    radius = radius[None, ...]
    costheta = costheta[..., None]
    for index in range(len(file_list)):
        # Get all the quantities needed for the calculation
        try:
            data_h5 = simulation.open_h5(file_list[index])
        except:
            NE220_full[:, index] = NE220_full[:, index - 1]
            NE220_inner[:, index] = NE220_inner[:, index - 1]
            NE220_convection[:, index] = NE220_convection[:, index - 1]
            NE220_outer[:, index] = NE220_outer[:, index - 1]
            continue
        rho = simulation.rho(data_h5)
        vR = simulation.radial_velocity(data_h5)
        vT = simulation.theta_velocity(data_h5)
        simulation.close_h5(data_h5)
        # Calculate the NE220
        # NE220 of the full space
        NE220 = dV * radius * rho *  ( vR * \
            ( 3 * costheta ** 2 - 1 ) - 3 * vT * costheta * np.sqrt( 1 - costheta ** 2 ) )
        # Create masks
        mask_inner = radius <= simulation.ghost.remove_ghost_cells_radii(convective_radius[..., index],
                                                                        simulation.dim, **gh_conv)[..., None]
        mask_convection = (radius <= simulation.ghost.remove_ghost_cells_radii(inner_radius[..., index],
                                                                            simulation.dim, **gh_inner)[..., None] + 2e6) & \
                            ( np.logical_not(mask_inner) )
        mask_outer = np.logical_not(mask_inner + mask_convection)
        # Calculate the NE220s
        NE220_full[:, index] = np. sum( NE220, axis = 0 )
        NE220_inner[:, index] = np. sum( NE220 * mask_inner, axis = 0 )
        NE220_convection[:, index] = np. sum( NE220 * mask_convection, axis = 0 )
        NE220_outer[:, index] = np. sum( NE220 * mask_outer, axis = 0 )
    
    const = ( u.G * 16 * np.pi ** 0.5 / (np.sqrt( 15 ) * u.speed_light ** 4 ) )

    return time, IDL_derivative(time, NE220_full * const), calculate_strain(NE220_full, time), \
            calculate_strain(NE220_inner, time), calculate_strain(NE220_convection, time), \
            calculate_strain(NE220_outer, time)
            

def calculate_strain(NE220, time):
    """
    Calculates the strain from the NE220
    """
    # Put again the correct units
    NE220 *= ( u.G * 16 * np.pi ** 0.5 / (np.sqrt( 15 ) * u.speed_light ** 4 ) )
    # Calculate the strain
    strain = IDL_derivative(time, NE220)
    return strain.sum(axis=0)
