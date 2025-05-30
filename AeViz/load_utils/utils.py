import os
import numpy as np
from AeViz.units.units import units

u = units()

def return_neutrino_flux(sim, name, file, **kwargs):
    comp = name.split('_')[-1]
    if comp == 'e':
        out = sim.neutrino_energy_density_grey(file)
    else: 
        out = sim.neutrino_momenta_grey(file)
    if 'nue' in name:
        if sim.dim == 1 or comp == 'e':
            out = out[..., 0]
        elif comp == 'x':
            out = out[..., 0, 0]
        elif comp == 'y':
            out = out[..., 0, 1]
        elif comp == 'z':
            out = out[..., 0, 2]
    elif 'nua' in name:
        if sim.dim == 1 or comp == 'e':
            out = out[..., 1]
        elif comp == 'x':
            out = out[..., 1, 0]
        elif comp == 'y':
            out = out[..., 1, 1]
        elif comp == 'z':
            out = out[..., 1, 2]
    elif 'nux' in name:
        if sim.dim == 1 or comp == 'e':
            out = out[..., 2]
        elif comp == 'x':
            out = out[..., 2, 0]
        elif comp == 'y':
            out = out[..., 2, 1]
        elif comp == 'z':
            out = out[..., 2, 2]
    return out

def return_neutrino_mean_ene(sim, name, file, **kwargs):
    if 'nue' in name:
        out = sim.neutrino_mean_energy(file)[..., 0]
    elif 'nua' in name:
        out = sim.neutrino_mean_energy(file)[..., 1]
    elif 'nux' in name:
        out = sim.neutrino_mean_energy(file)[..., 2]
    return out

def return_integrated_neutrinos(sim, name, **kwargs):
    if 'ene' in name:
        out = sim.global_neutrino_mean_energies(**kwargs)
    elif 'lum' in name:
        out = sim.global_neutrino_luminosity(**kwargs)[:, :4]
    if 'all' in name:
        return out
    elif 'nue' in name:
        return out[:, :2]
    elif 'nua' in name:
        return np.stack((out[:, 0], out[:, 2]), axis=1)
    elif 'nux' in name:
        return np.stack((out[:, 0], out[:, 3]), axis=1)

def return_angular_momentum(sim, name, **kwargs):
    data = sim.PNS_angular_mom(**kwargs)
    if 'all' in name:
        return data
    elif 'Lx' in name:
        return [data[0], data[1]]
    elif 'Ly' in name:
        return [data[0], data[2]]
    elif 'Lz' in name:
        return [data[0], data[3]]
    elif 'Ltot' in name:
        return [data[0], data[4]]
    
def return_PNS_kick(sim, name, **kwargs):
    data = list(sim.PNS_kick_velocity(**kwargs))
    data[1] = u.convert_to_km(data[1])
    data[2] = u.convert_to_km(data[2])
    data[3] = u.convert_to_km(data[3])
    if 'all' in name:
        return data
    elif 'nu' in name:
        return [data[0], data[3]]
    elif 'modulus' in name:
        return [data[0], data[1]]
    elif 'hydro' in name:
        return [data[0], data[2]]
    
def check_file_to_load(path):
    if path.endswith('.hdf5') or path.endswith('.h5') or path.endswith('.hdf'):
        return 'hdf5'
    elif not os.path.exists(path) or os.path.isdir(path):
        return 'sim'

def return_index(hydro_dict, name):
    convert_dict = {
        'rho': {'type': 'hydro/data', 'index': hydro_dict['hydro']['I_RH']},
        'MHD_energy': {'type': 'hydro/data', 
                       'index': hydro_dict['hydro']['I_EN']},
        'VX': {'type': 'hydro/data', 'index': hydro_dict['hydro']['I_VX']},
        'VY': {'type': 'hydro/data', 'index': hydro_dict['hydro']['I_VY']},
        'VZ': {'type': 'hydro/data', 'index': hydro_dict['hydro']['I_VZ']},
        'Ye': {'type': 'hydro/data', 'index': hydro_dict['hydro']['I_YE']},
        'YN': {'type': 'hydro/data', 'index': hydro_dict['hydro']['I_YN']},
        'error': {'type': 'thd/data', 'index': hydro_dict['thd']['I_EOSERR']},
        'lorentz': {'type': 'thd/data', 'index': hydro_dict['thd']['I_LRTZ']},
        'DENS': {'type': 'thd/data', 'index': hydro_dict['thd']['I_DENS']},
        'internal_energy': {'type': 'thd/data',
                            'index': hydro_dict['thd']['I_EINT']},
        'enthalphy': {'type': 'thd/data',
                      'index': hydro_dict['thd']['I_ENTH']},
        'PELE': {'type': 'thd/data', 'index': hydro_dict['thd']['I_PELE']},
        'TELE': {'type': 'thd/data', 'index': hydro_dict['thd']['I_TELE']},
        'NELE': {'type': 'thd/data', 'index': hydro_dict['thd']['I_NELE']},
        'PION': {'type': 'thd/data', 'index': hydro_dict['thd']['I_PION']},
        'TION': {'type': 'thd/data', 'index': hydro_dict['thd']['I_TION']},
        'NION': {'type': 'thd/data', 'index': hydro_dict['thd']['I_NION']},
        'radial_velocity': {'type': 'thd/data',
                            'index': hydro_dict['thd']['I_VELX']},
        'theta_velocity': {'type': 'thd/data',
                           'index': hydro_dict['thd']['I_VELY']},
        'phi_velocity': {'type': 'thd/data',
                         'index': hydro_dict['thd']['I_VELZ']},
        'temperature': {'type': 'thd/data',
                        'index': hydro_dict['thd']['I_TMPR']},
        'entropy': {'type': 'thd/data', 'index': hydro_dict['thd']['I_ENTR']},
        'adiabatic_index': {'type': 'thd/data',
                            'index': hydro_dict['thd']['I_GAMM']},
        'HEAT': {'type': 'thd/data', 'index': hydro_dict['thd']['I_HEAT']},
        'DELP': {'type': 'thd/data', 'index': hydro_dict['thd']['I_DELP']},
        'JX': {'type': 'thd/data', 'index': hydro_dict['thd']['I_SMOMX']},
        'JY': {'type': 'thd/data', 'index': hydro_dict['thd']['I_SMOMY']},
        'JZ': {'type': 'thd/data', 'index': hydro_dict['thd']['I_SMOMZ']},
        'gas_pressure': {'type': 'thd/data',
                         'index': hydro_dict['thd']['I_PGAS']},
        'soundspeed': {'type': 'thd/data',
                       'index': hydro_dict['thd']['I_CSND']},
        'neutron_fraction': {'type': 'thd/data', 'index': hydro_dict['thd']['I_COMP'][0]},
        'proton_fraction': {'type': 'thd/data', 'index': hydro_dict['thd']['I_COMP'][1]},
        'alpha_fraction': {'type': 'thd/data', 'index':
                    hydro_dict['thd']['I_COMP'][2]},
        'heavy_fraction': {'type': 'thd/data', 'index': hydro_dict['thd']['I_COMP'][3]},
        'Abar': {'type': 'thd/data', 'index': hydro_dict['thd']['I_COMP'][4]},
        'Zbar': {'type': 'thd/data', 'index': hydro_dict['thd']['I_COMP'][5]},
        'electron_chemical_potential': {'type': 'thd/data', 'index': 
                   hydro_dict['thd']['I_CPOT'][0]},
        'neutron_chemical_potential': {'type': 'thd/data', 'index': 
                   hydro_dict['thd']['I_CPOT'][1]},
        'proton_chemical_potential': {'type': 'thd/data', 'index': 
                   hydro_dict['thd']['I_CPOT'][2]},
        'neutrino_chemical_potential': {'type': 'thd/data', 'index': 
                    hydro_dict['thd']['I_CPOT'][3]},
        'BHEX': {'type': 'thd/data', 'index': hydro_dict['thd']['I_BHEX']},
        'nue_moment_e': {'type': 'neutrinogrey/egrey', 'index': [0, 0] },
        'nue_moment_x': {'type': 'neutrinogrey/egrey', 'index': [0, 1] },
        'nue_moment_y': {'type': 'neutrinogrey/egrey', 'index': [0, 2] },
        'nue_moment_z': {'type': 'neutrinogrey/egrey', 'index': [0, 3] },
        'nua_moment_e': {'type': 'neutrinogrey/egrey', 'index': [1, 0] },
        'nua_moment_x': {'type': 'neutrinogrey/egrey', 'index': [1, 1] },
        'nua_moment_y': {'type': 'neutrinogrey/egrey', 'index': [1, 2] },
        'nua_moment_z': {'type': 'neutrinogrey/egrey', 'index': [1, 3] },
        'nux_moment_e': {'type': 'neutrinogrey/egrey', 'index': [2, 0] },
        'nux_moment_x': {'type': 'neutrinogrey/egrey', 'index': [2, 1] },
        'nux_moment_y': {'type': 'neutrinogrey/egrey', 'index': [2, 2] },
        'nux_moment_z': {'type': 'neutrinogrey/egrey', 'index': [2, 3] },
        'BX': {'type': 'mag_vol/data', 'index': 0},
        'BY': {'type': 'mag_vol/data','index': 1},
        'BZ': {'type': 'mag_vol/data', 'index': 2},
        }
    return convert_dict[name]