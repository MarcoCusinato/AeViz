import os

def check_file_to_load(path):
    if path.endswith('.hdf5') or path.endswith('.h5') or path.endswith('.hdf'):
        return 'hdf5'
    elif not os.path.exists(path):
        return 'sim'

def convert_name(name):
    name_dict = {
      'rho': 'rho', 'rh': 'rho', 'rho_max': 'rho_max', 'rhom': 'rho_max',
      'mhd_energy': 'MHD_energy', 'mhd': 'MHD_energy', 'mhd_e': 'MHD_energy',
      'internal_energy': 'MHD_energy', 'internal': 'MHD_energy',
      'int_ene': 'MHD_energy',
      'radial_velocity': 'radial_velocity', 'vr': 'radial_velocity', 
      'vx': 'radial_velocity',
      'theta_velocity': 'theta_velocity', 'vtheta': 'theta_velocity',
      'vy': 'theta_velocity',
      'phi_velocity': 'phi_velocity', 'vphi': 'phi_velocity', 
      'vz': 'phi_velocity',
      'souundspeed': 'soundspeed', 'cs': 'soundspeed', 'csound': 'soundspeed',
      'omega': 'omega', 'w': 'omega',
      'pressure': 'gas_pressure', 'p': 'gas_pressure', 'pgas': 'gas_pressure',
      'temperature': 'temperature', 't': 'temperature', 'temp': 'temperature',
      'enthalpy': 'enthalpy', 'h': 'enthalpy', 'enth': 'enthalpy',
      'adiabatic_index': 'adiabatic_index', 'gamma': 'adiabatic_index',
      'lorentz_factor': 'lorentz', 'lorentz': 'lorentz',
      'gravitational_potential': 'gravitational_potential', 
      'phi': 'gravitational_potential', 'gpot': 'gravitational_potential',
      
    }

def return_index(hydro_dict, name):
    convert_dict = {
        'RHO': {'type': 'hydro/data', 'index': hydro_dict['hydro']['I_RH']},
        'ENE': {'type': 'hydro/data', 'index': hydro_dict['hydro']['I_EN']},
        'VX': {'type': 'hydro/data', 'index': hydro_dict['hydro']['I_VX']},
        'VY': {'type': 'hydro/data', 'index': hydro_dict['hydro']['I_VY']},
        'VZ': {'type': 'hydro/data', 'index': hydro_dict['hydro']['I_VZ']},
        'YE': {'type': 'hydro/data', 'index': hydro_dict['hydro']['I_YE']},
        'YN': {'type': 'hydro/data', 'index': hydro_dict['hydro']['I_YN']},
        'ERR': {'type': 'thd/data', 'index': hydro_dict['thd']['I_EOSERR']},
        'LRTZ': {'type': 'thd/data', 'index': hydro_dict['thd']['I_LRTZ']},
        'DENS': {'type': 'thd/data', 'index': hydro_dict['thd']['I_DENS']},
        'EINT': {'type': 'thd/data', 'index': hydro_dict['thd']['I_EINT']},
        'ENTH': {'type': 'thd/data', 'index': hydro_dict['thd']['I_ENTH']},
        'PELE': {'type': 'thd/data', 'index': hydro_dict['thd']['I_PELE']},
        'TELE': {'type': 'thd/data', 'index': hydro_dict['thd']['I_TELE']},
        'NELE': {'type': 'thd/data', 'index': hydro_dict['thd']['I_NELE']},
        'PION': {'type': 'thd/data', 'index': hydro_dict['thd']['I_PION']},
        'TION': {'type': 'thd/data', 'index': hydro_dict['thd']['I_TION']},
        'NION': {'type': 'thd/data', 'index': hydro_dict['thd']['I_NION']},
        'VELX': {'type': 'thd/data', 'index': hydro_dict['thd']['I_VELX']},
        'VELY': {'type': 'thd/data', 'index': hydro_dict['thd']['I_VELY']},
        'VELZ': {'type': 'thd/data', 'index': hydro_dict['thd']['I_VELZ']},
        'T': {'type': 'thd/data', 'index': hydro_dict['thd']['I_TMPR']},
        'ENTR': {'type': 'thd/data', 'index': hydro_dict['thd']['I_ENTR']},
        'GAMM': {'type': 'thd/data', 'index': hydro_dict['thd']['I_GAMM']},
        'HEAT': {'type': 'thd/data', 'index': hydro_dict['thd']['I_HEAT']},
        'DELP': {'type': 'thd/data', 'index': hydro_dict['thd']['I_DELP']},
        'JX': {'type': 'thd/data', 'index': hydro_dict['thd']['I_SMOMX']},
        'JY': {'type': 'thd/data', 'index': hydro_dict['thd']['I_SMOMY']},
        'JZ': {'type': 'thd/data', 'index': hydro_dict['thd']['I_SMOMZ']},
        'PGAS': {'type': 'thd/data', 'index': hydro_dict['thd']['I_PGAS']},
        'VSOUND': {'type': 'thd/data', 'index': hydro_dict['thd']['I_CSND']},
        'X_n': {'type': 'thd/data', 'index': hydro_dict['thd']['I_COMP'][0]},
        'X_p': {'type': 'thd/data', 'index': hydro_dict['thd']['I_COMP'][1]},
        'X_alpha': {'type': 'thd/data', 'index':
                    hydro_dict['thd']['I_COMP'][2]},
        'X_h': {'type': 'thd/data', 'index': hydro_dict['thd']['I_COMP'][3]},
        'Abar': {'type': 'thd/data', 'index': hydro_dict['thd']['I_COMP'][4]},
        'Zbar': {'type': 'thd/data', 'index': hydro_dict['thd']['I_COMP'][5]},
        'CPOT_e': {'type': 'thd/data', 'index': 
                   hydro_dict['thd']['I_CPOT'][0]},
        'CPOT_n': {'type': 'thd/data', 'index': 
                   hydro_dict['thd']['I_CPOT'][1]},
        'CPOT_p': {'type': 'thd/data', 'index': 
                   hydro_dict['thd']['I_CPOT'][2]},
        'CPOT_nu': {'type': 'thd/data', 'index': 
                    hydro_dict['thd']['I_CPOT'][3]},
        'BHEX': {'type': 'thd/data', 'index': hydro_dict['thd']['I_BHEX']},
        'NUEE': {'type': 'neutrinogrey/egrey', 'index': [0, 0] },
        'NUEX': {'type': 'neutrinogrey/egrey', 'index': [0, 1] },
        'NUEY': {'type': 'neutrinogrey/egrey', 'index': [0, 2] },
        'NUEZ': {'type': 'neutrinogrey/egrey', 'index': [0, 3] },
        'NUAE': {'type': 'neutrinogrey/egrey', 'index': [1, 0] },
        'NUAX': {'type': 'neutrinogrey/egrey', 'index': [1, 1] },
        'NUAY': {'type': 'neutrinogrey/egrey', 'index': [1, 2] },
        'NUAZ': {'type': 'neutrinogrey/egrey', 'index': [1, 3] },
        'NUXE': {'type': 'neutrinogrey/egrey', 'index': [2, 0] },
        'NUXX': {'type': 'neutrinogrey/egrey', 'index': [2, 1] },
        'NUXY': {'type': 'neutrinogrey/egrey', 'index': [2, 2] },
        'NUXZ': {'type': 'neutrinogrey/egrey', 'index': [2, 3] },
        'BX': {'type': 'mag_vol/data', 'index': 0},
        'BY': {'type': 'mag_vol/data','index': 1},
        'BZ': {'type': 'mag_vol/data', 'index': 2},
        }
    return convert_dict[name]