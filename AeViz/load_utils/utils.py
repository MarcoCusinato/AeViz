import os

def check_file_to_load(path):
    if path.endswith('.hdf5') or path.endswith('.h5') or path.endswith('.hdf'):
        return 'hdf5'
    elif not os.path.exists(path):
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