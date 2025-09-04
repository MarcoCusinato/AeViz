import os
import f90nml

def load_parfile(file_name, path_folder):
    parfile = f90nml.read(os.path.join(path_folder, file_name))
    return parfile

def get_indices_from_parfile(parfile):
    """
    Reads the simulation parfile and return a dictionary with the indices
    of the hydro and thd quantities.
    keywords: 
        first dictionary:
            'hydro' : hydrodynamics quantities
            'thd' : thermodynamics quantities
        second dictionary:

                                      HYDRO
            |----------------------------------------------------------|
            |'I_RH' : density           |   'I_EN' : conserved energy  |
            |'I_VX' : radial velocity   |   'I_YE' : electron fraction |
            |'I_VY' : polar velocity    |   'I_YN' :                   |
            |'I_VZ' : azimutal velocity |   'I_LRTZ' : lorentz factor  |
            ------------------------------------------------------------
                                       THD
    """
    indices = {'hydro':{}, 'thd':{}}
    if 'I_RH' in parfile["IINDICES"]:
        indices['hydro']['I_RH'] = parfile["IINDICES"]["I_RH"] - 1
    if 'I_EN' in parfile["IINDICES"]:
        indices['hydro']['I_EN'] = parfile["IINDICES"]["I_EN"] - 1
    if 'I_VX' in parfile["IINDICES"]:
        indices['hydro']['I_VX'] = parfile["IINDICES"]["I_VX"] - 1
    if 'I_VY' in parfile["IINDICES"]:
        indices['hydro']['I_VY'] = parfile["IINDICES"]["I_VY"] - 1
    if 'I_VZ' in parfile["IINDICES"]:
        indices['hydro']['I_VZ'] = parfile["IINDICES"]["I_VZ"] - 1
    if 'I_YE' in parfile["IINDICES"]:
        indices['hydro']['I_YE'] = parfile["IINDICES"]["I_YE"] - 1
    if 'I_YN' in parfile["IINDICES"]:
        indices['hydro']['I_YN'] = parfile["IINDICES"]["I_YN"] - 1
    if 'I_EOSERR' in parfile["IINDICES"]:
        indices['thd']['I_EOSERR'] = parfile["IINDICES"]["I_EOSERR"] - 1
    if 'I_LRTZ' in parfile["IINDICES"]:
        indices['thd']['I_LRTZ'] = parfile["IINDICES"]["I_LRTZ"] - 1
    if 'I_DENS' in parfile["IINDICES"]:
        indices['thd']['I_DENS'] = parfile["IINDICES"]["I_DENS"] - 1    
    if 'I_EINT' in parfile["IINDICES"]:
        indices['thd']['I_EINT'] = parfile["IINDICES"]["I_EINT"] - 1
    if 'I_ENTH' in parfile["IINDICES"]:
        indices['thd']['I_ENTH'] = parfile["IINDICES"]["I_ENTH"] - 1
    if 'I_PELE' in parfile["IINDICES"]:
        indices['thd']['I_PELE'] = parfile["IINDICES"]["I_PELE"] - 1    
    if 'I_TELE' in parfile["IINDICES"]:
        indices['thd']['I_TELE'] = parfile["IINDICES"]["I_TELE"] - 1
    if 'I_NELE' in parfile["IINDICES"]:
        indices['thd']['I_NELE'] = parfile["IINDICES"]["I_NELE"] - 1
    if 'I_PION' in parfile["IINDICES"]:
        indices['thd']['I_PION'] = parfile["IINDICES"]["I_PION"] - 1
    if 'I_TION' in parfile["IINDICES"]:
        indices['thd']['I_TION'] = parfile["IINDICES"]["I_TION"] - 1    
    if 'I_NION' in parfile["IINDICES"]:
        indices['thd']['I_NION'] = parfile["IINDICES"]["I_NION"] - 1
    if 'I_VELX' in parfile["IINDICES"]:
        indices['thd']['I_VELX'] = parfile["IINDICES"]["I_VELX"] - 1
    if 'I_VELY' in parfile["IINDICES"]:
        indices['thd']['I_VELY'] = parfile["IINDICES"]["I_VELY"] - 1    
    if 'I_VELZ' in parfile["IINDICES"]:
        indices['thd']['I_VELZ'] = parfile["IINDICES"]["I_VELZ"] - 1
    if 'I_TMPR' in parfile["IINDICES"]:
        indices['thd']['I_TMPR'] = parfile["IINDICES"]["I_TMPR"] - 1
    if 'I_ENTR' in parfile["IINDICES"]:
        indices['thd']['I_ENTR'] = parfile["IINDICES"]["I_ENTR"] - 1
    if 'I_GAMM' in parfile["IINDICES"]:
        indices['thd']['I_GAMM'] = parfile["IINDICES"]["I_GAMM"] - 1    
    if 'I_HEAT' in parfile["IINDICES"]:
        indices['thd']['I_HEAT'] = parfile["IINDICES"]["I_HEAT"] - 1
    if 'I_DELP' in parfile["IINDICES"]:
        indices['thd']['I_DELP'] = parfile["IINDICES"]["I_DELP"] - 1
    if 'I_SMOMX' in parfile["IINDICES"]:
        indices['thd']['I_SMOMX'] = parfile["IINDICES"]["I_SMOMX"] - 1    
    if 'I_SMOMY' in parfile["IINDICES"]:
        indices['thd']['I_SMOMY'] = parfile["IINDICES"]["I_SMOMY"] - 1
    if 'I_SMOMZ' in parfile["IINDICES"]:
        indices['thd']['I_SMOMZ'] = parfile["IINDICES"]["I_SMOMZ"] - 1
    if 'I_PGAS' in parfile["IINDICES"]:
        indices['thd']['I_PGAS'] = parfile["IINDICES"]["I_PGAS"] - 1
    if 'I_CSND' in parfile["IINDICES"]:
        indices['thd']['I_CSND'] = parfile["IINDICES"]["I_CSND"] - 1    
    if 'I_COMP' in parfile["IINDICES"]:
        list_comp_ind = parfile["IINDICES"]["I_COMP"]
        list_comp_ind = [x -1 for x in list_comp_ind]
        indices['thd']['I_COMP'] = list_comp_ind
    if 'I_CPOT' in parfile["IINDICES"]:
        list_cpot_ind = parfile["IINDICES"]["I_CPOT"]
        list_cpot_ind = [x -1 for x in list_cpot_ind]
        indices['thd']['I_CPOT'] = list_cpot_ind
    if 'I_BHEX' in parfile["IINDICES"]:
        indices['thd']['I_BHEX'] = parfile["IINDICES"]["I_BHEX"] - 1    
    for  k1 in indices['hydro']:
        if type(indices['hydro'][k1]) == list:
            for i in range(0,len(indices['hydro'][k1])):
                if indices['hydro'][k1][i] < 0:
                    indices['hydro'][k1][i] = None
        elif indices['hydro'][k1] < 0:
            indices['hydro'][k1] = None
    for  k1 in indices['thd']:
        if type(indices['thd'][k1]) == list:
            for i in range(0, len(indices['thd'][k1])):
                if indices['thd'][k1][i] < 0:
                    indices['thd'][k1][i] = None
        elif indices['thd'][k1] < 0:
            indices['thd'][k1] = None
    return indices

def get_stencils(parfile):
    """
    Returns the number of ghost cells per local output
    """
    return parfile["GRIDPARS"]["STENCIL"]

def get_simulation_info(parfile):
    """
    Returns the geometry of the simulation and a dictionary with the
    list of quantities in the simulations
    """
    varspars = parfile["VARSPARS"]
    geometry = varspars["GEOMETRY"]
    dim = varspars["EVOLVE_X"] + varspars["EVOLVE_Y"] + varspars["EVOLVE_Z"]
    lapse_form = int(parfile["GRAVPARS"]["LAPSE_FORM"])
    
    qts = {
        key.casefold(): varspars[key] for key in varspars.keys() if key not in
        ['geometry', 'evolve_x', 'evolve_y', 'evolve_z']
    }
    return geometry, dim, parfile["PHYSSYST"]["RELATIVISTIC"], qts, lapse_form

def get_initial_parameters(path_folder):
    initial_parameters = {"omgadd": 0.0,
                          "omgmult": 0.0,
                          "bfact": 0.0,
                          "btorfact": 0.0,
                          "b0": 0.0,
                          "bt": 0.0
                          }
    start_pars = f90nml.read(os.path.join(path_folder, 'start.pars'))
    if start_pars["PHYSSYST"]["RELATIVISTIC"]:
        initial_parameters['gravity'] = 'Pseudo-relativistic'
    else:
        initial_parameters['gravity'] = 'Newtonian'
    initial_parameters["gravitational_potential"] = start_pars["GRAVPARS"]["MDPOT"]
    initial_parameters["lapse_function"] = start_pars["GRAVPARS"]["LAPSE_FORM"]
    initial_parameters["NS_EOS"] = (os.path.basename(start_pars["SHENEOSPARS"]["SHEN_TBFILE"]))
    initial_parameters["NS_EOS"] = initial_parameters["NS_EOS"].split('.')[:-1]
    if type(initial_parameters["NS_EOS"]) == list:
        if len(initial_parameters["NS_EOS"]) == 1:
            initial_parameters["NS_EOS"] = initial_parameters["NS_EOS"][0]
        else:
            tmp = initial_parameters["NS_EOS"][0]
            for i in initial_parameters["NS_EOS"][1:]:
                tmp += '.' + i
            initial_parameters["NS_EOS"] = tmp
    files = os.listdir(path_folder)
    files = [x for x in files if not x == 'start.pars' and not x == '.run']
    files.sort()
    Found_heger = False
    for file in files:
        namelist = f90nml.read(os.path.join(path_folder, file))
        if 'Hegerpars' in namelist:
            Found_heger = True
            initial_parameters["Heger_model"] = namelist["Hegerpars"]["Heger_model"]
            if "omgadd" in namelist["Hegerpars"]:
                initial_parameters["omgadd"] = namelist["Hegerpars"]["omgadd"]
            if "omgmult" in namelist["Hegerpars"]:
                initial_parameters["omgmult"] = namelist["Hegerpars"]["omgmult"]
            if "bfield" in namelist["Hegerpars"]:
                if namelist["Hegerpars"]["bfield"]:
                    if "bfact" in namelist["Hegerpars"]:
                        initial_parameters["bfact"] = namelist["Hegerpars"]["bfact"]
                    if "btorfact" in namelist["Hegerpars"]:
                        initial_parameters["btorfact"] = namelist["Hegerpars"]["btorfact"]
        if "axivecpotpars" in namelist:
            if "b0" in namelist["axivecpotpars"]:
                initial_parameters["b0"] = namelist["axivecpotpars"]["b0"]
            if "bt" in namelist["axivecpotpars"]:
                initial_parameters["bt"] = namelist["axivecpotpars"]["bt"]
        if Found_heger:
            break
    if not Found_heger:
        run_path = os.path.join(path_folder, '.run')
        run_files = os.listdir(run_path)
        run_files.sort()
        for file in run_files:
            namelist = f90nml.read(os.path.join(run_path, file))
            if 'Hegerpars' in namelist and namelist["Hegerpars"]["Heger_model"]:
                Found_heger = True
                initial_parameters["Heger_model"] = namelist["Hegerpars"]["Heger_model"]
                if "omgadd" in namelist["Hegerpars"]:
                    initial_parameters["omgadd"] = namelist["Hegerpars"]["omgadd"]
                if "omgmult" in namelist["Hegerpars"]:
                    initial_parameters["omgmult"] = namelist["Hegerpars"]["omgmult"]
                if "bfield" in namelist["Hegerpars"]:
                    if namelist["Hegerpars"]["bfield"]:
                        if "bfact" in namelist["Hegerpars"]:
                            initial_parameters["bfact"] = namelist["Hegerpars"]["bfact"]
                        if "btorfact" in namelist["Hegerpars"]:
                            initial_parameters["btorfact"] = namelist["Hegerpars"]["btorfact"]
            if "axivecpotpars" in namelist:
                if "b0" in namelist["axivecpotpars"]:
                    initial_parameters["b0"] = namelist["axivecpotpars"]["b0"]
                if "bt" in namelist["axivecpotpars"]:
                    initial_parameters["bt"] = namelist["axivecpotpars"]["bt"]
            if Found_heger:
                break
    if not Found_heger:
        initial_parameters["omgadd"] = None
        initial_parameters["omgmult"] = None
        initial_parameters["Heger_model"] = None
        initial_parameters["bfact"] = None
        initial_parameters["btorfact"] = None
        initial_parameters["b0"] = None
        initial_parameters["bt"] = None
    return initial_parameters
