from AeViz.units import (u, apply_monkey_patch, remove_monkey_patch,
                         aerray, aeseries)
import AeViz.utils.snewpy
from AeViz.utils.snewpy import detectors_list
import numpy as np
import os
from typing import Literal, get_args
try:
    from snewpy import snowglobes
except:
    remove_monkey_patch()
    from snewpy import snowglobes
    apply_monkey_patch()

"""
Interface between AeViz and SNEWPY and hence SNOwGLoBES. Two caveats:
1- It uses my version of SNEWPY, which contains the reader class for
   the Aenus_files
2- Supposes my naming convention of simulations: 
            `progenitor-EOS-rotation-Bfield-otherstuff'
"""

av_transformations = Literal['NoTransformation',
                      'AdiabaticMSW_NMO',
                      'AdiabaticMSW_IMO',
                      'NonAdiabaticMSWH_NMO',
                      'NonAdiabaticMSWH_IMO',
                      'TwoFlavorDecoherence',
                      'ThreeFlavorDecoherence',
                      'NeutrinoDecay_NMO',
                      'NeutrinoDecay_IMO',
                      'QuantumDecoherence_NMO',
                      'QuantumDecoherence_IMO']

def get_snewpy_metadata(simulation, los:Literal['avg', 'eq', 'pol']='avg'):
    ## Create the sn dictionary for SNEWPY and extract the path of the model
    sim_name = simulation.simulation_name.split('-')
    progenitor = sim_name[0]
    eos = sim_name[1]
    rotrate = sim_name[2]
    ostuff = ''
    if (sim_name[3] == 'orig' or sim_name[3].endswith('B')):
        bfield = sim_name[3]
        if len(sim_name) > 4:
            ostuff = '-' + '-'.join(sim_name[4:])
    else:
        bfield = '-'.join(sim_name[3:5])
        if len(sim_name) > 5:
            ostuff = '-' + '-'.join(sim_name[5:])
    snmodel_dict = {'los':los,
                    'Bfield': bfield,
                    'rotation': rotrate,
                    'eos':eos,
                    'suffix': ostuff}
    model_path = os.path.join(simulation.storage_path, 'snewpy_files',
                              progenitor)
    return snmodel_dict, model_path

def generate_event_rate_series(simulation, los:Literal['avg', 'eq', 'pol']='avg',
                               distance=(10*u.kpc), dt=(1*u.ms), 
                               transformation:av_transformations='NoTransformation',
                               detector:detectors_list='superk'):
    assert transformation in get_args(av_transformations), \
        f"Transformation not recognized: choose from {get_args(av_transformations)}"
    assert detector in get_args(detectors_list), \
        f"Detector not recognized: choose from {get_args(detectors_list)}"
    snmodel_dict, model_path = get_snewpy_metadata(simulation, los)
    ## get the detector
    det = getattr(AeViz.utils.snewpy, detector)
    ## Generate a timeseries for the fluence with SNOwGLoBES / SNEWPY
    if isinstance(dt, aerray):
        dt = dt.to(u.s).value
    if isinstance(distance, aerray):
        distance = distance.to(u.kpc).value
    ##Convert time to astropy unit
    remove_monkey_patch()
    dt_bck = dt
    dt = u.s * dt
    print('Generating time series ...')
    ts = snowglobes.generate_time_series(model_path,
                                         'Aenus_models',
                                         transformation,
                                         distance,
                                         None,
                                         deltat=dt,
                                         snmodel_dict=snmodel_dict)
    print("Simulating detector effects with SNOwGLoBES ...")
    snow_res = snowglobes.simulate(None, ts, detector_input=det.det_type,
                                   save_file=False)
    print("Collating results ...")
    evs = snowglobes.collate(snow_res, True)
    apply_monkey_patch()
    # we want the total event rate per single bin
    nbins = (len(evs.keys()) - 1) // 2
    ev_rate = np.zeros(nbins)
    for i in range(nbins):
        ## We just want the smeared events
        key_name = f'Collated_aenus_model_{i}_{det.det_type}_events_smeared_weighted.dat'
        ev_rate[i] = np.sum(evs[key_name]['data'][1:, :])
    
    ## Account for dection mass
    ev_rate *= det.eff_fact * u.dimensionless_unscaled
    ## Find the initial time
    tm = np.loadtxt(os.path.join(simulation.storage_path, 'snewpy_files',
                                 f'{simulation.simulation_name}_nue_{los}.txt'))[0, 0]
    time = (dt_bck * np.arange(nbins) + tm) * u.s
    time.set(name='time', label=r'$t-t_\mathrm{b}$',
             limits=[-0.005, time[-1].value], log=False)
    ev_rate = ev_rate / (dt_bck * u.s)
    ev_rate.set(f'Event_rate {det.name}', r'Event rate', log=False,
                limits=[0, ev_rate.value.max() * 1.1])
    return aeseries(ev_rate, time=time)

def total_events_number(simulation, los:Literal['avg', 'eq', 'pol']='avg',
                        distance=(10*u.kpc), time_range=[None, None], 
                        transformation:av_transformations='NoTransformation',
                        detector:detectors_list='superk'):
    assert transformation in get_args(av_transformations), \
        f"Transformation not recognized: choose from {get_args(av_transformations)}"
    assert detector in get_args(detectors_list), \
        f"Detector not recognized: choose from {get_args(detectors_list)}"
    snmodel_dict, model_path = get_snewpy_metadata(simulation, los)
    ## get the detector
    det = getattr(AeViz.utils.snewpy, detector)
    if isinstance(distance, aerray):
        distance = distance.to(u.kpc).value
    ## Find the time
    tm = np.loadtxt(os.path.join(simulation.storage_path, 'snewpy_files',
                                 f'{simulation.simulation_name}_nue_{los}.txt'))[:, 0]
    if time_range[0] is None:
        tstart = tm[0]
    else:
        if isinstance(time_range[0], aerray):
            tstart = time_range[0].to(u.s).value
        else:    
            tstart = time_range[0]
        if tstart < tm[0]:
            tstart = tm[0]
        
    if time_range[1] is None:
        tend = tm[-1]
    else:
        if isinstance(time_range[1], aerray):
            tend = time_range[1].to(u.s).value
        else:    
            tend = time_range[1]
        if tend > tm[-1]:
            tend = tm[-1]
    
    remove_monkey_patch()
    tstart = tstart * u.s
    tend = tend * u.s
    print('Generating fluence ...')
    fluence = snowglobes.generate_fluence(model_path,
                                          'Aenus_models',
                                          transformation,
                                          distance,
                                          None,
                                          tstart=tstart,
                                          tend=tend,
                                          snmodel_dict=snmodel_dict)
    print("Simulating detector effects with SNOwGLoBES ...")
    snow_res = snowglobes.simulate(None, fluence, detector_input=det.det_type,
                                   save_file=False)
    print("Collating results ...")
    evs = snowglobes.collate(snow_res, True)
    apply_monkey_patch()
    key_name = f'Collated_aenus_model_{det.det_type}_events_smeared_weighted.dat'
    print(f"Number of events for {det.name} at {distance} [kpc]")
    total_events = 0
    for i, channel in enumerate(evs[key_name]['header'].split()):
        if i == 0:
            continue
        n_events = sum(evs[key_name]['data'][i]) * det.eff_fact
        total_events += n_events
        print(f"{channel:10}: {n_events:.3f} events")
    print(f'Total events: {total_events}')
    return total_events, det



    
    


    
    

    
    
