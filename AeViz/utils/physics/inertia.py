from AeViz.utils.utils import (check_existence, progressBar, checkpoints)
from AeViz.utils.files.file_utils import save_hdf, create_series
from AeViz.grid.grid import grid
import os, h5py
import numpy as np
from AeViz.units.aerray import aerray
from AeViz.units import u

def compute_inertia(simulation, inertia_density, radius, radius_limit, gcell):
    if simulation.dim == 1:
        mask = (radius <= radius_limit)
    else:
        mask = (radius <= \
                simulation.ghost.remove_ghost_cells_radii(radius_limit,
                                                          simulation.dim,
                                                          **gcell)[..., None])
    return np.sum(inertia_density[mask])

def calculate_moment_inertia(simulation, save_checkpoints=True):
    if check_existence(simulation, 'inertia_moment.h5'):
        time, iso14, iso13, iso12, iso11, iso10, iso9, iso8, \
            processed_hdf = \
            read_inertia_moment(simulation)
        if processed_hdf[-1].decode("utf-8") == simulation.hdf_file_list[-1]:
            return create_series(time, iso14, iso13, iso12, iso11, iso10, iso9,
                                 iso8)
        else:
            start_point = len(processed_hdf)
            processed_hdf = [ff.decode("utf-8") for ff in processed_hdf]
            print('Checkpoint found for the moment of inertia file, starting' \
                  ' from checkpoint.\nPlease wait...')
    else:
        start_point = 0
        processed_hdf = []
        print('No checkpoint found for the moment of inertia file, starting' \
              ' from the beginning.\nPlease wait...')
    if (checkpoints[simulation.dim] == False) or (not save_checkpoints):
        checkpoint = len(simulation.hdf_file_list)
    else:
        checkpoint = checkpoints[simulation.dim]
    
    ## Set up the points indices
    findex = start_point
    check_index = 0
    progress_index = 0
    total_points = len(simulation.hdf_file_list) - start_point

    ## Radii
    r1e08, giso = simulation.isodensities_lines(rad='full', comp='1e+08')
    r1e09, _    = simulation.isodensities_lines(rad='full', comp='1e+09')
    r1e10, _    = simulation.isodensities_lines(rad='full', comp='1e+10')
    r1e11, gPNS = simulation.isodensities_lines(rad='full', comp='1e+11')
    r1e12, _    = simulation.isodensities_lines(rad='full', comp='1e+12')
    r1e13, _    = simulation.isodensities_lines(rad='full', comp='1e+13')
    r1e14, _    = simulation.isodensities_lines(rad='full', comp='1e+14')
    
    
    r1e08 = r1e08.data
    r1e09 = r1e09.data
    r1e10 = r1e10.data
    r1e11 = r1e11.data
    r1e12 = r1e12.data
    r1e13 = r1e13.data
    r1e14 = r1e14.data
    

    ## Get the elements
    dV = simulation.cell.dVolume_integration(simulation.ghost)
    radius = simulation.cell.radius(simulation.ghost)
    sinth = np.sin(simulation.cell.theta(simulation.ghost))
    if simulation.dim == 1:
        rt = radius ** 2 * dV
    elif simulation == 2:
        rt = (radius[None, :] * sinth[:, None]) ** 2 * dV
    else:
        rt = (radius[None, None, :] * sinth[None, :, None]) ** 2 * dV
    
    for file in simulation.hdf_file_list[start_point:]:
        progressBar(progress_index, total_points,
                    suffix='Computing inertia moment...')
        rho = simulation.rho(file) * rt
        tm  = simulation.time(file)
        ## compute the inertia moments at specific isolines
        i08 = compute_inertia(simulation, rho, radius, r1e08[..., findex],giso)
        i09 = compute_inertia(simulation, rho, radius, r1e09[..., findex],giso)
        i10 = compute_inertia(simulation, rho, radius, r1e10[..., findex],giso)
        i11 = compute_inertia(simulation, rho, radius, r1e11[..., findex],gPNS)
        i12 = compute_inertia(simulation, rho, radius, r1e12[..., findex],giso)
        i13 = compute_inertia(simulation, rho, radius, r1e13[..., findex],giso)
        i14 = compute_inertia(simulation, rho, radius, r1e14[..., findex],giso)

        try:
            time = np.concatenate((time, tm))
            iso8 = np.concatenate((iso8, i08))
            iso9 = np.concatenate((iso9, i09))
            iso10 = np.concatenate((iso10, i10))
            iso11 = np.concatenate((iso11, i11))
            iso12 = np.concatenate((iso12, i12))
            iso13 = np.concatenate((iso13, i13))
            iso14 = np.concatenate((iso14, i14))            
        except Exception as e:
            print(e)
            time = tm
            iso8 = i08
            iso9 = i09
            iso10 = i10
            iso11 = i11
            iso12 = i12
            iso13 = i13
            iso14 = i14
        processed_hdf.append(file)
        if (check_index >= checkpoint) and save_checkpoints:
            print('Checkpoint reached, saving...\n')
            save_hdf(os.path.join(simulation.storage_path,
                    'inertia_moment.h5'),
                     ['time', 'I8', 'I9', 'I10', 'I11', 'I12',
                      'I13', 'I14', 'processed'],
                     [time, iso8, iso9, iso10, iso11, iso12, iso13, iso14,
                      processed_hdf])
            
            check_index = 0
        check_index += 1
        progress_index += 1
        findex += 1   
    print('Computation completed, saving...')
    save_hdf(os.path.join(simulation.storage_path,
                    'inertia_moment.h5'),
                     ['time', 'I8', 'I9', 'I10', 'I11', 'I12',
                      'I13', 'I14', 'processed'],
                     [time, iso8, iso9, iso10, iso11, iso12, iso13, iso14,
                      processed_hdf])
    time, iso14, iso13, iso12, iso11, iso10, iso9, iso8, processed_hdf = \
            read_inertia_moment(simulation)
    return create_series(time, iso14, iso13, iso12, iso11, iso10, iso9, iso8)


def read_inertia_moment(simulation):
    inertia_file = h5py.File(os.path.join(simulation.storage_path, 
                                            'inertia_moment.h5'), 'r')
    data = [
        aerray(inertia_file['time'][...], u.s, 'time',
               r'$t-t_\mathrm{b}$', None, [-0.005, inertia_file['time'][-1]]),
        aerray(inertia_file['I14'][...], u.cm ** 2 * u.g, 'Inertia_moment_1e14',
               r'$I_{1e14}$', None, [10, 1e47], False),
        aerray(inertia_file['I13'][...], u.cm ** 2 * u.g, 'Inertia_moment_1e13',
               r'$I_{1e13}$', None, [10, 1e47], False),
        aerray(inertia_file['I12'][...], u.cm ** 2 * u.g, 'Inertia_moment_1e12',
               r'$I_{1e12}$', None, [10, 1e47], False),
        aerray(inertia_file['I11'][...], u.cm ** 2 * u.g, 'Inertia_moment_PNS',
               r'$I_\mathrm{PNS}$', None, [10, 1e47], False),
        aerray(inertia_file['I10'][...], u.cm ** 2 * u.g, 'Inertia_moment_1e10',
               r'$I_{1e10}$', None, [10, 1e47], False),
        aerray(inertia_file['I9'][...], u.cm ** 2 * u.g, 'Inertia_moment_1e09',
               r'$I_{1e09}$', None, [10, 1e47], False),
        aerray(inertia_file['I8'][...], u.cm ** 2 * u.g, 'Inertia_moment_1e08',
               r'$I_{1e08}$', None, [10, 1e47], False),
    ]
    data.append(inertia_file['processed'][...])
    inertia_file.close()
    return data
    