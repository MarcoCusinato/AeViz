from AeViz.utils.files.path_utils import (pltf, local_storage_folder,
                              save_paths_dictionary, get_paths_dictionary)


local_storage_folder(pltf())
sim_dict = get_paths_dictionary(True)
save_paths_dictionary(sim_dict)

from .simulation import Simulation