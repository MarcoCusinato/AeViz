from AeViz.quantities_plotting.plotting import Plotting
from typing import Literal
from AeViz.utils.files.file_utils import list_module_functions
from AeViz.utils.decorators.figure import fig_window_open
import os
import types

class AeViz(Plotting):
    def __init__(self):
        Plotting.__init__(self)
        self.__loaded_methods = []
        self.__setup_method_lists()
    
    
    def Load(self, path, simulation_path=None, dim=None):
        ## Clear loaded methods
        self.__delete_methods()
        ## Load the simulation
        super().Load(path, simulation_path, dim)
        ## Load the methods
        if self.data_type == 'hdf':
            self.__load_hdf_methods()
        elif self.data_type == 'sim':
            self.__load_simulation_methods()
        else:
            raise ValueError('Data type not recognized.')

    @fig_window_open
    def movie(self, qt1=None, qt2=None, qt3=None, qt4=None, top_qt=None,
              fields: Literal['velocity', 'Bfield', 'all']=None,
              plane: Literal['xy', 'yz']='xz', 
              start_time=None, end_time=None, lims=None):
        if fields == 'all':
            vf, bf = True, True
        elif fields == 'velocity':
            vf, bf = True, False
        elif fields == 'Bfield':
            vf, bf = False, True
        else:
            vf, bf = False, False

        self.make_movie(qt1=qt1, qt2=qt2, qt3=qt3, qt4=qt4, top=top_qt,
              plane=plane, start_time=start_time, lims=lims,
              end_time=end_time, vfield=vf, Bfield=bf, top_time=True)

    @fig_window_open
    def save_plot(self, name, **kwargs):
        savedir = kwargs.pop('savedir', self.save_path)
        if name.endswith('.png') or name.endswith('.pdf') or \
            name.endswith('.jpg'):
            self.fig.savefig(os.path.join(savedir, name))
        else:
            self.fig.savefig(os.path.join(savedir, name + '.png'))
    
    def get_loaded_methods(self):
        return self.__loaded_methods
    
    def __setup_method_lists(self):
        import AeViz.AeVizMethods.hydro
        import AeViz.AeVizMethods.composition
        import AeViz.AeVizMethods.thermo
        import AeViz.AeVizMethods.aux
        import AeViz.AeVizMethods.neutrinos
        import AeViz.AeVizMethods.magnetic_fields
        import AeViz.AeVizMethods.instabilities
        import AeViz.AeVizMethods.GWs
        import AeViz.AeVizMethods.other_spherical_sym
        import AeViz.AeVizMethods.supernova
        self.__hdf_methods = list_module_functions(AeViz.AeVizMethods.hydro)
        self.__hdf_methods.extend(list_module_functions(
            AeViz.AeVizMethods.composition))
        self.__hdf_methods.extend(list_module_functions(
            AeViz.AeVizMethods.thermo))
        self.__hdf_methods.extend(list_module_functions(
            AeViz.AeVizMethods.aux))
        self.__hdf_methods.extend(list_module_functions(
            AeViz.AeVizMethods.neutrinos))
        self.__hdf_methods.extend(list_module_functions(
            AeViz.AeVizMethods.magnetic_fields))
        self.__simulation_methods = self.__hdf_methods.copy()
        self.__supernova_methods = self.__simulation_methods.copy()
        self.__supernova_methods.extend(list_module_functions(
            AeViz.AeVizMethods.supernova))
        self.__supernova_methods.extend(list_module_functions(
            AeViz.AeVizMethods.other_spherical_sym))
        self.__supernova_methods.extend(list_module_functions(
            AeViz.AeVizMethods.instabilities))
        self.__supernova_methods.extend(list_module_functions(
            AeViz.AeVizMethods.GWs))

    def __delete_methods(self):
        for name in self.__loaded_methods:
            if hasattr(self, name):
                delattr(self, name)
    
    def __load_hdf_methods(self):
        for name, obj in self.__hdf_methods:
            if name in ['omega', 'add_field'] :
                continue
            setattr(self, name, types.MethodType(obj, self))
            self.__loaded_methods.append(name)
    
    def __load_simulation_methods(self):
        ignore_list = []
        if self.loaded_data.GEOM != 2 and self.loaded_data.dim < 2:
            ignore_list.extend(['omega', 'GWs', 'IMFs', 'epicyclic_frequency'])
        if self.loaded_data.evolved_qts['magdim'] < 1:
            ignore_list.extend(['alfven_velocity', 'magnetic_fields',
                                'magnetic_energy'])
        if self.loaded_data.evolved_qts['neudim'] < 1:
            ignore_list.extend(['nue_moment', 'nua_moment', 'nux_moment',
                                'nu_mean_ene', 'nu_integrated'])
        if self.loaded_data.evolved_qts['tempratr'] < 1:
            ignore_list.extend(['temperature'])
        if self.loaded_data.evolved_qts['entropie'] < 1:
            ignore_list.extend(['entropy'])
        if self.loaded_data.evolved_qts['y_edim'] < 1:
            ignore_list.extend(['Ye'])
        if not self.loaded_data.relativistic:
            ignore_list.extend(['lorentz_factor'])
        if self.loaded_data.evolved_qts['comp_dim'] < 1:
            ignore_list.extend(['Xp', 'Xn', 'Xalpha', 'Xheavy', 'Abar',
                                'Zbar'])
        if self.loaded_data.evolved_qts['cpot_dim'] < 1:
            ignore_list.extend(['cpot_e', 'cpot_n', 'cpot_p', 'cpot_nu'])
        if self.loaded_data.GEOM != 2:
            for name, obj in self.__simulation_methods:
                if name in ignore_list:
                    continue
                setattr(self, name, types.MethodType(obj, self))
                self.__loaded_methods.append(name)
        else:
            try:
                import PyEMD
            except:
                ignore_list.append('IMFs')
            for name, obj in self.__supernova_methods:
                if name in ignore_list:
                    continue
                setattr(self, name, types.MethodType(obj, self))
                self.__loaded_methods.append(name)

        

