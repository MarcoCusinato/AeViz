from AeViz.quantities_plotting.plotting import Plotting
from typing import Literal
from AeViz.utils.aeviz_utils import plot_qt
import os

class AeViz(Plotting):
    def __init__(self):
        Plotting.__init__(self)
    
    def rho(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        qt = 'rho'
        if projection == '1D':
            if plane == 'time':
                qt = 'rho_max'
            self.plot1D(file, qt, plane, index1, index2)
        elif projection == '2D':
            if plane == 'time':
                self.plotProfile(qt)
            else:
                self.plot2D(file, plane, index1, qt)
    
    def MHD_energy(self, file=None, projection:Literal['1D', '2D']='1D',
                   index1=None, index2=None, plane:Literal['xy', 'xz',
                        'radius', 'theta', 'phi', 'time']='radius'):
        qt = 'MHD_energy'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def velocity(self, file=None, projection:Literal['1D', '2D']='1D', 
                 index1=None, index2=None, comp:Literal['all', 'r', 'th',
                 'ph']='all', plane:Literal['xy', 'xz', 'radius', 'theta',
                                            'phi', 'time']='radius'):
        if comp == 'all':
            qt = ['radial_velocity', 'theta_velocity', 'phi_velocity']
        elif comp == 'r':
            qt = ['radial_velocity']
        elif comp == 'th':
            qt = ['theta_velocity']
        elif comp == 'ph':
            qt = ['phi_velocity']
        plot_qt(self, file, qt, projection, index1, index2, plane)

    def soundspeed(self, file=None, projection:Literal['1D', '2D']='1D',
                   index1=None, index2=None, plane:Literal['xy', 'xz',
                   'radius', 'theta', 'phi', 'time']='radius'):
        qt = 'soundspeed'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def alfven_velocity(self, file=None, projection:Literal['1D', '2D']='1D',
                        index1=None, index2=None, plane:Literal['xy', 'xz',
                                'radius', 'theta', 'phi', 'time']='radius'):
        qt = 'alfven_velocity'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def kick_velocity(self,
                      comp:Literal['modulus', 'nu', 'hydro', 'all']='all'):
        plot_qt(self, None, 'kick_velocity_'+comp, '1D', None, None, 'time')

    def convective_velocity(self, file=None,
                         comp:Literal['all', 'convective', 'turbulent']='all'):
        if comp == 'all':
            qt = ['convective_velocity', 'turbulent_velocity']
        elif comp == 'convective':
            qt = ['convective_velocity']
        elif comp == 'turbulent':
            qt = ['turbulent_velocity']

        for q in qt:
            self.plot1D(file, q, 'radius', None, None)
    
    def omega(self, file=None, projection:Literal['1D', '2D']='1D', 
              index1=None, index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        qt = 'omega'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def convective_flux(self, file=None,
                        plane:Literal['time', 'radius']='radius'):
        qt = 'convective_flux'
        if plane == 'time':
            self.plotProfile(qt)
        else:
            self.plot1D(file, qt, plane, None, None)
    
    def Rossby_number(self, file=None, projection:Literal['1D', '2D']='1D',
                      index1=None, index2=None, plane:Literal['xy', 'xz',
                      'radius', 'theta', 'phi', 'time']='radius'):
        qt = 'Rossby_number'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def BV_frequency(self, file=None, projection:Literal['1D', '2D']='1D',
                        index1=None, index2=None, plane:Literal['xy', 'xz',
                        'radius', 'theta', 'phi', 'time']='radius'):
        qt = 'BV_frequency'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def gas_pressure(self, file=None, projection:Literal['1D', '2D']='1D',
                     index1=None, index2=None, plane:Literal['xy', 'xz',
                     'radius', 'theta', 'phi', 'time']='radius'):
        qt = 'gas_pressure'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def temperature(self, file=None, projection:Literal['1D', '2D']='1D',
                    index1=None, index2=None, plane:Literal['xy', 'xz',
                    'radius', 'theta', 'phi', 'time']='radius'):
        qt = 'temperature'
        plot_qt(self, file, qt, projection, index1, index2, plane)

    def enthalphy(self, file=None, projection:Literal['1D', '2D']='1D',
                  index1=None, index2=None, plane:Literal['xy', 'xz',
                  'radius', 'theta', 'phi', 'time']='radius'):
        qt = 'enthalphy'
        plot_qt(self, file, qt, projection, index1, index2, plane)

    def entropy(self, file=None, projection:Literal['1D', '2D']='1D',
                index1=None, index2=None, plane:Literal['xy', 'xz', 'radius',
                'theta', 'phi', 'time']='radius'):
        qt = 'entropy'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def adiabatic_index(self, file=None, projection:Literal['1D', '2D']='1D',
                        index1=None, index2=None, plane:Literal['xy', 'xz',
                        'radius', 'theta', 'phi', 'time']='radius'):
        qt = 'adiabatic_index'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def lorentz_factor(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        qt = 'lorentz'
        plot_qt(self, file, qt, projection, index1, index2, plane)

    def grav_pot(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        qt = 'gravitational_potential'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def grav_ene(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        qt = 'gravitational_energy'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def int_ene(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        qt = 'internal_energy'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def nu_heat(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        qt = 'nu_heat'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def Ye(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        qt = 'Ye'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def Xn(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        qt = 'neutron_fraction'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def Xp(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        qt = 'proton_fraction'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def Xalpha(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        qt = 'alpha_fraction'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def Xheavy(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        qt = 'heavy_fraction'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def Abar(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        qt = 'Abar'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def Zbar(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        qt = 'Zbar'
        plot_qt(self, file, qt, projection, index1, index2, plane)

    def cpot_e(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        qt = 'electron_chemical_potential'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def cpot_n(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        qt = 'neutron_chemical_potential'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def cpot_p(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        qt = 'proton_chemical_potential'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def cpot_nu(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        qt = 'neutrino_chemical_potential'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def error(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        qt = 'error'
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def magnetic_fields(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, comp:Literal['all', 'r', 'th', 'ph', 'torpol', 
                                      'poloidal', 'toroidal']='all',
            plane:Literal['xy', 'xz', 'radius', 'theta', 'phi', 'time']='radius'):
        if comp == 'all':
            qt = ['BX', 'BY', 'BZ']
        elif comp == 'r':
            qt = ['BX']
        elif comp == 'th':
            qt = ['BY']
        elif comp == 'ph':
            qt = ['BZ']
        elif comp == 'poloidal':
            qt = ['poloidal_magnetic_fields']
        elif comp == 'toroidal':
            qt = ['toroidal_magnetic_fields']
        elif comp == 'torpol':
            qt = ['poloidal_magnetic_fields', 'toroidal_magnetic_fields']
        plot_qt(self, file, qt, projection, index1, index2, plane)

    def magnetic_energy(self, file=None, projection:Literal['1D', '2D']='1D', 
                        index1=None, index2=None, comp:Literal['all', 'tot', 
                                                'poloidal', 'toroidal']='all', 
                        plane:Literal['xy', 'xz', 'radius', 'theta', 'phi',
                                     'time']='radius'):
        if comp == 'all':
            qt = ['total_magnetic_energy', 'poloidal_magnetic_energy',
                  'toroidal_magnetic_energy']
        elif comp == 'tot':
            qt = ['total_magnetic_energy']
        elif comp == 'poloidal':
            qt = ['poloidal_magnetic_energy']
        elif comp == 'toroidal':
            qt = ['toroidal_magnetic_energy']
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def nue_moment(self, file=None, projection:Literal['1D', '2D']='1D',
                  comp:Literal['x', 'y', 'z', 'e', 'all']='all', index1=None,
                  index2=None, plane:Literal['xy', 'xz', 'radius',
                  'theta', 'phi', 'time']='radius'):
        if comp == 'all':
            qt = ['nue_moment_x', 'nue_moment_y', 'nue_moment_z',
                  'nue_moment_e']
        else:
            qt = ['nue_moment_' + comp]
        plot_qt(self, file, qt, projection, index1, index2, plane)
        
    def nua_moment(self, file=None, projection:Literal['1D', '2D']='1D',
                  comp:Literal['x', 'y', 'z', 'e', 'all']='all', index1=None,
                  index2=None, plane:Literal['xy', 'xz', 'radius',
                  'theta', 'phi', 'time']='radius'):
        if comp == 'all':
            qt = ['nua_moment_x', 'nua_moment_y', 'nua_moment_z',
                  'nua_moment_e']
        else:
            qt = ['nua_moment_' + comp]
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def nux_moment(self, file=None, projection:Literal['1D', '2D']='1D',
                  comp:Literal['x', 'y', 'z', 'e', 'all']='all', index1=None,
                  index2=None, plane:Literal['xy', 'xz', 'radius',
                  'theta', 'phi', 'time']='radius'):
        if comp == 'all':
            qt = ['nux_moment_x', 'nux_moment_y', 'nux_moment_z',
                  'nux_moment_e']
        else:
            qt = ['nux_moment_' + comp]
        plot_qt(self, file, qt, projection, index1, index2, plane)
        
    def nu_mean_ene(self, file=None, projection:Literal['1D', '2D']='1D',
                    comp:Literal['nue', 'nua', 'nux', 'all']='all', index1=None,
                    index2=None, plane:Literal['xy', 'xz', 'radius', 'theta',
                                               'phi', 'time']='radius'):
        if comp == 'all':
            qt = ['nue_mean_ene', 'nua_mean_ene', 'nux_mean_ene']
        else:
            qt = [comp + '_mean_ene']
        plot_qt(self, file, qt, projection, index1, index2, plane)
    
    def nu_integrated(self, flavour:Literal['all', 'nue', 'nua', 'nux']='all',
                      comp: Literal['ene', 'lum']='lum'):
        qt = 'nu_integrated_' + comp + '_' + flavour
        self.plot1D(None, qt, 'time', None, None)
        
    def PNS_radius(self, comp:Literal['all', 'min', 'max', 'avg']='avg'):
        self.plot1D(None, 'PNS_radius_' + comp, 'time', None, None)
    
    def shock_radius(self, comp:Literal['all', 'min', 'max','avg']='avg'):
        self.plot1D(None, 'shock_radius_' + comp, 'time', None, None)
    
    def neutrino_spheres(self, comp:Literal['all', 'min', 'max', 'avg']='avg'):
        for flavour in ['nue', 'nua', 'nux']:
            self.plot1D(None, 'neutrino_spheres_' + flavour + '_' + comp,
                        'time', None, None)
    
    def gain_radius(self,  comp:Literal['all', 'min', 'max', 'avg']='avg'):
        self.plot1D(None, 'gain_radius_' + comp, 'time', None, None)
    
    def innercore_radius(self, comp:Literal['all', 'min', 'max','avg']='avg'):
        self.plot1D(None, 'innercore_radius_' + comp, 'time', None, None)
    
    def PNS_nucleus_radius(self, comp:Literal['all', 'min', 'max',
                                              'avg']='avg'):
        self.plot1D(None, 'PNS_nucleus_radius_' + comp, 'time', None, None)
    
    def explosion(self, comp: Literal['all', 'mass', 'ene', 'kin', 'mag']):
        if comp == 'all':
            qt = ['mass', 'ene', 'kin', 'mag']
        else:
            qt = [comp]
        for q in qt:
            self.plot1D(None, 'explosion_' + q, 'time', None, None)
    
    def gain(self, comp: Literal['all', 'mass', 'ene']):
        if comp == 'all':
            qt = ['mass', 'ene']
        else:
            qt = [comp]
        for q in qt:
            self.plot1D(None, 'gain_' + q, 'time', None, None)
    
    def innercore(self, comp: Literal['mass', 'ene', 'kin', 'mag', 'rot',
                                      'grav', 'T/W']):
        self.plot1D(None, 'innercore_' + comp, 'time', None, None)
    
    def PNS(self, comp: Literal['mass', 'ene', 'kin', 'mag', 'rot', 'grav',
                               'conv']):
        self.plot1D(None, 'PNS_' + comp, 'time', None, None)
    
    def PNS_angular_mom(self, comp: Literal['all', 'Lx', 'Ly', 'Lz', 'Ltot']):
        
        self.plot1D(None, 'PNS_angular_mom_' + comp, 'time', None, None)
    
    def mass_accretion(self):
        self.plot1D(None, 'mass_accretion_500km', 'time', None, None)
    
    def GWs(self, comp:Literal['all', 'h+eq', 'h+pol', 'hxeq', 'hxpol']='h+eq',
            decomposition=False, spectrogram=False):
        if self.sim_dim == 1:
            pass
        else:
            if comp == 'all' and self.sim_dim == 3:
                qt = ['h+eq', 'h+pol', 'hxeq', 'hxpol']
            elif self.sim_dim == 2:
                qt = ['h+eq']
            else:
                qt = [comp]
            if decomposition:
                for q in qt:
                    self.plotGWDecomposition(q)
                    
            elif spectrogram:
                for q in qt:
                    self.plotGWspectrogram(q)
            else:
                for q in qt:
                    self.plot1D(None, 'GW_Amplitudes_' + q, 'time', None, None)
    
    def add_field(self, file, plot, comp: Literal['velocity', 'Bfield'],
                  plane: Literal['xy', 'yz']='xz', index1=None):
        """
        Adds velocity or magnetic field to the selected plot.
        """
        self.add_2Dfield(file, plot, comp, plane, index1)

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

    def save_plot(self, name):
        if name.endswith('.png') or name.endswith('.pdf') or \
            name.endswith('.jpg'):
            self.fig.savefig(os.path.join(self.save_path, name))
        else:
            self.fig.savefig(os.path.join(self.save_path, name + '.png'))