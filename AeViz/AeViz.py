from AeViz.quantities_plotting.plotting import Plotting
from typing import Literal
from AeViz.utils.aeviz_utils import plot_qt
from AeViz.utils.decorators import fig_window_open
import os

class AeViz(Plotting):
    def __init__(self):
        Plotting.__init__(self)
    
    @fig_window_open
    def rho(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                          'theta', 'phi', 'time']='radius', **kwargs):
        qt = 'rho'
        if projection == '1D':
            if plane == 'time':
                qt = 'rho_max'
            self.plot1D(file, qt, plane, index1, index2, **kwargs)
        elif projection == '2D':
            if plane == 'time':
                self.plotProfile(qt)
            else:
                self.plot2D(file, plane, qt, **kwargs)
    
    @fig_window_open
    def MHD_energy(self, file=None, projection:Literal['1D', '2D']='1D',
                   index1=None, index2=None, plane:Literal['xz', 'yz', 'xy',
                        'radius', 'theta', 'phi', 'time']='radius', **kwargs):
        qt = 'MHD_energy'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def velocity(self, file=None, projection:Literal['1D', '2D']='1D', 
                 index1=None, index2=None, comp:Literal['all', 'r', 'th',
                 'ph']='all', plane:Literal['xz', 'yz', 'xy', 'radius', 'theta',
                                            'phi', 'time']='radius', **kwargs):
        if comp == 'all':
            qt = ['radial_velocity', 'theta_velocity', 'phi_velocity']
        elif comp == 'r':
            qt = ['radial_velocity']
        elif comp == 'th':
            qt = ['theta_velocity']
        elif comp == 'ph':
            qt = ['phi_velocity']
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

    @fig_window_open
    def soundspeed(self, file=None, projection:Literal['1D', '2D']='1D',
                   index1=None, index2=None, plane:Literal['xz', 'yz', 'xy',
                   'radius', 'theta', 'phi', 'time']='radius', **kwargs):
        qt = 'soundspeed'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def alfven_velocity(self, file=None, projection:Literal['1D', '2D']='1D',
                        index1=None, index2=None, plane:Literal['xz', 'yz', 'xy',
                                'radius', 'theta', 'phi', 'time']='radius',
                        **kwargs):
        qt = 'alfven_velocity'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def kick_velocity(self,
                      comp:Literal['modulus', 'nu', 'hydro', 'all']='all',
                      **kwargs):
        plot_qt(self, None, 'kick_velocity_'+comp, '1D', None, None, 'time')

    @fig_window_open
    def convective_velocity(self, file=None,
                         comp:Literal['all', 'convective', 'turbulent']='all',
                         **kwargs):
        if comp == 'all':
            qt = ['convective_velocity', 'turbulent_velocity']
        elif comp == 'convective':
            qt = ['convective_velocity']
        elif comp == 'turbulent':
            qt = ['turbulent_velocity']

        for q in qt:
            self.plot1D(file, q, 'radius', None, None, **kwargs)
    
    @fig_window_open
    def omega(self, file=None, projection:Literal['1D', '2D']='1D', 
              index1=None, index2=None, plane:Literal['xz', 'yz', 'xy',
                                                      'radius',
                          'theta', 'phi', 'time']='radius', **kwargs):
        qt = 'omega'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def convective_flux(self, file=None,
                        plane:Literal['time', 'radius']='radius', **kwargs):
        qt = 'convective_flux'
        if plane == 'time':
            self.plotProfile(qt)
        else:
            self.plot1D(file, qt, plane, None, None, **kwargs)
    
    @fig_window_open
    def Rossby_number(self, file=None, projection:Literal['1D', '2D']='1D',
                      index1=None, index2=None, plane:Literal['xz', 'yz', 'xy',
                      'radius', 'theta', 'phi', 'time']='radius', **kwargs):
        qt = 'Rossby_number'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def BV_frequency(self, file=None, projection:Literal['1D', '2D']='1D',
                        index1=None, index2=None,
                        plane:Literal['xz', 'yz', 'xy', 'radius', 'theta',
                                      'phi', 'time']='radius', **kwargs):
        qt = 'BV_frequency'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def gas_pressure(self, file=None, projection:Literal['1D', '2D']='1D',
                     index1=None, index2=None, plane:Literal['xz', 'yz', 'xy',
                     'radius', 'theta', 'phi', 'time']='radius', **kwargs):
        qt = 'gas_pressure'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def temperature(self, file=None, projection:Literal['1D', '2D']='1D',
                    index1=None, index2=None, plane:Literal['xz', 'yz', 'xy',
                    'radius', 'theta', 'phi', 'time']='radius', **kwargs):
        qt = 'temperature'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

    @fig_window_open
    def enthalphy(self, file=None, projection:Literal['1D', '2D']='1D',
                  index1=None, index2=None, plane:Literal['xz', 'yz', 'xy',
                  'radius', 'theta', 'phi', 'time']='radius', **kwargs):
        qt = 'enthalphy'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

    @fig_window_open
    def entropy(self, file=None, projection:Literal['1D', '2D']='1D',
                index1=None, index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                'theta', 'phi', 'time']='radius', **kwargs):
        qt = 'entropy'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def adiabatic_index(self, file=None, projection:Literal['1D', '2D']='1D',
                        index1=None, index2=None, plane:Literal['xz', 'yz', 'xy',
                        'radius', 'theta', 'phi', 'time']='radius', **kwargs):
        qt = 'adiabatic_index'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def lorentz_factor(self, file=None, projection:Literal['1D', '2D']='1D',
                       index1=None, index2=None,
                       plane:Literal['xz', 'yz', 'xy', 'radius', 'theta', 
                                     'phi', 'time']='radius', **kwargs):
        qt = 'lorentz'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

    @fig_window_open
    def grav_pot(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                          'theta', 'phi', 'time']='radius', **kwargs):
        qt = 'gravitational_potential'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def grav_ene(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                          'theta', 'phi', 'time']='radius', **kwargs):
        qt = 'gravitational_energy'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def int_ene(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                          'theta', 'phi', 'time']='radius', **kwargs):
        qt = 'internal_energy'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def nu_heat(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                          'theta', 'phi', 'time']='radius', **kwargs):
        qt = 'nu_heat'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def Ye(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                          'theta', 'phi', 'time']='radius', **kwargs):
        qt = 'Ye'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def Xn(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                          'theta', 'phi', 'time']='radius', **kwargs):
        qt = 'neutron_fraction'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def Xp(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                          'theta', 'phi', 'time']='radius', **kwargs):
        qt = 'proton_fraction'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def Xalpha(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                          'theta', 'phi', 'time']='radius', **kwargs):
        qt = 'alpha_fraction'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def Xheavy(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                          'theta', 'phi', 'time']='radius', **kwargs):
        qt = 'heavy_fraction'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def Abar(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                          'theta', 'phi', 'time']='radius', **kwargs):
        qt = 'Abar'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def Zbar(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                          'theta', 'phi', 'time']='radius', **kwargs):
        qt = 'Zbar'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

    @fig_window_open
    def cpot_e(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                          'theta', 'phi', 'time']='radius', **kwargs):
        qt = 'electron_chemical_potential'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def cpot_n(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                          'theta', 'phi', 'time']='radius', **kwargs):
        qt = 'neutron_chemical_potential'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def cpot_p(self, file=None, projection:Literal['1D', '2D']='1D',
               index1=None, index2=None,
               plane:Literal['xz', 'yz', 'xy', 'radius', 'theta', 'phi',
                             'time']='radius', **kwargs):
        qt = 'proton_chemical_potential'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def cpot_nu(self, file=None, projection:Literal['1D', '2D']='1D',
                index1=None, index2=None,
                plane:Literal['xz', 'yz', 'xy', 'radius', 'theta', 'phi',
                              'time']='radius', **kwargs):
        qt = 'neutrino_chemical_potential'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def error(self, file=None, projection:Literal['1D', '2D']='1D',
              index1=None, index2=None,
              plane:Literal['xz', 'yz', 'xy', 'radius', 'theta', 'phi',
                            'time']='radius', **kwargs):
        qt = 'error'
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def magnetic_fields(self, file=None, projection:Literal['1D', '2D']='1D',
                        index1=None, index2=None,
                        comp:Literal['all', 'r', 'th', 'ph', 'torpol',
                                     'poloidal', 'toroidal']='all',
                        plane:Literal['xz', 'yz', 'xy', 'radius', 'theta',
                                      'phi', 'time']='radius', **kwargs):
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
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)

    @fig_window_open
    def magnetic_energy(self, file=None, projection:Literal['1D', '2D']='1D', 
                        index1=None, index2=None,
                        comp:Literal['all', 'tot', 'poloidal',
                                     'toroidal']='all', 
                        plane:Literal['xz', 'yz', 'xy', 'radius', 'theta',
                                      'phi', 'time']='radius', **kwargs):
        if comp == 'all':
            qt = ['total_magnetic_energy', 'poloidal_magnetic_energy',
                  'toroidal_magnetic_energy']
        elif comp == 'tot':
            qt = ['total_magnetic_energy']
        elif comp == 'poloidal':
            qt = ['poloidal_magnetic_energy']
        elif comp == 'toroidal':
            qt = ['toroidal_magnetic_energy']
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def nue_moment(self, file=None, projection:Literal['1D', '2D']='1D',
                  comp:Literal['x', 'y', 'z', 'e', 'all']='all', index1=None,
                  index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                  'theta', 'phi', 'time']='radius', **kwargs):
        if comp == 'all':
            qt = ['nue_moment_x', 'nue_moment_y', 'nue_moment_z',
                  'nue_moment_e']
        else:
            qt = ['nue_moment_' + comp]
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
        
    @fig_window_open
    def nua_moment(self, file=None, projection:Literal['1D', '2D']='1D',
                  comp:Literal['x', 'y', 'z', 'e', 'all']='all', index1=None,
                  index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                  'theta', 'phi', 'time']='radius', **kwargs):
        if comp == 'all':
            qt = ['nua_moment_x', 'nua_moment_y', 'nua_moment_z',
                  'nua_moment_e']
        else:
            qt = ['nua_moment_' + comp]
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def nux_moment(self, file=None, projection:Literal['1D', '2D']='1D',
                  comp:Literal['x', 'y', 'z', 'e', 'all']='all', index1=None,
                  index2=None, plane:Literal['xz', 'yz', 'xy', 'radius',
                  'theta', 'phi', 'time']='radius', **kwargs):
        if comp == 'all':
            qt = ['nux_moment_x', 'nux_moment_y', 'nux_moment_z',
                  'nux_moment_e']
        else:
            qt = ['nux_moment_' + comp]
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
        
    @fig_window_open
    def nu_mean_ene(self, file=None, projection:Literal['1D', '2D']='1D',
                    comp:Literal['nue', 'nua', 'nux', 'all']='all',
                    index1=None, index2=None,
                    plane:Literal['xz', 'yz', 'xy', 'radius', 'theta', 'phi',
                                  'time']='radius', **kwargs):
        if comp == 'all':
            qt = ['nue_mean_ene', 'nua_mean_ene', 'nux_mean_ene']
        else:
            qt = [comp + '_mean_ene']
        plot_qt(self, file, qt, projection, index1, index2, plane, **kwargs)
    
    @fig_window_open
    def nu_integrated(self, flavour:Literal['all', 'nue', 'nua', 'nux']='all',
                      comp: Literal['ene', 'lum']='lum', **kwargs):
        qt = 'nu_integrated_' + comp + '_' + flavour
        self.plot1D(None, qt, 'time', None, None, **kwargs)
        
    @fig_window_open
    def PNS_radius(self, comp:Literal['all', 'min', 'max', 'avg']='avg',
                   **kwargs):
        self.plot1D(None, 'PNS_radius_' + comp, 'time', None, None, **kwargs)
    
    @fig_window_open
    def shock_radius(self, comp:Literal['all', 'min', 'max','avg']='avg',
                     **kwargs):
        self.plot1D(None, 'shock_radius_' + comp, 'time', None, None, **kwargs)
    
    @fig_window_open
    def neutrino_spheres(self, comp:Literal['all', 'min', 'max', 'avg']='avg',
                         **kwargs):
        for flavour in ['nue', 'nua', 'nux']:
            self.plot1D(None, 'neutrino_spheres_' + flavour + '_' + comp,
                        'time', None, None, **kwargs)
    
    @fig_window_open
    def gain_radius(self,  comp:Literal['all', 'min', 'max', 'avg']='avg',
                    **kwargs):
        self.plot1D(None, 'gain_radius_' + comp, 'time', None, None, **kwargs)
    
    @fig_window_open
    def innercore_radius(self, comp:Literal['all', 'min', 'max','avg']='avg',
                         **kwargs):
        self.plot1D(None, 'innercore_radius_' + comp, 'time', None, None, **kwargs)
    
    @fig_window_open
    def PNS_nucleus_radius(self, comp:Literal['all', 'min', 'max',
                                              'avg']='avg', **kwargs):
        self.plot1D(None, 'PNS_nucleus_radius_' + comp, 'time', None, None, **kwargs)
    
    @fig_window_open
    def explosion(self, comp: Literal['all', 'mass', 'ene', 'kin', 'mag'],
                  **kwargs):
        if comp == 'all':
            qt = ['mass', 'ene', 'kin', 'mag']
        else:
            qt = [comp]
        for q in qt:
            self.plot1D(None, 'explosion_' + q, 'time', None, None, **kwargs)
    
    @fig_window_open
    def gain(self, comp: Literal['all', 'mass', 'ene'], **kwargs):
        if comp == 'all':
            qt = ['mass', 'ene']
        else:
            qt = [comp]
        for q in qt:
            self.plot1D(None, 'gain_' + q, 'time', None, None, **kwargs)
    
    @fig_window_open
    def innercore(self, comp: Literal['mass', 'ene', 'kin', 'mag', 'rot',
                                      'grav', 'T/W'], **kwargs):
        self.plot1D(None, 'innercore_' + comp, 'time', None, None, **kwargs)
    
    @fig_window_open
    def PNS(self, comp: Literal['mass', 'ene', 'kin', 'mag', 'rot', 'grav',
                               'conv'], **kwargs):
        self.plot1D(None, 'PNS_' + comp, 'time', None, None, **kwargs)
    
    @fig_window_open
    def PNS_angular_mom(self, comp: Literal['all', 'Lx', 'Ly', 'Lz', 'Ltot'],
                        **kwargs):
        
        self.plot1D(None, 'PNS_angular_mom_' + comp, 'time', None, None,
                    **kwargs)
    
    @fig_window_open
    def mass_accretion(self, **kwargs):
        self.plot1D(None, 'mass_accretion_500km', 'time', None, None, **kwargs)
    
    @fig_window_open
    def GWs(self, comp:Literal['all', 'h+eq', 'h+pol', 'hxeq', 'hxpol']='h+eq',
            decomposition=False, spectrogram=False, **kwargs):
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
                    self.plotGWDecomposition(q, **kwargs)
                    
            elif spectrogram:
                for q in qt:
                    self.plotGWspectrogram(q)
            else:
                for q in qt:
                    self.plot1D(None, 'GW_Amplitudes_' + q, 'time', None, None,
                                **kwargs)
    
    @fig_window_open
    def IMFs(self, strain:Literal['h+eq', 'hxeq', 'h+pol', 'hxpol']='h+eq',
             mode:Literal['EMD', 'EEMD']='EMD', min_imfs=0, max_imfs=10,
             spectro=False, **kwargs):
        loc = locals()
        for q in loc.keys():
            if q not in ['self', 'kwargs']:
                kwargs[q] = loc[q]
        if spectro:
            self.plotHHT(**kwargs)
        else:
            raise NotImplementedError('IMFs plotting is not implemented yet.')
            self.plotIMFs(**kwargs)
    
    @fig_window_open
    def add_field(self, file, plot, comp: Literal['velocity', 'Bfield'],
                  plane: Literal['xy', 'yz', 'xz']='xz', **kwargs):
        """
        Adds velocity or magnetic field to the selected plot.
        """
        self.add_2Dfield(file, plot, comp, plane)

    @fig_window_open
    def rho_decomposition(self, l=0, m=None, correct_zero=True, **kwargs):
        loc = locals()
        for q in loc.keys():
            if q not in ['self', 'kwargs']:
                kwargs[q] = loc[q]
        qt = 'rho_spherical_harmonics'
        self.plotProfile(qt, **kwargs)
    
    @fig_window_open
    def rho_decomposition_barcode(self, lmin=None, lmax=None, msum=False,
                                  r=None, rhomin=None, rhomax=None,
                                  zero_norm=True, **kwargs):
        loc = locals()
        for q in loc.keys():
            if q not in ['self', 'kwargs']:
                kwargs[q] = loc[q]
        self.barcode(**kwargs)
        
        
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
    def save_plot(self, name):
        if name.endswith('.png') or name.endswith('.pdf') or \
            name.endswith('.jpg'):
            self.fig.savefig(os.path.join(self.save_path, name))
        else:
            self.fig.savefig(os.path.join(self.save_path, name + '.png'))