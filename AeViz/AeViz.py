from AeViz.quantities_plotting.plotting import Plotting
import AeViz.utils.radii_utils
from typing import Literal
from AeViz.utils.GW_utils import GW_spectrogram

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
                self.plot2D(file, qt)
    
    def MHD_energy(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        if projection == '1D':
            self.plot1D(file, 'MHD_energy', plane, index1, index2)
        elif projection == '2D':
            self.plot2D(file,'MHD_energy')

    
    def velocity(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, comp:Literal['all', 'r', 'th', 'ph']='all',
            plane:Literal['xy', 'xz', 'radius', 'theta', 'phi', 
                          'time']='radius'):
        pass
    def soundspeed(self, file=None, projection:Literal['1D', '2D']='1D', index1=None,
            index2=None, plane:Literal['xy', 'xz', 'radius', 'theta', 'phi',
                                              'time']='radius'):
            pass
    
    def alfven_velocity(self, file=None, projection:Literal['1D', '2D']='1D', index1=None,
            index2=None, plane:Literal['xy', 'xz', 'radius', 'theta', 'phi',
                          'time']='radius'):
        pass
    
    def convection_velocity(self, file=None, projection:Literal['1D', '2D']='1D', index1=None,
            index2=None, comp:Literal['all', 'convection', 'turbulent']='all',
            plane:Literal['xy', 'xz', 'radius', 'theta', 'phi', 
                          'time']='radius'):
        pass
    
    def omega(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def gas_pressure(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def temperature(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass

    def enthalphy(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass

    def entropy(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def adiabatic_index(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def lorentz_factor(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass

    def grav_pot(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def grav_ene(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def int_ene(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def nu_heat(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def Ye(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def Xn(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def Xp(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def Xalpha(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def Xheavy(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def Abar(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def Zbar(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass

    def cpot_e(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def cpot_n(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def cpot_p(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def cpot_nu(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def error(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def magnetic_fields(self, file=None, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, comp:Literal['all', 'r', 'th', 'ph', 'poloidal',
                                      'toroidal']='all',
            plane:Literal['xy', 'xz', 'radius', 'theta', 'phi', 'time']='radius'):
        pass
    
    def magnetic_energy(self, file=None, projection:Literal['1D', '2D']='1D', 
                        index1=None, index2=None, comp:Literal['tot', 
                                                'poloidal', 'toroidal']='tot', 
                        plane:Literal['xy', 'xz', 'radius', 'theta', 'phi',
                                     'time']='radius'):
        pass
    
    def PNS_radius(self, comp:Literal['min', 'max', 'avg']='avg'):
        pass
    
    def shock_radius(self, comp:Literal['min', 'max','avg']='avg'):
        pass
    
    def neutrino_spheres(self, comp:Literal['min', 'max', 'avg']='avg'):
        pass
    
    def gain_radius(self,  comp:Literal['min', 'max', 'avg']='avg'):
        pass
    
    def innercore_radius(self, comp:Literal['min', 'max','avg']='avg'):
        pass
    
    def PNS_nucleus_radius(self, comp:Literal['min', 'max', 'avg']='avg'):
        pass
    
    def explosion(self, comp: Literal['mass', 'ene', 'kin', 'mag']):
        pass
    
    def gain(self, comp: Literal['mass', 'ene']):
        pass
    
    def innercore(self, comp: Literal['mass', 'ene', 'kin', 'mag', 'rot',
                                      'grav', 'tot', 'T_W']):
        pass
    
    def mass_accretion(self):
        pass
    
    def GWs(self, comp:Literal['all', 'h+eq', 'h+pol', 'hxeq', 'hxpol']='h+eq',
            decomposition=False, GW_spectrogram=False):
        pass