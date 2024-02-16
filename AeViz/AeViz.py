from AeViz.quantities_plotting.plotting import Plotting
from typing import Literal

class AeViz(Plotting):
    def __init__(self):
        Plotting.__init__(self)
    
    def rho(self, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def MHD_energy(self, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def velocity(self, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, comp:Literal['all', 'r', 'th', 'ph', 'sound']='all',
            plane:Literal['xy', 'xz', 'radius', 'theta', 'phi', 'time']='radius'):
        pass
    
    def omega(self, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def gas_pressure(self, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def temperature(self, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass

    def enthalphy(self, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass

    def entropy(self, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def adiabatic_index(self, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def lorentz_factor(self, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass

    def grav_pot(self, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def grav_ene(self, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def int_ene(self, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def nu_heat(self, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def Ye(self, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def Xn(self, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def Xp(self, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def Xalpha(self, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass
    
    def Xheavy(self, projection:Literal['1D', '2D']='1D', index1=None, 
            index2=None, plane:Literal['xy', 'xz', 'radius',
                          'theta', 'phi', 'time']='radius'):
        pass