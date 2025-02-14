from astropy import units as u
import numpy as np
from AeViz.units.units import units


class aerray(np.ndarray):
    
    def __new__(cls,
                value,
                unit=None,
                name=None,
                label=None,
                cmap=None,
                limits=None):
        obj = np.asarray(value).view(cls)
        obj.unit = unit if unit else u.dimensionless_unscaled
        obj.name = name
        obj.label = label
        obj.cmap = cmap
        obj.limits = limits
        return obj
    
    def __array_finalize__(self, obj):
        if obj is None: return
        self.unit = getattr(obj, 'unit', u.dimensionless_unscaled)
        self.name = getattr(obj, 'name', None)
        self.label = getattr(obj, 'label', None)
        self.cmap = getattr(obj, 'cmap', None)
        self.limits = getattr(obj, 'limits', None)
        
    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """ Handle NumPy functions like np.sin, np.exp. """
        input_values = []
        for i in inputs:
            if isinstance(i, aerray):
                old_unit = i.unit
                if ufunc in [np.sin, np.cos, np.tan] and not \
                    old_unit.is_equivalent(u.radian):
                    raise ValueError(f"Trigonometric functions require "
                                     "dimensionless quantities or angles, "
                                     "got unit {old_unit}")
                elif ufunc in [np.sin, np.cos, np.tan]:
                    i = i.to(u.radian)
                input_values.append(i.view(np.ndarray))
            else:
                input_values.append(i)
                    
        result = getattr(ufunc, method)(*input_values, **kwargs)

        if ufunc in [np.sin, np.cos, np.tan]:  # Trig functions return dimensionless values
            new_unit = u.dimensionless_unscaled
        elif ufunc in [np.sqrt, np.exp, np.log]:  # sqrt, exp, log keep unit
            new_unit = u.dimensionless_unscaled
        else:
            new_unit = old_unit
        
        return aerray(result, unit=new_unit, name=None, label=None, cmap=None,
                      limits=None)
    
    def __mul__(self, other):
        """ Handle multiplication with Astropy units and scalars. """
        if isinstance(other, u.UnitBase):  # Handle aerray * unit
            try:
                conv = other.to(self.unit)
                other = self.unit
            except:
                conv = 1
            return aerray(self.value * conv, unit=self.unit * other,
                          name=self.name, label=self.label, cmap=self.cmap,
                          limits=self.limits)
        elif isinstance(other, aerray): #Handle aerray * aerray
            try:
                other = other.to(self.unit)
            except:
                pass
            return aerray(self.value * other.value, unit=self.unit * other.unit)
        elif isinstance(other, (int, float, np.ndarray)):  # Handle aerray * scalar
            return aerray(self.value * other, unit=self.unit, name=self.name,
                          label=self.label, cmap=self.cmap, limits=self.limits)
            
        return NotImplemented
    
    def __rmul__(self, other):
        return self.__mul__(other)
    
    def __truediv__(self, other):
        """ Handle division with scalars, units, and other `arrrray` objects."""
        if isinstance(other, u.UnitBase):  # Unit division (removes unit)
            return aerray(self.value, unit=self.unit / other, name=self.name,
                          label=self.label, cmap=self.cmap, limits=self.limits)
        elif isinstance(other, aerray):  # Array division (unit-aware)
            return aerray(self.value / other.value, unit=self.unit / other.unit)
        elif isinstance(other, (int, float, np.ndarray)):  # Scalar division
            return aerray(self.value / other, unit=self.unit, name=self.name,
                          label=self.label, cmap=self.cmap, limits=self.limits)
        return NotImplemented

    def __rtruediv__(self, other):
        """ Handle division when `other` is on the left: `scalar / arrrray`. """
        if isinstance(other, u.UnitBase):  # Unit division (removes unit)
            return aerray(1 / self.value, unit=other / self.unit, name=self.name,
                          label=self.label, cmap=self.cmap, limits=self.limits)
        if isinstance(other, (int, float, np.ndarray)):  # Scalar divided by arrrray
            return aerray(other / self.value,
                          unit=u.dimensionless_unscaled / self.unit,
                          name=self.name,
                          label=self.label, cmap=self.cmap, limits=self.limits)

        return NotImplemented


    def __repr__(self):
        return f"aerray({super().__repr__()}, unit={self.unit}, " + \
            f"name={self.name})"

        
    @property
    def value(self):
        """ Returns the plain NumPy array without units. """
        return np.asarray(self)
    
    
    def to(self, unit):
        convererersion = self.unit.to(unit)
        return aerray(self.value * convererersion, unit=unit, name=self.name,
                      label=self.label, cmap=self.cmap, limits=self.limits)
        