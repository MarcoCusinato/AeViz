from AeViz.units import u
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
        """
        Handle NumPy functions like np.sin, np.exp.
        """
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

        if ufunc in [np.sin, np.cos, np.tan]:
            new_unit = u.dimensionless_unscaled
        elif ufunc in [np.sqrt, np.exp, np.log]:
            new_unit = u.dimensionless_unscaled
        else:
            new_unit = old_unit
        
        return aerray(result, unit=new_unit, name=None, label=None, cmap=None,
                      limits=None)
    
    def __mul__(self, other):
        """
        Handle multiplication with Astropy units and scalars.
        """
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
        elif isinstance(other, (int, float, np.ndarray)):
            return aerray(self.value * other, unit=self.unit, name=self.name,
                          label=self.label, cmap=self.cmap, limits=self.limits)
            
        return NotImplemented
    
    def __rmul__(self, other):
        return self.__mul__(other)
    
    def __truediv__(self, other):
        """
        Handle division with scalars, units, and other `aerray` objects.
        """
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
        """
        Handle division when `other` is on the left: `scalar / aerray`.
        """
        if isinstance(other, u.UnitBase):  # Unit division (removes unit)
            return aerray(1 / self.value, unit=other / self.unit, name=self.name,
                          label=self.label, cmap=self.cmap, limits=self.limits)
        if isinstance(other, (int, float, np.ndarray)):  # Scalar divided by aerray
            return aerray(other / self.value,
                          unit=u.dimensionless_unscaled / self.unit,
                          name=self.name,
                          label=self.label, cmap=self.cmap, limits=self.limits)

        return NotImplemented


    def __repr__(self):
        return f"aerray({self.value}, unit={self.unit}, " + \
            f"name={self.name})"
    
    def __add__(self, other):
        """
        Handle addition with another aerray or quantity.
        """
        if isinstance(other, aerray):
            if self.unit.is_equivalent(other.unit):  # Convert if necessary
                other_converted = other.to(self.unit)
                return aerray(self.value + other_converted.value,
                              unit=self.unit)
            else:
                raise ValueError(f"Cannot add incompatible units {self.unit}" \
                                 f" and {other.unit}")
                
        elif isinstance(other, (int, float, np.ndarray))  and \
                        self.unit == u.dimensionless_unscaled:
            return aerray(self.value + other, unit=self.unit)
        raise TypeError("Addition only works between compatible units" \
                        " or dimensionless values.")

    def __sub__(self, other):
        """
        Handle subtraction like addition but with `-`.
        """
        if isinstance(other, aerray):
            if self.unit.is_equivalent(other.unit):  # Convert if necessary
                other_converted = other.to(self.unit)
                return aerray(self.value - other_converted.value,
                              unit=self.unit)
            else:
                raise ValueError(f"Cannot add incompatible units {self.unit}" \
                                 f" and {other.unit}")
                
        elif isinstance(other, (int, float, np.ndarray)) and \
                        self.unit == u.dimensionless_unscaled:
            return aerray(self.value - other, unit=self.unit,
                              name=self.name)
        raise TypeError("Subtraction only works between compatible units" \
                        " or dimensionless values.")

    def __radd__(self, other):
        return self.__add__(other)

    def __rsub__(self, other):
        """
        Handle `other - self` correctly.
        """
        if isinstance(other, (int, float, np.ndarray)) and \
            self.unit == u.dimensionless_unscaled:
            return aerray(other - self.value, unit=self.unit)
        return NotImplemented
    
    def __iadd__(self, other):
        """
        In-place addition (+=) with unit handling.
        """
        if isinstance(other, aerray):
            if self.unit.is_equivalent(other.unit):
                other_converted = other.to(self.unit)
                if self.ndim == 0:
                    self.fill(self.item() + other_converted.item())
                else:
                    self[:] += other_converted.value
                return self
            else:
                raise ValueError(f"Cannot add incompatible units {self.unit} "\
                                 f"and {other.unit}")

        elif isinstance(other, (int, float, np.ndarray)) and \
            self.unit == u.dimensionless_unscaled:
            if self.ndim == 0:
                self.fill(self.item() + other)  # Handle 0-D case
            else:
                self[:] += other
            return self

        raise TypeError("In-place addition only works between compatible"\
                        " units or dimensionless values.")

    def __isub__(self, other):
        """
        In-place subtraction (-=) with unit handling.
        """
        if isinstance(other, aerray):
            if self.unit.is_equivalent(other.unit):
                other_converted = other.to(self.unit)
                if self.ndim == 0:
                    self.fill(self.item() - other_converted.item())  # Handle 0-D case
                else:
                    self[:] -= other_converted.value
                return self
            else:
                raise ValueError(f"Cannot add incompatible units {self.unit} "\
                                 f"and {other.unit}")

        elif isinstance(other, (int, float, np.ndarray)) and \
            self.unit == u.dimensionless_unscaled:
            if self.ndim == 0:
                self.fill(self.item() - other)  # Handle 0-D case
            else:
                self[:] -= other
            return self

        raise TypeError("In-place subtraction only works between compatible"\
                        " units or dimensionless values.")


    @property
    def value(self):
        """
        Returns the plain NumPy array without units.
        """
        return np.asarray(self)
    
    @value.setter
    def value(self, new_value):
        """
        Set the underlying array data dynamically.
        """
        new_value = np.asarray(new_value)
        if self.ndim == 0:
            self.fill(new_value)
        else:
            if new_value.shape != self.shape:
                raise ValueError(f"Shape mismatch: expected {self.shape}, " \
                                 f"got {new_value.shape}")
            self[:] = new_value
    
    def to(self, unit):
        convererersion = self.unit.to(unit)
        return aerray(self.value * convererersion, unit=unit, name=self.name,
                      label=self.label, cmap=self.cmap, limits=self.limits)

    def max(self):
        return aerray(np.max(self.value), unit=self.unit, name=self.name,
                      label=self.label, cmap=self.cmap, limits=self.limits)
    
    def min(self): 
        return aerray(np.min(self.value), unit=self.unit, name=self.name,
                      label=self.label, cmap=self.cmap, limits=self.limits)
    
    def set(self, name=None, label=None, cmap=None, limits=None):
        self.name = name
        self.label = label
        self.cmap = cmap
        self.limits = limits

# --- Monkey-Patch Astropy Units ---
def ae_multiply(self, other):
    """
    Redirect multiplication with an Astropy unit to create an `aerray`.
    """
    return aerray(other, unit=self)

u.UnitBase.__mul__ = ae_multiply
u.UnitBase.__rmul__ = ae_multiply