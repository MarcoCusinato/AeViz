from AeViz.units import u
from astropy.units import UnitBase, CompositeUnit, Unit
from astropy.utils.compat import COPY_IF_NEEDED
import numpy as np



class aerray(np.ndarray):
    
    def __new__(cls,
                value,
                unit=None,
                name=None,
                label=None,
                cmap=None,
                limits=None,
                log=False):
        obj = np.asarray(value).view(cls)
        obj.unit = unit if unit else u.dimensionless_unscaled
        obj.name = name
        obj.label = label
        obj.cmap = cmap
        obj.limits = limits
        obj.log = log
        return obj
    
    def __array_finalize__(self, obj):
        if obj is None: return
        self.unit = getattr(obj, 'unit', u.dimensionless_unscaled)
        self.name = getattr(obj, 'name', None)
        self.label = getattr(obj, 'label', None)
        self.cmap = getattr(obj, 'cmap', None)
        self.limits = getattr(obj, 'limits', None)
        self.log = getattr(obj, 'log', False)
        
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
        elif ufunc in [np.exp, np.log]:
            new_unit = u.dimensionless_unscaled
        elif ufunc == np.sqrt:
            new_unit = old_unit ** 0.5
        else:
            raise NotImplementedError
        
        return aerray(result, unit=new_unit)
    
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
                          limits=self.limits, log=self.log)
        elif isinstance(other, aerray): #Handle aerray * aerray
            try:
                other = other.to(self.unit)
            except:
                pass
            return aerray(self.value * other.value, unit=self.unit * other.unit)
        elif isinstance(other, (int, float, np.ndarray)):
            return aerray(self.value * other, unit=self.unit, name=self.name,
                          label=self.label, cmap=self.cmap, limits=self.limits,
                          log=self.log)
            
        return NotImplemented
    
    def __rmul__(self, other):
        return self.__mul__(other)
    
    def __imul__(self, other):
        """
        In-place multliplication (*=) with unit handling.
        """
        if isinstance(other, u.UnitBase):  # Handle aerray * unit
            try:
                conv = other.to(self.unit)
                other = self.unit
            except:
                conv = 1
            if self.ndim == 0:
                self.fill(self.item() * conv)  # Handle 0-D case
            else:
                self[:] *= conv
            self.unit *= other
            return self    
        elif isinstance(other, aerray): #Handle aerray * aerray
            try:
                other = other.to(self.unit)
            except:
                pass
            if self.ndim == 0:
                self.fill(self.item() * other.value)
            else:
                self[:] *= other.value
            self.unit *= other.unit
            return self
        elif isinstance(other, (int, float, np.ndarray)):
            if self.ndim == 0:
                self.fill(self.item() * other)
            else:
                self[:] *= other
            return self

        raise TypeError("In-place multilpication  only works between compatible"\
                        " units or dimensionless values.")
    
    def __truediv__(self, other):
        """
        Handle division with scalars, units, and other `aerray` objects.
        """
        if isinstance(other, u.UnitBase):  # Unit division (removes unit)
            return aerray(self.value, unit=self.unit / other, name=self.name,
                          label=self.label, cmap=self.cmap, limits=self.limits,
                          log=self.log)
        elif isinstance(other, aerray):  # Array division (unit-aware)
            return aerray(self.value / other.value, unit=self.unit / other.unit)
        elif isinstance(other, (int, float, np.ndarray)):  # Scalar division
            return aerray(self.value / other, unit=self.unit, name=self.name,
                          label=self.label, cmap=self.cmap, limits=self.limits,
                          log=self.log)
        return NotImplemented

    def __rtruediv__(self, other):
        """
        Handle division when `other` is on the left: `scalar / aerray`.
        """
        if isinstance(other, u.UnitBase):  # Unit division (removes unit)
            return aerray(1 / self.value, unit=other / self.unit, name=self.name,
                          label=self.label, cmap=self.cmap, limits=self.limits,
                          log=self.log)
        if isinstance(other, (int, float, np.ndarray)):  # Scalar divided by aerray
            return aerray(other / self.value,
                          unit=u.dimensionless_unscaled / self.unit,
                          name=self.name,
                          label=self.label, cmap=self.cmap, limits=self.limits,
                          log=self.log)

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

    def __pow__(self, exponent):
        """Handle exponentiation while adjusting units correctly."""
        if not np.isscalar(exponent):
            raise TypeError("Exponent must be a scalar value.")

        new_value = self.value ** exponent
        new_unit = self.unit ** exponent  # Properly scale the unit

        return aerray(new_value, new_unit)
    
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
                      label=self.label, cmap=self.cmap, limits=self.limits,
                      log=self.log)

    def max(self):
        return aerray(np.max(self.value), unit=self.unit, name=self.name,
                      label=self.label, cmap=self.cmap, limits=self.limits,
                      log=self.log)
    
    def min(self): 
        return aerray(np.min(self.value), unit=self.unit, name=self.name,
                      label=self.label, cmap=self.cmap, limits=self.limits,
                      log=self.log)
    
    def set(self, name=None, label=None, cmap=None, limits=None, log=False):
        if name is not None:
            self.name = name
        if label is not None:
            self.label = label
        if cmap is not None:
            self.cmap = cmap
        if limits is not None:
            self.limits = limits
        if log is not None:
            self.log = log



### Monkey patches for multiplication and division
def aetruediv(self, m):
    if isinstance(m, (bytes, str)):
        m = Unit(m)

    if isinstance(m, UnitBase):
        if m.is_unity():
            return self
        return CompositeUnit(1, [self, m], [1, -1], _error_check=False)

    try:
        return aerray(1, self) / m
    except TypeError:
        return NotImplemented

def aertruediv(self, m):
    if isinstance(m, (bytes, str)):
        return Unit(m) / self

    try:
        # Cannot handle this as Unit.  Here, m cannot be a Quantity,
        # so we make it into one, fasttracking when it does not have a
        # unit, for the common case of <array> / <unit>.

        if hasattr(m, "unit"):
            result = aerray(m)
            result /= self
            return result
        else:
            return aerray(m, self ** (-1))
    except TypeError:
        if isinstance(m, np.ndarray):
            raise
        return NotImplemented

def aemul(self, m):
    if isinstance(m, (bytes, str)):
        m = Unit(m)

    if isinstance(m, UnitBase):
        if m.is_unity():
            return self
        elif self.is_unity():
            return m
        return CompositeUnit(1, [self, m], [1, 1], _error_check=False)

    # Cannot handle this as Unit, re-try as Quantity.
    try:
        return aerray(1, unit=self) * m
    except TypeError:
        return NotImplemented

def aermul(self, m):
    if isinstance(m, (bytes, str)):
        return Unit(m) * self

    # Cannot handle this as Unit.  Here, m cannot be a Quantity,
    # so we make it into one, fasttracking when it does not have a unit
    # for the common case of <array> * <unit>.
    try:
        if hasattr(m, "unit"):
            result = aerray(m)
            result *= self
            return result
        else:
            return aerray(m, unit=self)
    except TypeError:
        if isinstance(m, np.ndarray):
            raise
        return NotImplemented

def aerlshift(self, m):
    try:
        return aerray(m, self, copy=COPY_IF_NEEDED, subok=True)
    except Exception:
        if isinstance(m, np.ndarray):
            raise
        return NotImplemented


u.UnitBase.__mul__ = aemul
u.UnitBase.__rmul__ = aermul
u.UnitBase.__truediv__ = aetruediv
u.UnitBase.__rtruediv__ = aertruediv
u.UnitBase.__rlshift__ = aerlshift