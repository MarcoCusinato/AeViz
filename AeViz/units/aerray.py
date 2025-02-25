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

    def __getitem__(self, indices):
        """
        Redefine the item getter
        """
        return aerray(self.value[indices], self.unit, self.name, self.label,
                      self.cmap, self.limits, self.log)
    
    def __repr__(self):
        outstring = f"aerray({self.value}"
        for nm in ['unit', 'name', 'label', 'cmap', 'limits', 'log']:
            if getattr(self, nm) is not None:
                outstring += f"\n\t{nm}: {getattr(self, nm)}"
        outstring += ")"
        return outstring  
    
    ## Operation ridefinition
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

    def __radd__(self, other):
        return self.__add__(other)
    
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
                    self.value = self.value + other_converted.value
                return self
            else:
                raise ValueError(f"Cannot add incompatible units {self.unit} "\
                                 f"and {other.unit}")

        elif isinstance(other, (int, float, np.ndarray)) and \
            self.unit == u.dimensionless_unscaled:
            if self.ndim == 0:
                self.fill(self.item() + other)  # Handle 0-D case
            else:
                self.value = self.value + other
            return self

        raise TypeError("In-place addition only works between compatible"\
                        " units or dimensionless values.")

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

    def __rsub__(self, other):
        """
        Handle `other - self` correctly.
        """
        if isinstance(other, (int, float, np.ndarray)) and \
            self.unit == u.dimensionless_unscaled:
            return aerray(other - self.value, unit=self.unit)
        return NotImplemented
    
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
                    self.value = self.value - other_converted.value
                return self
            else:
                raise ValueError(f"Cannot add incompatible units {self.unit} "\
                                 f"and {other.unit}")

        elif isinstance(other, (int, float, np.ndarray)) and \
            self.unit == u.dimensionless_unscaled:
            if self.ndim == 0:
                self.fill(self.item() - other)  # Handle 0-D case
            else:
                self.value = self.value - other
            return self
        raise TypeError("In-place subtraction only works between compatible"\
                        " units or dimensionless values.")

    def __mul__(self, other):
        """
        Handle multiplication with Astropy units and scalars.
        """
        if isinstance(other, UnitBase):  # Handle aerray * unit
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
        if isinstance(other, UnitBase):  # Handle aerray * unit
            try:
                conv = other.to(self.unit)
                other = self.unit
            except:
                conv = 1
            if self.ndim == 0:
                self.fill(self.item() * conv)  # Handle 0-D case
            else:
                self.value = self.__dict__[name] = value * conv
            self.unit = self.unit *  other
            return self    
        elif isinstance(other, aerray): #Handle aerray * aerray
            try:
                other = other.to(self.unit)
            except:
                pass
            if self.ndim == 0:
                self.fill(self.item() * other.value)
            else:
                self.value = self.value * other.value
            self.unit = self.unit * other.unit
            return self
        elif isinstance(other, (int, float, np.ndarray)):
            if self.ndim == 0:
                self.fill(self.item() * other)
            else:
                self.value = self.value * other
            return self
        raise TypeError("In-place multilpication  only works between compatible"\
                        " units or dimensionless values.")
    
    def __truediv__(self, other):
        """
        Handle division with scalars, units, and other `aerray` objects.
        """
        if isinstance(other, UnitBase):  # Unit division (removes unit)
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
        if isinstance(other, UnitBase):  # Unit division (removes unit)
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

    def __itruediv__(self, other):
        """
        In-place division (/=) with unit handling.
        """
        if isinstance(other, UnitBase):  # Handle aerray * unit
            try:
                conv = other.to(self.unit)
                other = self.unit
            except:
                conv = 1
            if self.ndim == 0:
                self.fill(self.item() / conv)  # Handle 0-D case
            else:
                self.value = self.value / conv
            self.unit = self.unit / other
            return self    
        elif isinstance(other, aerray): #Handle aerray * aerray
            try:
                other = other.to(self.unit)
            except:
                pass
            if self.ndim == 0:
                self.fill(self.item() / other.value)
            else:
                self.value = self.value / other.value
            self.unit = self.unit / other.unit
            return self
        elif isinstance(other, (int, float, np.ndarray)):
            if self.ndim == 0:
                self.fill(self.item() / other)
            else:
                self.value = self.value /other
            return self
        raise TypeError("In-place division  only works between compatible"\
                        " units or dimensionless values.")
    
    def __pow__(self, exponent):
        """Handle exponentiation while adjusting units correctly."""
        if not np.isscalar(exponent):
            raise TypeError("Exponent must be a scalar value.")

        new_value = self.value ** exponent
        new_unit = self.unit ** exponent  # Properly scale the unit

        return aerray(new_value, new_unit)
    
    ## Operators redefinition
    def __eq__(self, other):
        return np.equal(self, other)
    
    def __neq__(self, other):
        return np.not_equal(self, other)

    def __ge__(self, other):
        return np.greater_equal(self, other)
    
    def __le__(self, other):
        return np.less_equal(self, other)
    
    def __gt__(self, other):
        return np.greater(self, other)
    
    def __lt__(self, other):
        return np.less(self, other)
    
    def __ne__(self, other):
        return np.not_equal(self, other)
    
    ## Functions redefinition
    
    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """
        Handle NumPy functions like np.sin, np.exp.
        """
        ## First check the comparison functions
        if ufunc in [np.greater, np.greater_equal, np.less, np.less_equal,
                     np.equal, np.not_equal]:
            return self._comparison_functions(ufunc, *inputs)

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
            new_unit = old_unit
        
        return aerray(result, unit=new_unit)
    
    def __array_function__(self, func, types, args, kwargs):
        """Intercept NumPy functions."""
        if func == np.concatenate:
            return aerray._concatenate(*args, **kwargs)
        elif func == np.moveaxis:
            return aerray._moveaxis(*args, **kwargs)
        elif func in [np.argmax, np.nanargmax, np.argmin, np.nanargmin, 
                      np.where, np.argwhere, np.nonzero, np.flatnonzero,
                      np.count_nonzero]:
            return aerray._index_functions(func, *args, **kwargs)
        else:
            return aerray._other_functions(func, *args, **kwargs)
        return NotImplemented

    @staticmethod
    def _concatenate(arrays, axis=0):
        """Custom `concatenate` implementation for `aerray`."""
        # Ensure all elements are `aerray`
        if not all(isinstance(arr, aerray) for arr in arrays):
            raise TypeError("All inputs to concatenate must be aerray instances.")

        # Ensure all units match
        units = {arr.unit for arr in arrays}
        if len(units) > 1:
            raise ValueError(f"Cannot concatenate aerrays with different units: {units}")

        # Concatenate raw values and return a new `aerray`
        concatenated_values = np.concatenate([arr.view(np.ndarray) for arr in arrays], axis=axis)
        return aerray(concatenated_values, unit=arrays[0].unit)
    
    @staticmethod
    def _moveaxis(array, source, destination):
        """Custom moveaxis fiunction for the aerray"""
        if not isinstance(array, aerray):
            raise TypeError("Must be an aerray to move the axis")
        moved_axis = np.moveaxis(array.value, source, destination)
        return aerray(moved_axis, array.unit, array.name, array.label,
                      array.cmap, array.limits, array.log)
        
    @staticmethod
    def _index_functions(function, *args, **kwargs):
        arr = args[0]
        assert isinstance(arr, aerray), "Works only with aerray"
        return function(arr.value, *args[1:], **kwargs)
    
    @staticmethod
    def _comparison_functions(function, *args):
        assert isinstance(args[0], aerray) or isinstance(args[1], aerray), \
            "One of the two arguments must be an aerray"
        main = args[0]
        other = args[1]
        if isinstance(other, aerray):
            if not main.unit.is_equivalent(other.unit):
                raise TypeError(f"Cannot compare {main.unit} with {other.unit}")
            other = other.to(main.unit)
            return function(main.value, other.value)
        elif isinstance(other, (int, float, np.ndarray)):
            return function(main.value, other)

    @staticmethod
    def _other_functions(function, *args, **kwargs):
        arr = args[0]
        assert isinstance(arr, aerray), "Works only with aerray"
        return aerray(function(arr.value, *args[1:], **kwargs), unit=arr.unit,
                      label=arr.label, cmap=arr.cmap, limits=arr.limits,
                      log=arr.log)
    
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


UnitBase.__mul__ = aemul
UnitBase.__rmul__ = aermul
UnitBase.__truediv__ = aetruediv
UnitBase.__rtruediv__ = aertruediv
UnitBase.__rlshift__ = aerlshift