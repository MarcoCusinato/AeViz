from __future__ import annotations
from typing import TYPE_CHECKING
from AeViz.units import u
from astropy.units import UnitBase, CompositeUnit, Unit
from astropy.utils.compat import COPY_IF_NEEDED
import numpy as np
if TYPE_CHECKING:
    from AeViz.units.aeseries import aeseries

def apply_monkey_patch():
    UnitBase.__mul__ = aemul
    UnitBase.__rmul__ = aermul
    UnitBase.__truediv__ = aetruediv
    UnitBase.__rtruediv__ = aertruediv
    UnitBase.__rlshift__ = aerlshift

def remove_monkey_patch():
    UnitBase.__mul__      = UnitBase._original_mul
    UnitBase.__rmul__     = UnitBase._original_rmul
    UnitBase.__truediv__  = UnitBase._original_truediv
    UnitBase.__rtruediv__ = UnitBase._original_rtruediv
    UnitBase.__rlshift__  = UnitBase._original_rlshift

## Save original multiplication
if not hasattr(UnitBase, "_original_mul"):
    UnitBase._original_mul = UnitBase.__mul__

if not hasattr(UnitBase, "_original_rmul"):
    UnitBase._original_rmul = UnitBase.__rmul__

if not hasattr(UnitBase, "_original_truediv"):
    UnitBase._original_truediv = UnitBase.__truediv__

if not hasattr(UnitBase, "_original_rtruediv"):
    UnitBase._original_rtruediv = UnitBase.__rtruediv__

if not hasattr(UnitBase, "_original_rlshift"):
    UnitBase._original_rlshift = UnitBase.__rlshift__

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
                if self.name == other.name and self.label == other.label:
                    nm, lb, cm, lm, lg = self.name, self.label, self.cmap,\
                        self.limits, self.log
                else:
                    nm, lb, cm, lm, lg = None, None, None, None, False
                return aerray(self.value + other_converted.value,
                              self.unit, nm, lb, cm, lm, lg)
            else:
                raise ValueError(f"Cannot add incompatible units {self.unit}" \
                                 f" and {other.unit}")
                
        elif isinstance(other, (int, float, np.ndarray))  and \
                        self.unit == u.dimensionless_unscaled:
            return aerray(self.value + other, unit=self.unit)
        elif type(other).__name__ == "aeseries":
            return other.__add__(self)
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
        
        elif type(other).__name__ == "aeseries":
            raise TypeError("In place addition does not work between aerray and aeseries, try the opposite ;).")

        raise TypeError("In-place addition only works between compatible"\
                        " units or dimensionless values.")

    def __sub__(self, other):
        """
        Handle subtraction like addition but with `-`.
        """
        if isinstance(other, aerray):
            if self.unit.is_equivalent(other.unit):  # Convert if necessary
                other_converted = other.to(self.unit)
                if self.name == other.name and self.label == other.label:
                    nm, lb, cm, lm, lg = self.name, self.label, self.cmap,\
                        self.limits, self.log
                else:
                    nm, lb, cm, lm, lg = None, None, None, None, False
                return aerray(self.value - other_converted.value,
                              self.unit, nm, lb, cm, lm, lg)
            else:
                raise ValueError(f"Cannot add incompatible units {self.unit}" \
                                 f" and {other.unit}")
                
        elif isinstance(other, (int, float, np.ndarray)) and \
                        self.unit == u.dimensionless_unscaled:
            return aerray(self.value - other, unit=self.unit,
                              name=self.name)
        elif type(other).__name__ == "aeseries":
            raise TypeError("Subtraction only works between aeseries and aerray, not the opposite.")
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
        elif type(other).__name__ == "aeseries":
            raise TypeError("In-place subtraction only works between aeseries " \
                " and aerray, not the opposite.")
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
            # Fixing labels
            if isinstance(other, int):
                new_label = r'%d $\cdot$ %s' % (other, self.label)
            elif isinstance(other, float):
                new_label = r'%.1e $\cdot$ %s' % (other, self.label)
            else:
                new_label = self.label
            new_limits = [self.limits[0] * other, self.limits[1] * other]
            return aerray(self.value * other, unit=self.unit, name=self.name,
                          label=new_label, cmap=self.cmap, \
                          limits=new_limits, log=self.log)
        elif type(other).__name__ == "aeseries":
            return other.__mul__(self)
        
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
                self.value = self.value * conv
            self.unit = self.unit * other
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
        raise TypeError("Not implemented")
    
    def __truediv__(self, other):
        """
        Handle division with scalars, units, and other `aerray` objects.
        """
        if isinstance(other, UnitBase):  # Unit division (removes unit)
            return aerray(self.value, unit=self.unit / other, name=self.name,
                          label=self.label, cmap=self.cmap, limits=self.limits,
                          log=self.log)
        elif isinstance(other, aerray):  # Array division (unit-aware)
            try:
                other = other.to(self.unit)
            except:
                pass
            return aerray(self.value / other.value, unit=self.unit / other.unit)
        elif isinstance(other, (int, float, np.ndarray)):  # Scalar division
            # Fixing labels
            if isinstance(other, int):
                new_label = r'%s / %d' % (self.label, other)
            elif isinstance(other, float):
                new_label = r'%s / %.1e' % (self.label, other)
            else:
                new_label = self.label
            new_limits = [self.limits[0] / other, self.limits[1] / other]
            return aerray(self.value / other, unit=self.unit, name=self.name,
                          label=new_label, cmap=self.cmap, \
                          log=self.log, limits=self.limits / other)
        elif type(other).__name__ == "aeseries":
            raise TypeError("Division does not work between aerray and aeseries.")
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
            if isinstance(other, int):
                new_label = r'%d / %s' % (other, self.label)
            elif isinstance(other, float):
                new_label = r'%.1e / %s' % (other, self.label)
            else:
                new_label = self.label
            new_limits = [other / self.limits[0], other / self.limits[1]]
            return aerray(other / self.value,
                          unit=u.dimensionless_unscaled / self.unit,
                          name=self.name, limits=new_limits, \
                          label=new_label, cmap=self.cmap, \
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
        raise TypeError("Not implemented")
    
    def __pow__(self, exponent):
        """Handle exponentiation while adjusting units correctly."""
        if not np.isscalar(exponent):
            raise TypeError("Exponent must be a scalar value.")

        new_value = self.value ** exponent
        new_unit = self.unit ** exponent  # Properly scale the unit
        new_label = r'%s$^{%d}$' % (self.label, exponent)
        new_limits = [self.limits[0] ** exponent, self.limits[1] ** exponent]

        return aerray(new_value, unit=new_unit, label=new_label, log=self.log, \
                      cmap=self.cmap, limits=new_limits)
    
    ## Operators redefinition
    def __eq__(self, other):
        if type(other).__name__ == "aeseries":
            return other.__eq__(self)
        return np.equal(self, other)
    
    def __ne__(self, other):
        if type(other).__name__ == "aeseries":
            return other.__ne__(self)
        return np.not_equal(self, other)

    def __ge__(self, other):
        if type(other).__name__ == "aeseries":
            return other.__lt__(self)
        return np.greater_equal(self, other)
    
    def __le__(self, other):
        if type(other).__name__ == "aeseries":
            return other.__gt__(self)
        return np.less_equal(self, other)
    
    def __gt__(self, other):
        if type(other).__name__ == "aeseries":
            return other.__le__(self)
        return np.greater(self, other)
    
    def __lt__(self, other):
        if type(other).__name__ == "aeseries":
            return other.__ge__(self)
        return np.less(self, other)
    
    ## Functions redefinition
    
    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """
        Handle NumPy functions like np.sin, np.exp.
        """
        from AeViz.utils.files.string_utils import apply_symbol
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
            new_limits = [np.sqrt(self.limits[0]), np.sqrt(self.limits[1])]
            new_label = apply_symbol(self.label, '\\sqrt')
            #new_label = r'$\sqrt{$%s}$' % self.label
        elif ufunc == np.cbrt:
            new_unit = old_unit ** (1.0 / 3.0)
            new_label = self.label
            #new_label = r'$\cbrt{%s}$' % self.label
        else:
            new_unit = old_unit
        
        return aerray(result, unit=new_unit, label=new_label, log=self.log, \
                      cmap=self.cmap, limits=new_limits)
    
    def __array_function__(self, func, types, args, kwargs):
        """Intercept NumPy functions."""
        if func == np.concatenate:
            return aerray._concatenate(*args, **kwargs)
        elif func == np.moveaxis:
            return aerray._moveaxis(*args, **kwargs)
        elif func in [np.argmax, np.nanargmax, np.argmin, np.nanargmin, 
                      np.argwhere, np.nonzero, np.flatnonzero, np.count_nonzero]:
            return aerray._index_functions(func, *args, **kwargs)
        elif func == np.where:
            return aerray._where(*args, **kwargs)
        elif func == np.interp:
            return aerray._interp(*args, **kwargs)
        elif func == np.meshgrid:
            return self._meshgrid(*args, **kwargs)
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
            arrays = [arr.to(arrays[0].unit) for arr in arrays]
        # Concatenate raw values and return a new `aerray`
        concatenated_values = np.concatenate([[arr.value]
                                              if arr.ndim == 0 else
                                              arr.value
                                              for arr in arrays], axis=axis)
        return aerray(concatenated_values, arrays[0].unit, arrays[0].name,
                      arrays[0].label, arrays[0].cmap, arrays[0].limits,
                      arrays[0].log)
    
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
    def _where(condition, x=None, y=None):
        if x is None and y is None:
            return np.where(condition)
        if isinstance(x, aerray) and isinstance(y, aerray):
            if x.unit != y.unit:
                raise ValueError(f"Cannot mix units: {x.unit} and {y.unit}")
            return aerray(np.where(condition, x.value, y.value), unit=x.unit)
        elif isinstance(x, aerray):
            return aerray(np.where(condition, x.value, y), unit=x.unit,
                          name=x.name, label=x.label, cmap=x.cmap,
                          limits=x.limits, log=x.log)
        elif isinstance(y, aerray):
            return aerray(np.where(condition, x, y.value), unit=y.unit,
                          name=y.name, label=y.label, cmap=y.cmap,
                          limits=y.limits, log=y.log)
        return NotImplemented

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
    def _interp(x, xp, yp):
        
        assert isinstance(x, aerray), "Works only with aerray"
        return aerray(np.interp(x.value, xp.value, yp.value), unit=yp.unit,
                      name=yp.name, label=yp.label, cmap=yp.cmap,
                      limits=yp.limits, log=yp.log)
    
    @staticmethod
    def _meshgrid(*arrays, **kwargs):
        """Custom handling for np.meshgrid."""
        # Ensure all inputs are `aerray` and have compatible units
        assert all(isinstance(arr, aerray) for arr in arrays), "All inputs must be aerray."
        # Compute meshgrid on raw numpy arrays
        values = [arr.value for arr in arrays]
        grids = np.meshgrid(*values, **kwargs)

        # Convert back to aerray with correct unit
        return tuple(aerray(grid, arr.unit, arr.name, arr.label,
                            arr.cmap, arr.limits, arr.log) for (grid, arr)
                            in zip(grids, arrays))


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
        conversion = self.unit.to(unit)
        return aerray(self.value * conversion, unit=unit, name=self.name,
                      label=self.label, cmap=self.cmap, limits=self.limits,
                      log=self.log)

    def max(self, axis=None):
        return aerray(np.max(self.value, axis), unit=self.unit, name=self.name,
                      label=self.label, cmap=self.cmap, limits=self.limits,
                      log=self.log)
    
    def min(self, axis=None): 
        return aerray(np.min(self.value, axis), unit=self.unit, name=self.name,
                      label=self.label, cmap=self.cmap, limits=self.limits,
                      log=self.log)
    
    def mean(self, axis=None):
        return aerray(np.mean(self.value, axis), unit=self.unit, name=self.name,
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
        return UnitBase._original_truediv

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
        return UnitBase._original_rtruediv

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
        return UnitBase._original_mul

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
        return UnitBase._original_rmul

def aerlshift(self, m):
    try:
        return aerray(m, self, copy=COPY_IF_NEEDED, subok=True)
    except Exception:
        if isinstance(m, np.ndarray):
            raise
        return UnitBase._original_rlshift


