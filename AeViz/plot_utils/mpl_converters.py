from AeViz.units import aerray, aeseries, u
import numpy as np
from  matplotlib import ticker, units
from AeViz.utils.files.string_utils import merge_strings

def rad_fn(x, pos=None):
        n = int((x / np.pi) * 2.0 + 0.25)
        if n == 0:
            return "0"
        elif n == 1:
            return "π/2"
        elif n == 2:
            return "π"
        elif n % 2 == 0:
            return f"{n // 2}π"
        else:
            return f"{n}π/2"


class MplaerrayConverter(units.ConversionInterface):
    """Matplotlib Converter for aerray objects with unit support."""

    @staticmethod
    def convert(value, unit, axis):
        """
        Convert aerray to numeric data for plotting.
        """
        if unit is not None:
            lab, un = unit
            if isinstance(value, aerray):
                return value.to(un).value
            elif isinstance(value, list) and isinstance(value[0], aerray):
                return np.array([v.to(un).value if isinstance(v, aerray) 
                                else v for v in value])
        else:
            if isinstance(value, aerray):
                return value.value
            elif isinstance(value, list) and isinstance(value[0], aerray):
                return np.array([v.value if isinstance(v, aerray) 
                                else v for v in value])
        return value

    @staticmethod
    def axisinfo(unit, axis):
        """
        Provide axis information including labels with units.
        """
        if isinstance(unit, tuple):
            qtlabel, un = unit
        else:
            return
        ax_un = axis.get_units()
        if ax_un is None:
            axis.set_units(un)
        if un == u.dimensionless_unscaled:
            labl = qtlabel
        labl = merge_strings(qtlabel, f' [{un.to_string(format='latex')}]')
        if un == u.radian:
            return units.AxisInfo(
                majloc=ticker.MultipleLocator(base=np.pi / 2),
                majfmt=ticker.FuncFormatter(rad_fn),
                label=labl,
            )
        elif un == u.degree:
            return units.AxisInfo(
                majloc=ticker.AutoLocator(),
                majfmt=ticker.FormatStrFormatter("%g°"),
                label=labl,
            )
        return units.AxisInfo(
            label=labl,
            majloc=ticker.AutoLocator(),
            majfmt=ticker.ScalarFormatter()
        )

    @staticmethod
    def default_units(value, axis):
        """Return a tuple (label, units) for axis labeling."""
        if isinstance(value, aerray):
            return (value.label, value.unit)  # Store both label and units
        return None