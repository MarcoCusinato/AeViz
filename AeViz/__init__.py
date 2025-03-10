import sys
from AeViz.AeViz import AeViz
from AeViz.plot_utils.mpl_converters import MplaerrayConverter, aerray
from matplotlib import units

units.registry[aerray] = MplaerrayConverter()

ae = AeViz()
