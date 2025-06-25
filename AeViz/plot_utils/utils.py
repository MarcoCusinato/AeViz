import numpy as np
from AeViz.units import u, aerray

def get_1Dhist_data(xlims, nbins, xscale, xdata, ydata):
       assert xlims is None or len(xlims) == 2, "Too many limits"

       if xlims is None:
              xlims = xdata.limits
       if nbins is None:
              nbins = 20

       if isinstance(xlims[0], aerray) and isinstance(xlims[1], aerray):
              xlims = [xlims[0].value, xlims[1].value]
       if xscale == 'linear':
              bins = np.linspace(xlims[0], xlims[1], nbins)
       else:
              raise NotImplementedError("Logarithmic bins are not implemented yet")
       ybinned, bin_edges = np.histogram(xdata.value, bins=bins,
                                           weights=ydata.value)
       bin_centers = (bin_edges[:-1] + bin_edges[1:]) * 0.5

       ybinned = aerray(ybinned, ydata.unit, ydata.name, ydata.label,
                        None, [ybinned.min(), ybinned.max()], False)
       bin_centers = aerray(bin_centers, xdata.unit, xdata.name, xdata.label,
                            None, xlims, False)
       return bin_centers, ybinned, np.diff(bin_edges)


    