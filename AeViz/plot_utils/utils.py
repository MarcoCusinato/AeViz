import numpy as np
from AeViz.units import u, aerray
import re

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

def get_2Dhist_data(xlims, ylims, nbins, xscale, yscale, xdata, ydata, cdata):
       assert nbins is None or len(nbins) == 2, "Weird number of bins"
       if xlims is None:
              xlims = xdata.limits
       if ylims is None:
              ylims = ydata.limits
       if nbins is None:
              nxbins, nybins = 100, 100
       else:
              nxbins, nybins = nbins
       if isinstance(xlims[0], aerray) and isinstance(xlims[1], aerray):
              xlims = [xlims[0].value, xlims[1].value]
       if isinstance(ylims[0], aerray) and isinstance(ylims[1], aerray):
              ylims = [ylims[0].value, ylims[1].value]
       if xscale == 'linear':
              xbins = np.linspace(xlims[0], xlims[1], nxbins)
       else:
              raise NotImplementedError("Logarithmic bins are not implemented yet")
       if yscale == 'linear':
              ybins = np.linspace(ylims[0], ylims[1], nybins)
       else:
              raise NotImplementedError("Logarithmic bins are not implemented yet")
# Compute 2D histogram weighted by mass
       bins, xbinned, ybinned = np.histogram2d(
              xdata.value.flatten(), ydata.value.flatten(), bins=[xbins, ybins],
              weights=cdata.value.flatten()
              )
       bins = aerray(bins.T, cdata.unit, cdata.name, cdata.label,
                            cdata.cmap, [bins.max(), bins.min()], False)
       xbinned = aerray(xbinned, xdata.unit, xdata.name, xdata.label,
                            None, xlims, False)
       ybinned = aerray(ybinned, ydata.unit, ydata.name, ydata.label,
                        None, ylims, False)
       return xbinned, ybinned, bins

    