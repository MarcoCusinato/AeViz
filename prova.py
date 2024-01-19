import matplotlib.pyplot as plt
import numpy as np
from quantities_plotting.plotting import Plotting
from plot_utils.plotting_utils import PlottingUtils
from pprint import pprint

"""
fig = plt.figure(figsize=(10, 10))
axs = fig.subplot_mosaic()

fig.subplots_adjust(wspace=0.0)


axs["B"].sharex(axs["A"])
axs["B"].sharey(axs["A"])
axs["A"].tick_params(top=False, labeltop=False,
                                          bottom=True, labelbottom=True,
                                          left=True, labelleft=True,
                                          right=False, labelright=False)
axs["B"].tick_params(top=False, labeltop=False,
                                          bottom=True, labelbottom=True,
                                          left=False, labelleft=False,
                                          right=True, labelright=True)

axs["A"].set(xlim=(10, 20), ylim=(10, 20), aspect=1)
axs["B"].set(xlim=(10, 20), ylim=(10, 20), aspect=1)
for ax in axs.values():
    ax.plot(np.linspace(10, 20, 100), np.linspace(10, 20, 100), color='k')
print(axs["A"].get_xlim())
axs["A"].set_xlabel('pippo')
plt.show()
for ax in axs:
    print(ax)
    print(axs[ax].get_xlim())
    print(axs[ax].get_xlabel())
"""

Plot = Plotting()

Plot.Load('/almacen/marco/Simulations/sn2d/s16.5-SW14/s16.5-SFHo-1.0omg-5e+08-1e+09B/outp-hdf/h00045000.h5')

Plot.plot1D('PGAS', 'theta', None, None)
#Plot.plot2DwithPar('V')

input()