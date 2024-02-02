import matplotlib.pyplot as plt
import numpy as np
from AeViz.quantities_plotting.plotting import Plotting
from AeViz.plot_utils.plotting_utils import PlottingUtils
import AeViz.utils.radii_utils
from AeViz.simulation.simulation import Simulation
from AeViz.utils.math_utils import IDL_derivative
from scipy.interpolate import griddata
from AeViz.utils.load_save_radii_utils import calculate_radius
from AeViz.utils.profiles import calculate_profile

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


Plot = Plotting()

Plot.Load('/almacen/marco/Simulations/sn2d/s16.5-SW14/s16.5-SFHo-1.0omg-5e+08-1e+09B/outp-hdf/h00045000.h5')

Plot.plot1D('PGAS', 'radius', 32, None)

Plot.plot1D('ENTR', 'radius', 32, None)
#Plot.plot2DwithPar('V')

input()
"""

sim = Simulation('s16.5-SFHo-2.0omg-5e+08-1e+09B', '/almacen/marco/Simulations/sn2d/s16.5-SW14/')
# sim = Simulation('A26-1', '/home/marco/Escritorio/Simulations/martin/sn3d/Bonn/')
"""
radius = sim.cell.radius(sim.ghost)
P = sim.gas_pressure('h00045000.h5')
vx = sim.radial_velocity('h00045000.h5')

shock = shock_radius(sim, 'h00045000.h5')
print(np.any(shock==np.nan))
fig, axs = plt.subplots(1, 2, figsize=(10, 5))
axs[0].plot(sim.cell.theta(sim.ghost), shock)

axs[0].plot(sim.cell.radius(sim.ghost), P[32,...])
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].axvline(shock[32], color='k')
ax1 =axs[0].twinx()
ax1.plot(sim.cell.radius(sim.ghost), IDL_derivative(radius, P[32,:])/P[32, :] * radius, color='b')
axs[1].plot(sim.cell.radius(sim.ghost), vx[32,...])
axs[1].set_xscale('symlog')
ax2=axs[1].twinx()
ax2.plot(sim.cell.radius(sim.ghost), IDL_derivative(radius, vx[32,:])/np.abs(vx[32, :]) * radius, color='b')
"""



#plt.show()
time, radius, pr = sim.PNS_nucleus_radius()
#print(time.shape, radius.shape, pr.shape)