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
from AeViz.utils.GW_utils import Qdot_timeseries
from AeViz.utils.radii_utils import shock_radius, shock_radius_3D

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

Plot.plot1D('PGAS', 'radius',  12, None)

Plot.plot1D('ENTR', 'radius',  12, None)
#Plot.plot2DwithPar('V')

input()
"""

#sim = Simulation('s16.5-SFHo-1.0omg-5e+08-1e+09B-timestep0.1ms', '/almacen/marco/Simulations/sn2d/s16.5-SW14/')
sim = Simulation('A26-1', '/home/marco/Escritorio/Simulations/martin/sn3d/Bonn/')
sim.shock_radius()
h =Qdot_timeseries(sim, True)
plt.plot(h[0], h[2][0])
plt.show()

radius = sim.cell.radius(sim.ghost)
P = sim.gas_pressure('h00003300.h5')
vx = sim.radial_velocity('h00003300.h5')
s = sim.entropy('h00003300.h5')

shock = shock_radius(sim, 'h00003300.h5')
print(shock[ 100,  7])
#print(np.any(shock==np.nan))
fig, axs = plt.subplots(1, 3, figsize=(10, 5), sharex=True)
#axs[0].plot(sim.cell.theta(sim.ghost), shock)

axs[0].plot(sim.cell.radius(sim.ghost), P[ 100,  7, :])
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].axvline(shock[ 100,  7], color='k')
ax1 =axs[0].twinx()
ax1.plot(sim.cell.radius(sim.ghost), IDL_derivative(radius, P[ 100,  7, :])/P[ 100,  7, :] * radius, color='b')
axs[1].plot(sim.cell.radius(sim.ghost), vx[ 100,  7, :])
axs[1].set_xscale('symlog')
ax2=axs[1].twinx()
ax2.plot(sim.cell.radius(sim.ghost), IDL_derivative(radius, vx[ 100,  7, :])/np.abs(vx[ 100,  7, :]) * radius, color='b')
axs[2].plot(sim.cell.radius(sim.ghost), s[ 100,  7, :])
"""

sim.innercore_radius()
sim.PNS_nucleus_radius()
"""
#plt.show()
#data = sim.AE220()
#fig, axs = plt.subplots(3, 1, sharex=True)
#axs[0].plot(data[1], data[3])
#axs[1].plot(data[1], data[4])
#axs[2].plot(data[1], data[5])
plt.show()
#print(time.shape, radius.shape, pr.shape)