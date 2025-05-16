# AeViz

AeViz is a python package to explore data from core collapse supernova simulations made with Aenus-ALCAR ([Obergaulinger, 2008](https://ui.adsabs.harvard.edu/abs/2008PhDT.........1O/abstract); [Just _et al._, 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.453.3386J/abstract); [Just _et al._, 2018](https://ui.adsabs.harvard.edu/abs/2018MNRAS.481.4786J/abstract)). It contains analysis and visualization subroutines for 1, 2 and 3 dimensional simulations.


## Installation
### Package installation
First make sure to have installed Python>=3.8.<br/>
To install the latest AeViz version please clone the repository with
```
git clone https://github.com/MarcoCusinato/AeViz.git
```

Then install the requirements needed to run all the functionalities of AeViz. You can either `pip install` or `conda install` them one by one taking them from this list
```
numpy
scipy
h5py
pandas
f90nml
matplotlib
opencv-python
```
or navigate into the package folder and install them from the `requirements.txt` file
```
cd AeViz
pip install -r requirements.txt
```
After having installed the requirements, we can finally install the main package. If you are not yet in the package folder, navigate into it (`cd AeViz`) and simply run
```
pip install .
```
### First import
After the installation it is time for the first import of the package, after launching Python, let's tipe
```
import AeViz
```
After, this firts import a folder called `Aenus_simulation_postprocessing` will be created in your _home directory_, alongside several subfolders. These subfolders include three labelled as `*D` with `*` being 1, 2, or 3 and containing the postprocessing of the simulations according to their dimension, and a hidden folder called `.utils` containing the path of the simulations.
## Quick start
During AeViz installation, two system-wide scripts have been installed, namely `add_Aenus_sim_paths` and `run_postprocessing`.<br/>
First of all, let's use the former to save some simulation's path as
```
add_Aenus_sim_paths --add-paths /path/to/parent/
```
This will scan the parent folder for simulations. Be careful that if you feed it the simulation path it will not recognize it.

After a simulation path has been saved you can access it by simply using the name of the simulation folder.<br/>
Then we can run the second script with the additional keyword `--plot` as 
```
run_postprocessing --sim-name simulation_name --plot
```
where the `simulation_name` is the name of the simulation folder folder.
By doing this we will have created `h5` postprocessing files as well postprocessing plots insede the folder folder with the corresponding simulation name inside the `Aenus_simulation_postprocessing/DIMENSION` directory.

### Plotting
This is a brief guide on how to read and plot some results. It will refer to a live python session, but this can easily translated into scripts if needed.

#### Time evolution

**Note**: The zero of the time evolution will correspond to bounce time, so that the actual time will always be $t - t_b$.

First of all, let's import what you need:
```
from AeViz import ae
```
This will import the methods.

If you already add the path of your simulation thanks to `add_Aenus_sim_paths` script, you can load your data with
```
ae.Load('name-of-the-simulation')
```
If this is not the case, do not worry! You can always manually set the path when you load the simulation:
```
ae.Load('name-of-the-simulation', simulation_path='/path/to/the/simulation/')
```

In order to plot the time evolution of a quantity (for example, the entropy), just use
```
ae.entropy()
```
This will print the time evolution of the maximum value of the specified quantity. This works also for postprocessed quantities, such as the shock radius.
The command
```
ae.shock_radius('avg') # 'avg' (default), 'min', 'max'
```
will plot the time evolution of the average shock radius. If you want the minimum of maximum radius, just use `'min'` or `'max'` options instead.
Please note that if `run_postprocessing` was not previously run, this will compute and store the postprocessed quantity.

The limits of the x and y axes are modified with, respectively
```
ae.xlim([0.1,0.3])
```
and
```
ae.ylim([1e6, 5e8])
```

It is possible to plot also the time evolution of a quantity for each radius, obtaining 2D colormaps. In order to obtain this kind of plots, you just have to specify it.
```
ae.entropy(projection='2D', plane='time')
```
It is also possible to overplot lines, by simply using the same command as before. For example, for the average shock radius use `ae.shock_radius()`.

Say, that, for example, we want to plot the average shock radius with a while solid line, and the maximum shock radius with a yellow dashed line. You will write
```
ae.shock_radius('avg', color='white', ls='-')
ae.shock_radius('max', color='yellow', ls='--')
```
You can also show the legend for these two new lines, giving the labels that you want
```
ae.update_legend(['average', 'maximum'])
```
Please note that the length of the list must be equal to the number of lines.

It is possible to plot the radial profile of the different quantities. For example, for the density
```
ae.rho(1.6, plane='radius')
```
This would plot the radial porofile of the density at 1.6 seconds after bounce. In multi-D (2D or 3D), this same syntax would plot the average over the angles. If you want the radial profile
at a given angular coordinate, you should pass a touple in the following way:
```
ae.rho(1.6, plane=(None, theta_idx, phi_idx))
```
where `theta_idx` and `phi_idx` are the indices for theta and phi at which you want the profile. If you want the angular profile, just keep free (i.e. `None`) the coordinate over which
you want to plot.

#### 2D slicing
Moreover, for a given timestep, you can plot a 2D slice for a given plane. Using the density as example, the code is
```
ae.rho(0.6, projection='2D', plane='xz')
```
This will plot a 2D slice of the xz-plane at 0.6 seconds after bounce. You can modify the levels of the colorbar with
```
ae.cbar_levels([1e6, 1e13])
```
This will restrict the values of the colorbar between $10^6$ and $10^{13}$.

With this last plot active, you can add a second one. For example, let's add the 2D slice of the entropy.
```
ae.entropy(0.6, projection='2D', plane='xz')
```
AeViz is able to store the previous data for the plot of the density (with all the changes) and it will create a new splitted plot with
density on one side (or, in general, the quantity that was already plotted) and entropy (the new quantity) on the other. All the information
on the density plot are kept.

All the subplots are identified by an uppercase letter, `'A'`, '`B`', '`C`', and so on. So if now you want to change the settings of the colorbar
for the entropy, the line of code will be
```
ae.cbar_levels([12.5, 20], 'B')
```
or, equivalently
```
ae.cbar_levels([12.5, 20], axd_letter='B')
```
where `'B'` refer to the second plot (the entropy one, in this case).

**Note**: The colorbar associated to each plot is identified with the corresponding lowercase letter. So `'a'` will refer to the colorbar of plot `'A'`, `'b'` to plot `'B'`, and so on.

You can zoom in and out with the `ae.xlim()` method introduced before.

**Note**: Usually, the horizontal axis will entend in the domain [0, L], while the verical one is [-L,L].

If you are interested to overplot given lines, such as the shock radius, you can do it with
```
ae.shock_radius('full', plot='B')
```
This will plot the shock radius contour (ray-by-ray) on plot `'B'`.

It is possible to also plot fields. Useful ones can be the magnetic filed
```
ae.add_field(comp='Bfield', plot='A')
```
or the (radial) velocity field
```
ae.add_field(comp='velocity', plot='B')
```

It is possible to obtain Hammer projection by selectine the 2D projection and specifying the radius
```
ae.rho(0.6, projection='2D', plane=2e6)
```
This will produce the Hammer projection at 0.6 seconds after bounce at a radius of 20 km (2e6 cm).

#### Spectrogram
Every quantity with a `'time'` plane can be shown with the associated spectrogram. This only needs one line of code.
```
ae.rho(projection='2D', plane='time', spectrogram=True)
```
This will generate two subplots, the top one with a simple line showing the time evolution, and the bottom one with the spectrogram.

This is of particular interest in the case of gravitational waves. To obtain such a figure use
```
ae.GWs(projection='2D', spectrogram=True)
```
