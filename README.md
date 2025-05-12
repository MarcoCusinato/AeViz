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
