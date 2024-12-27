import sys
import os, subprocess, shutil
import platform

from distutils.command.clean import clean as _clean

from setuptools import Extension, setup, Command
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.install import install
from setuptools import find_packages

requires = []

def find_files(dirname, relpath=None):
    def find_paths(dirname):
        items = []
        for fname in os.listdir(dirname):
            path = os.path.join(dirname, fname)
            if os.path.isdir(path):
                items += find_paths(path)
            elif not path.endswith(".py") and not path.endswith(".pyc"):
                items.append(path)
        return items
    items = find_paths(dirname)
    if relpath is None:
        relpath = dirname
    return [os.path.relpath(path, relpath) for path in items]

def get_version():
    """Get the version info from the __version_ file"""

    with open(os.path.join('__version__')) as version_file:
        for line in version_file:
            if line.startswith('__version__'):
                return eval(line.split('=')[-1])

def get_requirements():
    """Get the requirements from the requirements file"""

    with open(os.path.join('requirements.txt')) as requirements_file:
        requirements = requirements_file.readlines()
    return [r.replace('\n', '') for r in requirements]

class CustomInstallCommand(install):
    def run(self):
        install.run(self)
        

VERSION = get_version()
setup_requires = []
install_requires = setup_requires + get_requirements()


setup(
    name = 'AeViz',
    version = VERSION,
    description = 'Library to analyze data from simulation of core ' \
        'collapse supernovae done with Aenus-ALCAR.',
    long_description = open('README.md').read(),
    long_description_content_type='text/markdown',
    author = 'Marco Cusinato',
    author_email = 'marco.cusinato@uv.es',
    url = 'https://github.com/MarcoCusinato/AeViz',
    download_url = f'https://github.com/MarcoCusinato/AeViz',
    keywords = [
        'Aenus-ALCAR',
        'astrophysics',
        'supernovae'
    ],
    setup_requires = setup_requires,
    install_requires = install_requires,
    scripts  = find_files('bin', relpath='./'),
    packages = find_packages(),
    package_data = {
        'AeViz.cell': find_files('AeViz/cell'),
        'AeViz.cell.cell_methods': find_files('AeViz/cell/cell_methods'),
        'AeViz.grid': find_files('AeViz/grid'),
        'AeViz.load_utils': find_files('AeViz/load_utils'),
        'AeViz.plot_utils': find_files('AeViz/plot_utils'),
        'AeViz.quantities_plotting': find_files('AeViz/quantities_plotting'),
        'AeViz.simulation': find_files('AeViz/simulation'),
        'AeViz.simulation.methods': find_files('AeViz/simulation/methods'),
        'AeViz.units': find_files('AeViz/units'),
        'AeViz.utils': find_files('AeViz/utils'),
        'AeViz.spherical_harmonics': find_files('AeViz/spherical_harmonics'),
        'AeViz.AeVizMethods': find_files('AeViz/AeVizMethods'),
    },
    python_requires='>=3.6',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Astrophysics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)