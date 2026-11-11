# Introduction to pyIACOB (v1.24 - 2026/11/10):
This package has been created mainly to reduce, review, modify, and analyze spectra from the IACOB spectroscopic database (see Simón-Díaz et al. 2020) within a python environment. The IACOB database currently gathers high-resolution spectra (R~25000-85000) of Galactic OBA stars. These spectra are taken with the FIES@NOT2.5m and HERMES@MERCATOR1.2m instruments at the Roque de los Muchachos Observatory, and the FEROS@MPG2.2m instrument at La Silla Observatory.

The main scientific goal pursued by the IACOB project is to perform single snap-shot and multi-epoch quantitative spectroscopic analyses providing a complete empirical overview of the physical properties of these stars.

This project is developed by Abel de Burgos Sierra since 2019, and has continuesly evolved with different contributions from the IACOB team. This particular branch of the package (IACOB-broad) includes a Python implementation of the original IDL tool developed by Simón-Díaz, S. & Herrero, A. (2014), and later translated to Python by Daniel Galán (PhD, IAC) together with other sub-modules maintained by Dr. Abel de Burgos Sierra (ESO).

# Installation
Clone the repository from GitHub (preferable for keeping the package up to date):

      git clone https://github.com/Abelink23/pyIACOB/tree/iacob-broad

Or download the zip file from GitHub and extract it in your preferred location.

Note: if your iacob-broad branch is older than the version indicated above, update it by running this command from the main folder of the package:

      git pull origin iacob-broad

## Python version:
pyIACOB requires **Python 3.11 or higher**. It is recommended to install a python environment via Conda or venv to avoid conflicts with other packages:

For Conda (or Miniconda) users:

      conda create -n <name_of_the_env> python=3.11
      conda activate <name_of_the_env>

Alternatively, for venv users with python 3.11 installed, go to your directory of choice (e.g., working directory or directory with other environments, but not the main folder of the package), and run:

      python -m venv <name_of_the_env>
      source <name_of_the_env>/bin/activate  (Linux/MacOS)
      <name_of_the_env>\Scripts\activate     (Windows)

## Installation script:
Run the following command in the terminal from the main folder of the package (pyIACOB):

      python setup.py

The required folders will be created and the package and its dependencies will be installed (see below).

## Python dependences:
Note: checked and installed automatically when running the setup.py script (see above).
- astropy 8.0.2
- astroquery 0.4.11
- matplotlib 3.11.2
- numpy 2.4.6
- progressbar 4.6.0
- scipy 1.17.1

If you use ipython or jupyter notebooks, it is recommended to install the corresponding packages.

## Downloading the IACOB spectra:
Spectra from the IACOB spectroscopic database stored locally. IACOB spectra can be downloaded from: https://research.iac.es/proyecto/iacob/iacobcat/

If you plan to use spectra in ascii format, you will need to put this files inside *data*/ASCII/


# Overview of the main modules:
- binarity.py contains functions review multi-epoch spectra and perform a preliminary analysis of the binary nature of the target.

- color_msg.py contains functions to print colored messages in the terminal.

- db.py is the main module of the package, which loads the default paths, and provides functions to search for lists, tables, files, spectra (by, e.g. S/N); generate master tables; make queries in Simbad or Gaia; check the integrity of the fits files; and other functions to interact with the IACOB database.

- iacob_broad.py contains the function to perform the line-broadening analysis of OB-type stars. It calculates the projected rotational velocity (vsini) and macroturbulent broadening (vmac) of a spectral line with the Fourier Transform + Goodness-of-Fit method.

- setup.py contains the code to install the package, the dependencies, and create the required folders.

- spec.py contains functions that manipulate an input spectrum (in FITS or ASCII format), plot it, make different measurements, and eventually export the new spectrum into a file. See also fitline_diagram.txt and fitline_readme.txt for a description of the line fitting procedure.

# Usage:

To use the package, simply import everything from the pyIACOB module you need to use. For example, to use the functions in `spec.py`, you would do:

      from iacob_broad import *
      s = spec('<name_of_the_star/spectrum>')
      s.plotspec()
      results = iacob_broad(s, snr='best', line='SiIII')              # interactive
      results = iacob_broad(s, snr='best', line='SiIII', cut=True)    # pre-processing only
      results = iacob_broad(s, snr='best', line='SiIII', noproc=True) # skip pre-processing if already done and GOF is computed

Note: You can use ascii spectra not coming from the IACOB database (e.g. synthetic spectra), but you will need to provide the wavelength and flux arrays in the correct format (see spec.py for details).

      s = spec('path/to/syn/spectrum.ascii', orig='ascii', outside_dir='path/to/syn')


# Sending corrections / comments:

- Preferably via pull requests / issues within the GitHub framework.
- Directly by email at adeburgo@eso.org
