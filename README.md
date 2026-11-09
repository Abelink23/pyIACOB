# Introduction to pyIACOB (v1.24 - 2026/11/10):
This package has been created mainly to reduce, review, modify, and analyze spectra from the IACOB spectroscopic database (see Simón-Díaz et al. 2020) within a python environment. The IACOB database currently gathers high-resolution spectra (R~25000-85000) of Galactic OBA stars. These spectra are taken with the FIES@NOT2.5m and HERMES@MERCATOR1.2m instruments at the Roque de los Muchachos Observatory, and the FEROS@MPG2.2m instrument at La Silla Observatory.

The main scientific goal pursued by the IACOB project is to perform single snap-shot and multi-epoch quantitative spectroscopic analyses providing a complete empirical overview of the physical properties of these stars.

Additionally, this package facilitates the use of MAUI (M.A. Urbaneja 2026), a tool that allows the user to perform a quantitative spectroscopic analysis of OB-type stars using grids of FASTWIND (Puls et al. 2005) stellar atmosphere models.

This project is developed by Abel de Burgos Sierra since 2019, and has continuesly evolved with different contributions from Sergio Simón-Díaz, Alba Casabuenas, Daniel Galan, Gonzalo Holgado, and Carlos Martinez-Sebastian, who are all members of the IACOB project. The package is currently maintained by Abel de Burgos Sierra.

# Installation
Clone the repository from GitHub:

      git clone https://github.com/Abelink23/pyIACOB/

Or download the zip file from GitHub and extract it in your preferred location.

## Python version:
pyIACOB requires **Python 3.11 or higher.**. It is recommended to install a python environment via Conda or venv to avoid conflicts with other packages:

For Conda (or Miniconda) users:

      conda create -n pyIACOB python=3.11
      conda activate pyIACOB

For venv users, go to the main folder of the package and run:

      python -m venv pyIACOB
      source pyIACOB/bin/activate  (Linux/MacOS)
      pyIACOB\Scripts\activate     (Windows)

## Installation script:
Run the following command in the terminal from the main folder of the package:

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

## Downloading the IACOB spectra:
Spectra from the IACOB spectroscopic database stored locally. IACOB spectra can be downloaded from: https://research.iac.es/proyecto/iacob/iacobcat/

If you plan to use spectra in ascii format, you will need to put this files inside *data*/ASCII/


# Overview of the main modules:
- binarity.py contains functions review multi-epoch spectra and perform a preliminary analysis of the binary nature of the target.

- color_msg.py contains functions to print colored messages in the terminal.

- db.py is the main module of the package, which loads the default paths, and provides functions to search for lists, tables, files, spectra (by, e.g. S/N); generate master tables; make queries in Simbad or Gaia; check the integrity of the fits files; and other functions to interact with the IACOB database.

- get_feros.py contains functions to download FEROS datasets from the ESO archive and only keep the files that contain the reduced spectra of the target star.

- hardcoded_maui.py contains hardcoded values, parameters, or masks used in the MAUI module.

- iacob_broad.py contains the function to perform the line-broadening analysis of OB-type stars using the IACOB-BROAD tool. This module is a Python implementation of the original IDL tool developed by Simón-Díaz, S. & Herrero, A. (2014), and later translated to Python by D. Galán (PhD, IAC) together with the sub-modules included in the broadening/ folder.

- isocrones.py contains functions to retrieve isochrones from different libraries [IN DEVELOPMENT].

- maui.py contains functions which helps the user to generate output and input tables for MAUI (Urbaneja 2026), with some additional features regarding the preparation of the analysis and the review or visualization of the results.

- measure.py contains functions to interactively or automatically obtain equivalent widths, full widths, full widths at 3/4 and 1/4 of the line depth, or radial velocities or a given input of stars and spectral lines.

- reduction.py contains functions to perform the reduction of the spectra, including the data extraction, normalization, and formatting to match the IACOB database format.

- rv.py contains functions to calculate the radial velocity offset via cross-correlation or via input list of lines. It also allows to create RV curves.

- setup.py contains the code to install the package, the dependencies, and create the required folders.

- spec.py contains functions that manipulate an input spectrum (in FITS or ASCII format), plot it, make different measurements, and eventually export the new spectrum into a file. See also fitline_diagram.txt and fitline_readme.txt for a description of the line fitting procedure.

- spec_posproc.py allows the user to perform a post-processing of the input spectra, including the radial velocity correction or the cosmetic defects and cosmic rays removal.

- tools.py contains utility functions which could also be interested while using spectroscopic data.

- tracks.py contains functions to use evolutionary tracks from different sources inclued in the tracks/ folder. Currently it includes the tracks from Brott et al. (2011) and Ekström et al. (2012), Sciarini et al. (2026), and the MIST tracks.

# Usage:

To use the package, simply import everything from the pyIACOB module you need to use. For example, to use the functions in `spec.py`, you would do:

      from spec import *
      s = spec('path/to/spectrum.fits')
      s.plotspec()


# Sending corrections / comments:

- Preferably via pull requests / issues within the GitHub framework.
- Directly by email at adeburgo@eso.org
