# Introduction to pyIACOB (v1.24 - 2026/11/10):

This package has been created mainly to reduce, review, modify, and analyze spectra from the IACOB spectroscopic database (see Simón-Díaz et al. 2020) within a python environment. The IACOB database currently gathers high-resolution spectra (R~25000-85000) of Galactic OBA stars. These spectra are taken with the FIES@NOT2.5m and HERMES@MERCATOR1.2m instruments at the Roque de los Muchachos Observatory, and the FEROS@MPG2.2m instrument at La Silla Observatory.

The main scientific goal pursued by the IACOB project is to perform single snap-shot and multi-epoch quantitative spectroscopic analyses providing a complete empirical overview of the physical properties of these stars.

Additionally, this package facilitates the use of MAUI (M.A. Urbaneja 2026), a tool that allows the user to perform a quantitative spectroscopic analysis of OB-type stars using grids of FASTWIND (Puls et al. 2005) stellar atmosphere models.

This project is developed by Abel de Burgos Sierra since 2019, and has continuesly evolved with different contributions from Sergio Simón-Díaz, Alba Casabuenas, Daniel Galan, Gonzalo Holgado, and Carlos Martinez-Sebastian, who are all members of the IACOB project. The package is currently maintained by Abel de Burgos Sierra.

# Python version required:
- Python 3.11 or higher.

It is recommended to install it via Conda or venv to avoid conflicts with other packages:

For Conda (or Miniconda) users:

      conda create -n pyIACOB python=3.11
      conda activate pyIACOB

For venv users, go to the main folder of the package and run:

      python -m venv pyIACOB
      source pyIACOB/bin/activate  (Linux/MacOS)
      pyIACOB\Scripts\activate     (Windows)


# Python dependences:
Note: checked and installed automatically when running the setup.py script (see below).
- astropy 8.0.2
- astroquery 0.4.11
- matplotlib 3.11.2
- numpy 2.4.6
- progressbar 4.6.0
- scipy 1.17.1

# Other Requirements:

Access to the spectroscopic data stored locally. Available spectra from the IACOB
database can be downloaded from https://research.iac.es/proyecto/iacob/iacobcat/

# Installation:

Run the following command in the terminal from the main folder of the package:

      python setup.py

The required folders will be created and the package will be installed.

      If you plan to use spectra in ascii format, you will need to put this files
      inside *data*/ASCII/

# Overview of the main modules:

- binarity.py contains functions review multi-epoch spectra and perform a preliminary
  analysis of the binary nature of the target.

- color_msg.py contains functions to print colored messages in the terminal.

- db.py contains functions to search spectra or create master tables within the IACOB
  database, to search tables and lists and make queries in Simbad.

- get_feros.py contains functions to download FEROS spectra from the ESO archive.

- hardcoded_maui.py contains hardcoded values, parameters, or masks used in the MAUI analysis.

- iacob_broad.py contains functions to perform the line-broadening analysis of OB-type
  stars using the IACOB-BROAD tool.

- isocrones.py contains functions to retrieve isochrones from different libraries (MIST, Bonn, Geneva...) [IN DEVELOPMENT].

- maui.py contains functions which helps the user to generate output and input tables for
  MAUI, with some additional features regarding the analysis of the results.

- measure.py contains functions to interactively or automatically obtain Equivalent Widths,
  Full Widths, and Radial Velocities for a given input of stars and spectral lines.

- reduction.py contains functions to perform the reduction of the spectra, including
  the data extraction, normalization, and formatting to match the IACOB database format.

- rv.py contains functions to calculate the radial velocity offset via cross-correlation
  or via input list of lines. It also allows to create RV curves.

- setup.py contains the code to install the package and create the required folders.

- spec.py contains functions that manipulate an input spectrum (in FITS or ASCII format),
  plot it, make different measurements, and eventually export the new spectrum into a file.

- spec_posproc.py allows the user to perform a post-processing of the input spectra,
  including the radial velocity correction or the cosmetic defects and cosmic rays removal.

- tools.py contains utility functions which could also be interested while using
  spectroscopic data.

- tracks.py contains functions to retrieve either evolutionary tracks or isochrones from
  different libraries (MIST, Bonn, Geneva...).

# Usage:

To use the package, simply import everything from the pyIACOB module you need to use. For example, to use the functions in `spec.py`, you would do:

      from spec import *
      s = spec('path/to/spectrum.fits')
      s.plotspec()


# Sending corrections / comments:

- Preferably via pull requests / issues within the GitHub framework.
- Directly by email at adeburgo@eso.org
