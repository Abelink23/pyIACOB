# Introduction to pyIACOB (v1.20 - 2026/05/05):

This package has been created mainly to manipulate data from the IACOB spectroscopic
database (see Simón-Díaz et al. 2011), which currently gathers high-resolution spectra
(R~25000-85000) of Galactic OBA stars taken with the FIES@NOT2.5m and HERMES@MERCATOR1.2m
facilities at the Roque de los Muchachos Observatory, and FEROS@MPG2.2m at La Silla
Observatory.

The main scientific goal pursued by the IACOB project is to use these data to perform
single snap-shot and multi-epoch quantitative spectroscopic analyses providing a complete
empirical overview of the physical properties of these stars.

To achieve this goal, several programs (written in IDL) have been developed by the team
of collaborators (mainly at the Instituto de Astrofísica de Canarias). This package has
modules that facilitate the use these programs either by preparing their input data or by
creating summary tables with their results. The current list of supported programs are:

- IACOB-Broad (S. Simón-Díaz & A. Herrero 2013)
- MAUI (M.A. Urbaneja 2026)

Created by: Abel de Burgos

# Python Requirements (installed automatically - see below):

- Python 3.12.2
- numpy 2.4.4
- scipy 1.17.1
- astropy 7.2.0
- astroquery 0.4.11
- matplotlib 3.10.8
- progressbar 2.5

# Other Requirements:

Access to the spectroscopic data stored locally. Available spectra from the IACOB
database can be downloaded from http://research.iac.es/proyecto/iacob/iacobcat/

# Installation:

Run the following command in the terminal in python 3.12.2 or above:

      python setup.py

The required folders will be created and the package will be installed.

NOTE: For the actual models needed in models.py to work, drop me an email and
      I will send you the needed files (~700mb).

      If you plan to use spectra in ascii format, you will need to put this files
      inside *data*/ASCII/

# Overview of the main modules:

- setup.py contains the code to install the package and create the required folders.

- db.py contains functions to search spectra or create master tables within the IACOB
  database, to search tables and lists and make queries in Simbad.

- reduction.py contains functions to perform the reduction of the spectra, including
  the data extraction, normalization, and formatting to match the IACOB database format.

- spec.py contains functions that manipulate an input spectrum (in FITS or ASCII format),
  plot it, make different measurements, and eventually export the new spectrum into a file.

- rv.py contains functions to calculate the radial velocity offset via cross-correlation
  or via input list of lines. It also allows to create RV curves.

- measure.py contains functions to interactively or automatically obtain Equivalent Widths,
  Full Widths, and Radial Velocities for a given input of stars and spectral lines.

- spec_posproc.py allows the user to perform a post-processing of the input spectra,
  including the radial velocity correction or the cosmetic defects and cosmic rays removal.

- models.py contains functions to retrieve either evolutionary tracks or isochrones from
  different libraries (MIST, Bonn, Geneva...).

- IACOBroad.py contains functions which helps the user to generate output and input
  tables for IACOB-Broad program.

- maui.py contains functions which helps the user to generate output and input tables for
  MAUI, with some additional features regarding the analysis of the results.

- binarity.py contains functions review multi-epoch spectra and perform a preliminary
  analysis of the binary nature of the target.

- get_feros.py contains functions to download FEROS spectra from the ESO archive.

- tools.py contains utility functions which could also be interested while using
  spectroscopic data.

# Sending corrections / comments:

- Preferably via pull requests / issues within the GitHub framework.
- Directly by email at adeburgo@eso.org

# Acknowledgments:

Sergio Simón-Díaz, Gonzalo Holgado, Alba Casabuenas, and Carlos Martinez-Sebastian

