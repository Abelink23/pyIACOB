# gaiadr3_zeropoint

This Python package contains the necessary tools to query the value of the parallax zero-point for Gaia EDR3  and
 Gaia DR3. Based on the functions described in Lindegren et al. 2020, the code returns the estimated parallax zero
 -point given the ecliptic latitude, magnitude and colour of any Gaia (E)DR3 source.

## Documentation

All classes and methods/functions are documented so use the python help() function to find out more. 


## Installation

This is a Python3 package (*issues may arise if executed with Python2*).

The required dependencies are:
* [numpy](https://numpy.org/)
* [pandas](https://pandas.pydata.org/) (only if you want to use the wrapper provided with the code)


To install the package:

### From source (recommended)
1. Clone the Github repository or download the source files
2. cd to the directory
3. Run `python setup.py install` or `python setup install --user` for installation in your own home directory

### With pip
```
pip install gaiadr3-zeropoint
```


## Basic usage

Once the package is installed, you can import it in Python:

```
from zero_point import zpt
```

Then, first load the coefficient tables by calling the `load_tables()` function. 

Once the tables are loaded, the
 parallax zero-point can be queried as:

```
zpt.get_zpt(phot_g_mean_mag, nu_eff_used_in_astrometry, pseudocolour, ecl_lat, astrometric_params_solved)
```

This function accepts both single values as well as iterables, and returns a float (or array of such) corresponding to the zero-point of the source(s) with those parameters.

**NOTE**: for 5-p solutions (ra-dec-parallax-pmra-pmdec), the field `astrometric_params_solved` equals 31 and the
 `pseudocolour` variable can take any arbitrary values (even *None*). On the other hand, for 6-p solutions (ra-dec
 -parallax-pmra-pmdec-pseudocolour), the field `astrometric_params_solved` equals 95 and the
  `nu_eff_used_in_astrometry` variable can take any arbitrary values (even *None*).

Finally, if you have a pandas DataFrame (DF) of sources with the columns `phot_g_mean_mag, nu_eff_used_in_astrometry, pseudocolour, ecl_lat, astrometric_params_solved`, you can simply use the pandas wrapper ```zpt_wrapper```:

``` 
zero_point = DF.apply(zpt_wrapper,axis=1) 
```


## Attribution

If you make use of this package for your research, please acknowledge the following papers: Lindegren+20.

## Help

If you encounter any problem with the software, please make use of the GitLub Issues page. Otherwise, contact p.ramos@unistra.fr.

Copyright: Pau Ramos, University of Barcelona

