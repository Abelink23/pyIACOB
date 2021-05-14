import numpy as np
import time
import scipy.constants as cte
from astropy.io import fits
from astropy.time import Time

import warnings
warnings.filterwarnings("ignore")

l_wl = None
r_wl = None
helcorr = 'hel'
width = 0.0
offset = 0.0

spectrum = "/Users/abelink/Desktop/HD129557_20130204_054417_F_V48000.fits"
#spectrum = "/Users/abelink/Desktop/HD2729_20131215_194517_M_V85000.fits"
#spectrum = "/Users/abelink/Desktop/HD2905_20081105_210738_N_V46000.fits"

hdu = fits.open(spectrum)# Open the fits image file
hdu0 = hdu[0]            # Load the header list of primary header
header0 = hdu0.header    # Read the values of the headers

x0 = header0['CRVAL1']          # Get the wavelenght of the first pixel
dx = header0['CDELT1']          # Step of increase in wavelength
pix0 = header0['CRPIX1']        # Reference pixel (generally 1, FEROS -49)
spec_length = header0['NAXIS1'] # Length of the spectrum

try: vbar = header0['I-VBAR'] # [km/s] Barycent. rv correction at midpoint
except: print('No helio/bary-centric correction is applied for ' + spectrum); vbar = 0

flux = []; wave = []
for j in range(0,spec_length):

    if l_wl != None and r_wl != None:
        if '_log' in spectrum:
            if np.exp(x0 + j*dx ) <= l_wl - width/2. or np.exp(x0 + j*dx) >= r_wl + width/2.: continue
        else:
            if float(x0 + j*dx) <= l_wl - width/2. or float(x0 + j*dx) >= r_wl + width/2.: continue

    if helcorr == 'hel' and not '_log' in spectrum:
        wave.append((x0 + j*dx)*(1 + 1000*vbar/cte.c) - offset)
    else:
        if '_log' in spectrum: wave.append(np.exp(x0 + j*dx) - offset)
        # Those with log are already corrected from helcorr
        else: wave.append(x0 + j*dx - offset)

    try: flux.append(hdu[0].data[0][j])
    except: flux.append(hdu[0].data[j])

wave = np.asarray(wave); flux = np.asarray(flux)


flux = []; wave = []
wave = x0 + dx*(np.arange(spec_length) - pix0 + 1)
if '_log' in spectrum: wave = np.exp(wave) - offset
elif helcorr == 'hel': wave = wave*(1 + 1000*vbar/cte.c) - offset
# Those with log are already corrected from helcorr

try: flux = hdu[0].data[0]
except: flux = hdu[0].data

if l_wl != None and r_wl != None:
    wave[(wave >= l_wl - width/2.) & (wave <= r_wl + width/2.)]
    flux = flux[(wave >= l_wl - width/2.) & (wave <= r_wl + width/2.)]

len(flux)

hdu.close()
