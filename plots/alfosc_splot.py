import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

file = sys.argv[1]
try: width = int(sys.argv[2])
except:
    print('No width selected, using +-2 as default...')
    width = 2

hdu = fits.open('/home/abelink/Desktop/%s' % file)  # Open the fits image file
header1 = hdu[1].header         # Read header of primary extension

x0 = header1['CRVAL2']          # Get the wavelenght of the first pixel
dx = header1['CDELT2']          # Step of increase in wavelength
pix0 = header1['CRPIX2']        # Reference pixel (generally 1, FEROS -49)
spec_length = header1['NAXIS2'] # Length of the spectrum

wave = x0 + dx*(np.arange(spec_length) - pix0 + 1)

rows = hdu[1].data
max_flux_col = [row.tolist().index(max(row)) for row in rows]
max_flux_col = [val for val in max_flux_col if not val == 0]
spec_col = int(round(np.median(max_flux_col),0))

flux = [sum(row[spec_col-width:spec_col+width+1]) for row in rows]

plt.figure(tight_layout=True,figsize=(10,4))
plt.plot(wave[1:],flux[1:],lw=.5)
plt.gca().invert_xaxis()
plt.show()
