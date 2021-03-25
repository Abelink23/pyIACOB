import os

from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord

with fits.open('RosettaSV/GALANTE_NGC_2244.fits') as hdu_list:
    hdu_list = Table.read(hdu_list)

ra = hdu_list['RA']
dec = hdu_list['DEC']

radec = SkyCoord(hdu_list['RA'],hdu_list['DEC'],unit='deg')

hdu_list.add_column(radec.galactic.l.deg,name='gal_long',index=5)
hdu_list.add_column(radec.galactic.l.deg,name='gal_lat',index=6)

hdu_list.write('RosettaSV/GALANTE_NGC_2244_phot.fits',overwrite=True)
