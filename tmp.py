from spec import *
from RV import *
from tools import *

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#table_selection = findtable('table_selection.fits')
#for i in table_selection:
#    plt.figure()
#    spec(i['Name'],SNR='bestMF').plotspec(4000,5000)
#    plt.show(block='False')

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#spectra = findstar('HD14543')
#for i,j in zip(spectra,range(len(spectra))):
#    x = spec(i)
#    x.waveflux(3900,6800)
#    x.degrade(20000)
#    x.resamp(5*x.dx)
#    print('\n',x.file_name,x.snr)
#    x.flux = x.flux+j*0.05
#    x.plotspec()
#    x.plotline(4552.622,width=8)
#    plt.legend().remove()
#
#plt.show(block=False)

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#[spec(i).plotline(4552.622) for i in findstar('HD46149')]

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# UPDATE FITS HEADER
#spectrum = findstar('sddsds.fits')[0]
#hdu = fits.open(spectrum); header = hdu[0].header
#header['OBJECT'] = 'XX'
#header['I-SPC   '] = 'XXX'
#header['I-SPCREF']= 'SIMBAD'
#hdu.writeto(spectrum,output_verify='ignore',overwrite=True)
