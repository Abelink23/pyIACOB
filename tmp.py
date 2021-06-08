from spec import *
from RV import *
from tools import *

x=findtable('ALL_OBs_n4+.txt',delimiter=',')

x['ion'] = [i.split(' ')[1] for i in x['spc']]
y = x[(x['elem'] == 'Si') & (x['ion'] =='III') & (x['configuration'] == '3s.4p-3s.4d') ]
z = setdiff(x,y,keys='wl_air')
y.sort('term')
terms = list(set(y['term']))
configs =list(set(y['configuration']))

he = x[['He' in i for i in x['spc']]].sort('wl_air')
he[(he['wl_air'] > 4500) & (he['wl_air'] < 5000) & (he['-lg(gf)'] > -2)]

x=findtable('ALL_OBs_n4+.txt',delimiter=',')
x[(x['wl_air'] > 4310) & (x['wl_air'] < 4370)]# & (x['-lg(gf)'] > -2)]


#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#table_selection = findtable('table_selection.fits')
#for i in table_selection:
#    plt.figure()
#    spec(i['ID'],SNR='bestMF').plotspec(4000,5000)
#    plt.show(block='False')

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#spectra = findstar('HD14543')
#for i,j in zip(spectra,range(len(spectra))):
#    x = spec(i)
#    x.waveflux(3900,6800)
#    x.degrade(20000)
#    x.resamp(5*x.dlam)
#    print('\n',x.filename,x.snr)
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
