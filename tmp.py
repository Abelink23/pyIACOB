from RV import *
from tools import *

table = findtable('IACOB_O9BAs_SNR20.fits') # file where quality flags are
table_REF = findtable('RVEWFWs_O9BAs.fits') # file where RVs, EWs and FWs are
table_REF.remove_columns(['mySpC','SpC','SpT_code','LC_code'])
table_IB = findtable('IB_results_ver5.txt') # file where vsini and vmac are
table_IB.remove_columns(['filename','line','snr'])
results = findtable('MAUI_results.fits')
gonzalo_raw = findtable('Gon_results.fits')

table_f = join(table,table_REF,keys='Name')
table_f = join(table_f,table_IB,keys='Name')
table_f = join(table_f,results,keys='Name')

gonzalo = setdiff(gonzalo_raw,table_f,keys='Name')

table_f = table_f[[i in ['d','<','>'] or j in ['d','<','>'] for i,j in table_f['l_Teff','l_lgf']]]

table_f.write(maindir+'tables/degenerated.fits',format='fits',overwrite=True)

from maui import *
maui_input('degenerated.fits')

#x = zp_edr3(table='Gaia.fits',ra='ra_epoch2000',dec='dec_epoch2000',search_radius=.5)


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
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

#[spec(i).plotline(4552.622) for i in findstar('HD195556',SNR='best')]

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#for row in table:
#    if row['FIES'] > 0 and '_N_' in findstar(row['Name'],'best')[0] and (row['mag_B']+row['mag_V'])/2 < 8.5:
#        if row['HERMES']+row['FEROS'] == 0: print('Check',row['Name']); continue
#        spectra = findstar(row['Name']); best_SNR = row['SNR_best']
#        check = True
#        for spectrum in spectra:
#            if spec(spectrum).snr > 100 and '_N_' not in spec(spectrum).file_name:
#                check = False
#        if check == True: print('Check',row['Name'])
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# UPDATE FITS HEADER
#spectrum = findstar('sddsds.fits')[0]
#hdu = fits.open(spectrum); header0 = hdu[0].header
#header0['OBJECT'] = 'XX'
#header0['I-SPC   '] = 'XXX'
#header0['I-SPCREF']= 'SIMBAD'
#hdu.writeto(spectrum,output_verify='ignore',overwrite=True)
