from spec import *
from RV import *

#table_selection = findtable('table_selection.fits')
#for i in table_selection:
#    plt.figure()
#    spec(i['Name'],SNR='bestMF').plotspec(4000,5000)
#    plt.show(block='False')

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#spectra = findstar('HD34447')
#for i,j in zip(spectra,range(len(spectra))):
#    x = spec(i)
#    x.waveflux(3900,6800)
#    x.degrade(20000)
#    x.resamp(5*x.dx)
#    print('\n',x.file_name,x.snr)
#    x.flux = x.flux+j*0.05
#    #x.plotspec(4500,4600)
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
#spectrum = findstar('HD193237')[0]
#hdu = fits.open(spectrum)  # Open the fits image file
#header0 = hdu[0].header         # Read header of primary extension
#header0['I-SNR'] = 61   # Instrument
#hdu.writeto(spectrum,output_verify='ignore',overwrite=True)
