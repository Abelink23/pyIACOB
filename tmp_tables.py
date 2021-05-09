from spec import *
from RV import *
from tools import *


#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#table = findtable('IACOB_O9BAs_SNR20.fits') # file where quality flags are
#table_REF = findtable('RVEWFWs_O9BAs.fits') # file where RVs, EWs and FWs are
#table_REF.remove_columns(['mySpC','SpC','SpT_code','LC_code','Comments'])
#table_IB = findtable('IB_results_ver5.txt') # file where vsini and vmac are
#table_IB.remove_columns(['filename','line','snr'])
#results = findtable('MAUI_results.fits')
#gonzalo_raw = findtable('Gon_results.fits')
#
#table_f = join(table,table_REF,keys='ID')
#table_f = join(table_f,table_IB,keys='ID')
#table_f = join(table_f,results,keys='ID')
#
#table_f = table_f[[i in ['d','<','>'] or j in ['d','<','>'] for i,j in table_f['l_Teff','l_lgf']]]
#
#table_f.write(maindir+'tables/degenerated.fits',format='fits',overwrite=True)
#
#from maui import *
#maui_input('degenerated.fits')

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#x = zp_edr3(table='Gaia.fits',ra='ra_epoch2000',dec='dec_epoch2000',search_radius=.5)

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#for row in table:
#    if row['FIES'] > 0 and '_N_' in findstar(row['ID'],'best')[0] and (row['mag_B']+row['mag_V'])/2 < 8.5:
#        if row['HERMES']+row['FEROS'] == 0: print('Check',row['ID']); continue
#        spectra = findstar(row['ID']); best_SNR = row['SNR_best']
#        check = True
#        for spectrum in spectra:
#            if spec(spectrum).snr > 100 and '_N_' not in spec(spectrum).file_name:
#                check = False
#        if check == True: print('Check',row['ID'])
