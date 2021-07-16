from db import *
from spec import *
from RV import *
from tools import *

# NEW FUNCTIONALITIES
table = findtable('IACOB_O9BAs_SNR20.fits') # file where quality flags are
table.add_index('ID')
table.loc['HD2905']

table_grouped = table.group_by('SpC_ref')
for group in table_grouped.groups:
    print(group,'\n')

table_grouped.groups[1].show_in_browser(jsviewer=True)

table['RAdeg_J2000'].format = ".8f"
table['DECdeg_J2000'].format = ".8f"
table['SpT_code'].format = ".3f"
# ---------------


#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# MAIN TABLES:
T_DB = findtable('IACOB_O9BAs_SNR20.fits') # file where quality flags are

T_IB = findtable('IB_results_ver5.txt') # file where vsini and vmac are
T_IB.remove_columns(['filename','line','snr'])

T_MAUI = findtable('MAUI_results.fits')
T_MAUI = T_MAUI[[str(i['Teff'])!='nan' and str(i['lgf'])!='nan' for i in T_MAUI]]
#T_deg = T_MAUI[[i in ['d','<','>'] or j in ['d','<','>'] for i,j in T_MAUI['l_Teff','l_lgf']]]
#T_deg.write(maindir + 'tables/degenerated.fits', format='fits', overwrite=True)
#from maui import *; #maui_input('degenerated.fits')
T_MAUI = T_MAUI[[i not in ['d','<','>'] and j not in ['d','<','>'] for i,j in T_MAUI['l_Teff','l_lgf']]]

T = join(T_DB, T_IB, keys='ID')
T = join(T, T_MAUI, keys='ID')

# OTHER TABLES:
T_DB_GH = findtable('IACOB_Os_Gon.fits')
T_Gbat_GH = findtable('Gbat_Gon.fits')
T_Gbat_GH = T_Gbat_GH[[i['eTeff'] != 0.0 and i['elgf'] != 0.0 for i in T_Gbat_GH]]

TGon_raw = join(T_DB_GH, T_Gbat_GH)

TGon = setdiff(TGon_raw, T, keys='ID') # Mine over Gon
TAll = join(TGon, T,  join_type='outer') # Mine over Gon



#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#T_i = T_DB[:10]
#x = zp_edr3(ra=T_i['RAdeg_J2000'], dec=T_i['DECdeg_J2000'], radius=0.5, IDs=T_i['ID'])

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
