from RV import *
import random

def gen_ascii_ML(input_table='OBAs_ML_new.fits',not_do=None):

    table = findtable(input_table)

    for row in table:

        if not_do is not None and row['ID'] in not_do: continue
        if row['SNR_best'] > 90: continue

        next = 'n'
        while next == 'n':

            best_star = spec(row['ID'],SNR='bestMF')

            skip = input("%s - Hit return to continue, type 's' to skip: " % row['ID'])
            if skip == 's': break

            if row['SpT_code'] <= 2.5: best_star.plotspec(4500,4600)
            else: best_star.plotspec(6337.11,6357.11)

            if row['SpT_code'] < 2: spt_list = 'rv_Os.lst'
            elif row['SpT_code'] >= 2 and row['SpT_code'] < 3: spt_list = 'rv_Bs.lst'
            elif row['SpT_code'] >= 3: spt_list = 'rv_As.lst'
            fun = '-'
            while fun not in ['g','l','v','r','vr']:
                fun = input('Choose function to fit between g/l/v/r/vr (default is g): ')
                if fun == '': fun = 'g'
            wid = '-'
            while type(wid) is not float:
                wid = input('Choose the initial width in angstroms (default is 15): ')
                if wid == '': wid = 15.
                else: wid = float(wid)

            plt.close()

            best_star.offset = RV0(spt_list,best_star.spectrum,ewcut=50,width=wid,tol=150,func=fun)#,plot='y')
            best_star.waveflux(3800,6900) # Applies the offset
            best_star.degrade(resol=5000)
            best_star.resamp(10*0.02564975,3800,6900)

            if row['SpT_code'] <= 2.5: best_star.plotspec(4500,4600,poslines='OB',ylim=(0.75,1.05))
            else: best_star.plotspec(6337.11,6357.11,poslines='OB',ylim=(0.75,1.05))

            next = input('Hit return to continue with the rest of spectra or any other key to repeat. ')

            plt.close('all')

        if next == '':
            goodspec = findstar(row['ID'],SNR=30) #90/30
            try: goodspec = random.sample(goodspec,4) #25/4
            except: pass

            for j,k in zip(goodspec,range(len(goodspec))):
                star = spec(j,SNR='best')

                star.offset = RV0(spt_list,star.spectrum,ewcut=50,width=wid,tol=150,func=fun)#,plot='y')
                star.waveflux(3800,6900) # Applies the offset
                star.degrade(resol=5000)
                star.resamp(10*0.02564975,3800,6900)

                star.export(tail='_RV_ML_%i' % k,extension='.ascii')


#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# Update tables with new data | / IMPLEMENT IN DB AT SOME POINT!!
#table1,table2 = findtable('OBAs_ML_raw.fits'),findtable('MAUI_results.fits')
#columns_to_update = table2.colnames[1:-9]+[table2.colnames[-1]]
#for row,i in zip(table1,range(len(table1))):
#    updated = False
#    if row['ID'] not in table2['ID']: continue
#    if row['Teff'] == table2[table2['ID'] == row['ID']]['Teff'].data[0]: continue
#    for col_name in columns_to_update:
#        table1[i][col_name] = table2[table2['ID'] == row['ID']][col_name].data[0]
#        if updated == False: print(row['ID']); updated = True
#table1.write(maindir+'tables/table1_updated.fits',format='fits',overwrite=True)

# This one is to empty bad data
#table1 = findtable('OBAs_ML_raw.fits')
#columns_to_update = [i for i in table1.columns[36:-1]]
#for row,i in zip(table1,range(len(table1))):
#    for col_name in columns_to_update:
#        if row[col_name] in ['d','<','>']:
#            for j in columns_to_update[columns_to_update.index(col_name)+1:columns_to_update.index(col_name)+4]:
#                table1[i][j] = np.nan
#table1.write(maindir+'tables/OBAs_ML_ver1.fits',format='fits',overwrite=True)

#table1 = findtable('OBAs_ML_ver1b.fits')
#for row,i in zip(table1,range(len(table1))):
#    if row['QSiIII'] < 3:
#        table1[i]['EWSiIII1'] = table1[i]['FWSiIII1'] = table1[i]['depSiIII1'] = np.nan
#        table1[i]['EWSiIII2'] = table1[i]['FWSiIII2'] = table1[i]['depSiIII2'] = np.nan
#        table1[i]['EWSiIII3'] = table1[i]['FWSiIII3'] = table1[i]['depSiIII3'] = np.nan
#    if row['QSiII'] < 3:
#        table1[i]['EWSiII'] = table1[i]['FWSiII'] = table1[i]['depSiII'] = np.nan
#    if row['QHb'] < 3:
#        table1[i]['EWHb'] = table1[i]['FWHb'] = table1[i]['FW14Hb'] = table1[i]['FW34Hb'] = table1[i]['depHb'] = table1[i]['gamma'] = np.nan
#table1.write(maindir+'tables/OBAs_ML_ver1.fits',format='fits',overwrite=True)

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# Replace values in table
#table = findtable('OBAs_ML_raw.fits')
#for row,i in zip(table,range(len(table))):
#    for column,j in zip(row,range(len(row))):
#        #if column == 'N': table[i][j] = '='
#        if str(column) == str(1e+20): table[i][j] = np.nan
#table.write(maindir+'tables/OBAs_ML_raw_.fits',format='fits',overwrite=True)
