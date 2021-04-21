from RV import *
import random

def gen_ascii_ML(input_table='OBAs_ML_new.fits',not_do=None):

    table = findtable(input_table)

    for row in table:

        if not_do is not None and row['Name'] in not_do: continue
        if row['SNR_best'] > 90: continue

        next = 'n'
        while next == 'n':

            best_star = spec(row['Name'],SNR='bestMF')

            skip = input("%s - Hit return to continue, type 's' to skip: " % row['Name'])
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
            goodspec = findstar(row['Name'],SNR=30) #90/30
            try: goodspec = random.sample(goodspec,4) #25/4
            except: pass

            for j,k in zip(goodspec,range(len(goodspec))):
                star = spec(j,SNR='best')

                star.offset = RV0(spt_list,star.spectrum,ewcut=50,width=wid,tol=150,func=fun)#,plot='y')
                star.waveflux(3800,6900) # Applies the offset
                star.degrade(resol=5000)
                star.resamp(10*0.02564975,3800,6900)

                star.export(tail='_RV_ML_%i' % k,extension='.ascii')
