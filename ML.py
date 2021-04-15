from spec import *
from RV import *
import random

table = findtable('OBAs_ML_f.fits')

for i in table:
    best_star = spec(i['Name'],SNR='bestMF')

    if i['SpT_code'] < 2: list = 'rv_Os.lst'
    elif i['SpT_code'] >= 2 and i['SpT_code'] < 3: list = 'rv_Bs.lst'
    elif i['SpT_code'] >= 3: list = 'rv_As.lst'

    best_star.offset = RV0(list,best_star.spectrum,func='g',ewcut=50,tol=150)#,plot='y')
    best_star.waveflux(3800,6900) # Applies the offset
    best_star.degrade(resol=5000)
    best_star.resamp(10*0.02564975,3800,6900)
    best_star.plotspec(4500,4600,poslines='OB')

    all = input('Continue with the rest of available spectra? [y/n]: '); plt.close()
    if all == 'y' or all == '':
        goodspec = findstar(i['Name'],SNR=100)
        try: goodspec = random.sample(goodspec,20)
        except: pass

        for j,k in zip(goodspec,range(len(goodspec))):
            star = spec(j,SNR='best')
            star.offset = RV0(list,star.spectrum,func='g',ewcut=50,tol=150)#,plot='y')
            star.waveflux(3800,6900) # Applies the offset
            star.degrade(resol=5000)
            star.resamp(10*0.02564975,3800,6900)

            star.export(tail='_RV_ML_%i' % k,extension='.ascii')
