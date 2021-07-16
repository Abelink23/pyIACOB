import sys
sys.path.append('../')

from spec import *
from RV import *

from tools import atokms
from scipy.signal import correlate,correlation_lags

#RV_cc('HD37128',windows=[(3950,4160),(4310,4360),(4370,4490),(4540,4690),(4840,4950)])
#RV1_cc(findstar('HD37128',SNR='best'),synthetic)
#RV('B0.lst','HD37128',linesRV0='rv_Bs.lst')

plt.close('all')

spec1 = spec(findstar('HD37128',SNR='best'))
spec1.waveflux(lwl=3900,rwl=5080)

synthetic = []
for file in os.listdir(datadir+'ASCII/Synthetic_MAUI/'):
    if spec1.id_star in file: synthetic.append(file)

if len(synthetic) == 0:
    print('No files found for %s.\n' % (spec1.id_star))
elif len(synthetic) == 1:
    synthetic = synthetic[0]
else:
    for name,i in zip(synthetic,range(len(synthetic))):
        print(name,i)
    which = input('Enter the number of the synthetic spectra you want to use: ')
    synthetic = synthetic[int(which)]

spec2 = spec(synthetic,txt=True)
spec2.txtwaveflux(lwl=3900,rwl=5080)
plt.plot(spec2.wave,spec2.flux,'r',lw=.5)

resol = 1/np.sqrt((1/spec1.resolution)**2-(1/85000)**2)
if not resol == np.inf: spec2.degrade(resol=resol)

#print(spec1.dlam,spec2.dlam); print(len(spec1.wave),len(spec2.wave))

spec2.resamp(dlam=spec1.dlam,lwl=spec1.wave[0],rwl=spec1.wave[-1])
spec1.resamp(dlam=spec1.dlam)

mask = [spec2.flux < .999]
spec1.wave = spec1.wave[mask]; spec1.flux = spec1.flux[mask]
spec2.wave = spec2.wave[mask]; spec2.flux = spec2.flux[mask]

plt.plot(spec1.wave,spec1.flux,'b',lw=.5)
plt.plot(spec2.wave,spec2.flux,'g',lw=.5)

#print(spec1.dlam,spec2.dlam); print(len(spec1.wave),len(spec2.wave))

windows=[(3950,4160),(4310,4360),(4370,4490),(4540,4690),(4840,4950)]
#windows=[(3950,4950)]
RVs_angs = []; RVs_kms = []; plt.figure()
for win in windows:
    flux2 = spec2.flux[(spec2.wave >= win[0]) & (spec2.wave <= win[1])]
    flux1 = spec1.flux[(spec1.wave >= win[0]) & (spec1.wave <= win[1])]
    wave1 = spec1.wave[(spec1.wave >= win[0]) & (spec1.wave <= win[1])]

    corr = correlate(flux1-1,flux2-1)
    corr /= np.max(corr)

    lags = correlation_lags(len(flux1),len(flux2))
    corr_shift = lags[np.argmax(corr)]
    plt.plot(lags,corr)

    #lags = np.arange(-len(flux1)+1,len(flux1),1)
    #corr_shift = corr.argmax()-len(flux1)+1
    #plt.plot(lags,corr)

    RVs_angs.append(corr_shift*spec1.dlam)
    RVs_kms.append(RVs_angs[-1]/np.mean(wave1)*cte.c/1000)

RV_angs = round(np.mean(RVs_angs),8)
RV_kms = round(np.mean(RVs_kms),4)

plt.title(spec1.id_star+' '+str(RV_angs)+'A '+str(RV_kms)+'km/s ')
plt.show(block=False)
