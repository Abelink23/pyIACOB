from spec import *
from RV import *

from tools import atokms
from scipy.signal import correlate,correlation_lags

RV_cc('HD37128',windows=[(3950,4160),(4310,4360),(4370,4490),(4540,4690),(4840,4950)])
#RV1_cc(findstar('HD37128',SNR='best'),synthetic)
#RV('B0.lst','HD37128',linesRV0='rv_Bs.lst')

plt.close('all')

spec1 = spec(findstar('HD37128',SNR='best'))
spec1.waveflux(lwl=3900,rwl=5080)

synthetic = []
for file in os.listdir(datadir+'ASCII/Synthetic_MAUI/'):
    if spec1.name_star in file: synthetic.append(file)

if len(synthetic) == 0:
    print('No files found for %s.\n' % (spec1.name_star))
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

#print(spec1.dx,spec2.dx); print(len(spec1.wave),len(spec2.wave))

spec2.resamp(dx=spec1.dx,lwl=spec1.wave[0],rwl=spec1.wave[-1])
spec1.resamp(dx=spec1.dx)

plt.plot(spec1.wave,spec1.flux,'b',lw=.5)
plt.plot(spec2.wave,spec2.flux,'g',lw=.5)

#print(spec1.dx,spec2.dx); print(len(spec1.wave),len(spec2.wave))

plt.figure()

corr = correlate(spec2.flux-1,spec1.flux-1)
corr /= np.max(corr)
lags = correlation_lags(len(spec1.flux),len(spec2.flux))
plt.plot(lags,corr)

plt.show(block=False)

print(atokms(-lags[np.argmax(corr)]*spec1.dx,spec1.wave.mean()))
