from spec import *

from tools import atokms
from scipy.signal import correlate,correlation_lags

plt.close('all')

spec1 = spec(findstar('HD7902')[-1])
spec1.waveflux(lwl=3900,rwl=5080)

spec2 = spec('HD7902_20121025_223409_M_V85000_red1.dat',txt=True)
spec2.txtwaveflux(lwl=3900,rwl=5080)
plt.plot(spec2.wave,spec2.flux,'r',lw=.5)

resol = 1/np.sqrt((1/spec1.resolution)**2-(1/spec2.resolution)**2)
spec2.degrade(resol=7000)

print(spec1.dx,spec2.dx); print(len(spec1.wave),len(spec2.wave))

spec2.resamp(dx=spec1.dx,lwl=spec1.wave[0],rwl=spec1.wave[-1])
spec1.resamp(dx=spec1.dx)

plt.plot(spec1.wave,spec1.flux,'b',lw=.5)
plt.plot(spec2.wave,spec2.flux,'g',lw=.5)

print(spec1.dx,spec2.dx); print(len(spec1.wave),len(spec2.wave))

plt.figure()
corr = correlate(spec2.flux-1,spec1.flux-1)
corr /= np.max(corr)
lags = correlation_lags(len(spec1.flux),len(spec2.flux))
plt.plot(lags,corr)

plt.show(block=False)

print(atokms(-lags[np.argmax(corr)]*spec1.dx,spec1.wave.mean()))
