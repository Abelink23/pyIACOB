import sys
sys.path.append('../')

from spec import *
from ML import *

plt.figure(figsize=(17,7))
#star = spec(findstar('HD41398_20201218_070204_M_V85000_log.fits'),SNR='bestMF')
star = spec('ALS15210',SNR='bestMF')
star.waveflux(3950,6850)
plt.plot(star.wave,star.flux,'orange',lw=.7,label='original')

star.flux = cosmicML(star.wave,star.flux,method='zscore',sigclip=3.0,iter=3)
flux_kern = cosmicML(star.wave,star.flux,method='zscore',sigclip=1.3,iter=3)
star.flux = np.where(star.flux > 1.01,flux_kern,star.flux)
plt.plot(star.wave,star.flux,'--b',lw=.7,label='finalML')

star.degrade(resol=5000)
star.resamp(10*0.02564975,3950,6850)
plt.plot(star.wave,star.flux+1,'b',lw=.7,label='degraded')

star.waveflux(3950,6850)
star.cosmic(method='kernel')
plt.plot(star.wave,star.flux,'--g',lw=.7,label='finalOri')

star.degrade(resol=5000)
star.resamp(10*0.02564975,3950,6850)
plt.plot(star.wave,star.flux+1.5,'g',lw=.7,label='degraded')

plt.rcParams['figure.dpi'] = 300
plt.legend(); plt.tight_layout()

plt.show()
