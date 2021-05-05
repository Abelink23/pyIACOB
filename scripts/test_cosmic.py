import sys; sys.path.append('../')

from spec import *

plt.close('all')
plt.figure()
star = spec('HD41398',SNR='bestMF')
wave,flux,_ = star.waveflux(4500,6380)

sigclip = 1.5

sigma = .2
# OR
sigma = 2*np.mean(star.wave)/(2.35482*float(star.resolution))

x = np.arange(-5*sigma,5*sigma+star.dx,star.dx)
gauss = f_gaussian(x,sigma)
kernel = gauss/np.trapz(gauss)

# computational convolution
convoluted = 1 + convolve(flux-1,kernel,mode='same')

plt.plot(x,kernel)

plt.figure()

plt.plot(wave,flux,label='original')
plt.plot(wave,convoluted-0.05,label='convolved')

flux_norm = flux/convoluted

plt.plot(wave,flux_norm+0.05,'gray',label='ori./conv.')

std = np.std(flux_norm)
flux_cleaned = np.where(flux_norm>1+sigclip*std,np.nan,flux)

nans = np.isnan(flux_cleaned); x = lambda z: z.nonzero()[0]
flux_cleaned[nans]= np.interp(x(nans),x(~nans),flux_cleaned[~nans])

plt.plot(wave,flux_cleaned-0.1,label='final')

plt.legend()
plt.show(block=False)
