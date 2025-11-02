from rv import *

import random

def findSB(ID, snr, zone='full', degrade=None, RV0lines='rv_Bs.lst', vspace=0.07, c=None):

    '''
    Function to plot all the available spectra from a star in order to visually
    identify spectroscopic binaries.

    Parameters
    ----------
    ID : str
        ID of the source to be analysed.

    snr : str/int
        SNR criteria for selecting the input spectra. See db.findstar().

    zone : str, optional
        Area of the spectrum where to look for binarity signs.
        Options are: full (default) / 4300 / Si / SiIII / Hb / HeI / Ha

    degrade : int/float, optional
        Resolution used to degrade all the spectra in order to do an even comparison.
        Default is no degradation.

    RV0lines : str, optional
        Input list of lines used to match all the radial velocities.
        Default is 'rv_Bs.lst'

    vspace : int/float, optional
        Vertical space between the spectra. Default is 0.07

    Returns
    -------
    Nothing, but the plot is created.
    '''

    spectra = findstar(ID, snr=snr)

    n_max = 15 # Maximum number of spectra to be plotted if there are more than 15

    if len(spectra) > n_max:
        n = input('Number of spectra is >15, do you want to take a random number of them? [#/no]: ')
        if n != 'no' and n != '':
            spectra = random.sample(spectra, int(n))
        elif n == '':
            print('Plotting 15 random spectra.')
            spectra = random.sample(spectra, n_max)

    plt.figure(figsize=(11, 5))
    for i,j in zip(spectra, range(len(spectra))):

        rv0 = RV0(RV0lines, i.split('/')[-1].split('\\')[-1], ewcut=50, width=15, tol=100, func='g')[0]

        sp = spec(i.split('/')[-1].split('\\')[-1], rv0=rv0)

        snr_fits = sp.snr; sp.snrcalc('V')
        print('\n', sp.filename, 'FITS:', snr_fits, ' V:', round(sp.snr))

        if zone == 'full': sp.waveflux(3900,6800) # Full spectra
        elif zone == '4300': sp.waveflux(4150,4300) # 4300
        elif zone == 'Si': sp.waveflux(4407,4600) # Si III
        elif zone == 'SiIII': sp.waveflux(4540,4580) # Si III zoom
        elif zone == 'Hb': sp.waveflux(4800,4920) # Hb
        elif zone == 'HeI': sp.waveflux(5865,5885) # HeI
        elif zone == 'Ha': sp.waveflux(6530,6600) # Halpha

        if degrade is not None:
            sp.degrade(degrade)
        #sp.resamp(5*sp.dlam)

        ax = plt.gca()
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
        plt.tick_params(axis='x', which='both', top=True, direction='in')
        plt.plot(sp.wave, sp.flux + j*vspace, lw=.5, c=c)

    plt.tight_layout()
    plt.show(block=False)
