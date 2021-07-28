from spec import *
from RV import *

import random

def findSB(table=None, RVcorr=True):

    '''
    Function to plot all the available spectra from a star in order to visually
    identify spectroscopic binaries.

    Parameters
    ----------
    table : str
        Name of the input table contaning the list of stars to analyze.

    RVcorr : boolean
        True if RV0 corrections is applied before plotting the spectra.

    redo : str, optional
        Coma separated string with the list of stars for which repeat the analysis.

    Returns
    -------
    Nothing, but the plot is created.
    '''

    if name == None: table = findtable('IACOB_O9BAs_SNR20.fits')

    quit = ''
    for row in table:

        if quit == 'quit': break

        name = row['Name'].strip()

        skip = input('%s - Hit return to continue, type "s" to skip: ' % name)
        if skip == 's': continue

        if row['SpT_code'] < 2: list = 'rv_Os.lst'
        elif row['SpT_code'] >= 2: list = 'rv_Bs.lst'

        spectra = findstar(name,SNR=20)
        if len(spectra) > 15: spectra = random.sample(spectra,15)

        for spectrum,i in zip(spectra,range(len(spectra))):
            star = spec(spectrum,SNR='best')

            if RVcorr == True:
                star.rv0 = RV0(list,star.spectrum,func='g',ewcut=30,tol=150)
            star.waveflux(4530,6700) # Applies the rv0 correction
            star.cosmic(sigclip=1.2)
            star.degrade(20000)
            star.flux = star.flux+i*0.05
            star.plotline('4552.622,5411.52,5875.64,6562.80',width=25)
            plt.legend().remove()

        print('\n\t',name,'\n')
        next = input('\nHit return to move to the next star.')
        plt.close()
