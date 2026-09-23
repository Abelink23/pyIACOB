from db import *

def isomist(myr=None, logmyr=None, av=1.0, vini_vcrit=0.4):

    '''
    Function to retrieve a specific isochrone from those downloaded from MIST webtool.

    NOTE: Isochrones from MIST are downloaded from:
          https://waps.cfa.harvard.edu/MIST/interp_isos.html
          MIST version used: 1.2
          Two "list of ages" with Log10 Scale are used to create the available ages (see function below)
          Synthetic Photometry: UBV(RI)c + 2MASS + Kepler + Hipparcos + Gaia (DR2/MAW/EDR3) + Tess

    NOTE: Remove the first lines before the column names on every raw table downloaded from MIST.

    Parameters
    ----------
    myr : int/float, optional
        Enter the age in Myr of the isochrone you want to retrieve.

    logmyr : int/float, optional
        Enter the age as log10(Myr) of the isochrone you want to retrieve.

    av : float, optional
        Enter the extinction (Av) of the isochrone to retrieve. Default is 1.0.

    vr : float, optional
        Enter the initial v/v_crit value [0.0/0.4]. Default is 0.4.

    Returns
    -------
    MIST isochrone.
    '''

    logmyr_list = [6.0, 6.301, 6.477, 6.602, 6.699, 6.778, 6.845, 6.903, 6.954, 7.0, 7.041, 7.079, \
    7.114, 7.146, 7.176, 7.204, 7.23, 7.255, 7.279, 7.301, 7.342, 7.38, 7.415, 7.447, 7.477, 7.505, \
    7.544, 7.58, 7.613, 7.653, 7.699, 7.778, 7.845, 7.903, 7.954, 8.0, 8.041, 8.079, 8.114, 8.146, \
    8.176, 8.204, 8.23, 8.255, 8.279, 8.301, 8.322, 8.342, 8.362, 8.38, 8.398, 8.415, 8.431, 8.447, \
    8.462, 8.477, 8.491, 8.505, 8.519, 8.531, 8.544, 8.556, 8.568, 8.58, 8.591, 8.602]

    myr_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 26, 28, 30, \
    32, 35, 38, 41, 45, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, \
    220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400]

    if av < 1.0:
        Av = str(av).replace('.','')
    else:
        Av = int(av*10)

    vini_vcrit = str(vini_vcrit).replace('.','')

    if myr != None:
        if not myr in myr_list:
            print('Age not in list %s' % myr_list)
            myr = int(input('Pick a new age from the list: '))

        logage = round(np.log10(myr*1e6),3)

    if logmyr != None:
        logage = min(logmyr_list, key=lambda x:abs(x-logmyr))
        if abs(logage-logmyr) > 0.3:
            print('Difference to closes isochrone is grater than 0.3 (~2Myr)')

    if logage < 7.676: ranage = '1-45'
    else: ranage = '50-300'

    iso_mist = Table.read(modeldir + 'MIST/ISOCHRONES/ISOC_FeH0_%sMyr_Av%s_V%s.fits' % \
        (ranage,Av,vini_vcrit),format='fits')

    iso_mist = iso_mist[iso_mist['log10_isochrone_age_yr'] == logage]

    return iso_mist

