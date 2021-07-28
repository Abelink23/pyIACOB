from db import *

def isomist(myr=None, logmyr=None, av=1.0, vr=0.4):

    '''
    Function to retrieve a specific isochrone from the MIST.

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

    myr_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,\
    32,35,38,41,45,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,\
    220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400]
    logmyr_list = [6.0,6.301,6.477,6.602,6.699,6.778,6.845,6.903,6.954,7.0,7.041,7.079,\
    7.114,7.146,7.176,7.204,7.23,7.255,7.279,7.301,7.342,7.38,7.415,7.447,7.477,7.505,\
    7.544,7.58,7.613,7.653,7.699,7.778,7.845,7.903,7.954,8.0,8.041,8.079,8.114,8.146,\
    8.176,8.204,8.23,8.255,8.279,8.301,8.322,8.342,8.362,8.38,8.398,8.415,8.431,8.447,\
    8.462,8.477,8.491,8.505,8.519,8.531,8.544,8.556,8.568,8.58,8.591,8.602]

    if av < 1.0: Av = str(av).replace('.','')
    else: Av = int(av*10)

    vr = str(vr).replace('.','')

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

    # NOTE: Edit every new table from MIST removing the first lines before the column names.
    t_mist = Table.read(mistdir + 'ISOCHRONES/ISOC_FeH0_%sMyr_Av%s_V%s.fits' % \
        (ranage,Av,vr),format='fits')

    t_mist = t_mist[t_mist['log10_isochrone_age_yr'] == logage]

    return t_mist


def trackmist(mass=None, av=0.0, vr=0.4):

    '''
    Function to retrieve a specific track from the MIST.

    Parameters
    ----------
    mass : int/float, optional
        Enter the mass in M/M_sun of the track you want to retrieve.
        If None as input, all the tracks will be selected.

    av : float, optional
        Enter the extinction (Av) of the isochrone to retrieve. Default is 1.0.

    vr : float, optional
        Enter the initial v/v_crit value [0.0/0.4]. Default is 0.4.

    Returns
    -------
    MIST isochrone.
    '''

    if av < 1.0: Av = str(av).replace('.','')
    else: Av = int(av*10)

    vr = str(vr).replace('.','')

    mass_list = [.8,.9,1,1.1,1.2,1.3,1.5,1.7,2,2.5,3,4,5,7,9,12,15,20,25,32,40,60,85,120]

    if mass != None and mass in mass_list:
        mass = str(float(mass)).replace('.',''); digit = 4-len(mass)
        mass = '0'*digit+mass

        t_mist = Table.read(mistdir + 'TRACKS/TRAC_FeH0_Av%s_V%s/TRAC_FeH0_%sMsol_Av%s_V%s.fits' % \
            (Av,vr,mass,Av,vr), format='fits')

    else:
        t_mist = Table.read(mistdir + 'TRACKS/TRAC_FeH0_Av%s_V%s/TRAC_FeH0_08-120_Av%s_V%s.fits' % \
            (Av,vr,Av,vr))

    return t_mist

## To generate individual tables from MIST imput files (tracks), and a master FITS
## with all the tracks
## NOTE: Edit every new table from MIST removing the first lines before the column names.
#t_master = Table()
#path = os.path.expanduser('~')+'/Documents/MIST/'
#for file in os.listdir(path):
#    if file.endswith('.cmd'):
#        try: t_mist = Table.read(path+file,format='ascii')
#        except: print(file,' could not be read. Check the file.')
#        mass = str(round(t_mist['star_mass'][0],1)).replace('.','')
#        digit = 4-len(mass); mass = '0'*digit+mass
#        t_mist = t_mist[(t_mist['phase']>=0) & (t_mist['phase']<=4)]
#
#        hdu = fits.BinTableHDU(data=t_mist.filled(np.nan))
#        hdu.writeto(path+'TRAC_FeH0_%sMsol_Av00_V00.fits' % mass,overwrite=True)
#
#        t_master = vstack([t_master,t_mist],join_type='outer')
#
#hdu = fits.BinTableHDU(data=t_master.filled(np.nan))
#hdu.writeto(path+'TRAC_FeH0_08-120_Av00_V00.fits' ,overwrite=True)
