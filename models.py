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


def trackmist(mass=None, av=0.0, vini_vcrit=0.4):

# To generate individual tables from MIST track files and a master FITS with all tracks
# 1) Edit every new table from MIST removing the first lines before the column names.
# 2) Run the following code:
#    t_master = Table()
#    path = os.path.expanduser('~')+'/Documents/MIST/'
#    for file in os.listdir(path):
#        if file.endswith('.cmd'):
#            try: t_mist = Table.read(path+file, format='ascii')
#            except: print(file,' could not be read. Check the file.')
#            mass = str(round(t_mist['star_mass'][0],1)).replace('.','')
#            digit = 4-len(mass); mass = '0'*digit+mass
#            t_mist = t_mist[(t_mist['phase']>=0) & (t_mist['phase']<=4)]
#    
#            hdu = fits.BinTableHDU(data=t_mist.filled(np.nan))
#            hdu.writeto(path+'TRAC_FeH0_%sMsol_Av00_V00.fits' % mass,overwrite=True)
#    
#            t_master = vstack([t_master,t_mist],join_type='outer')
#    
#    hdu = fits.BinTableHDU(data=t_master.filled(np.nan))
#    hdu.writeto(path+'TRAC_FeH0_08-120_Av00_V00.fits' ,overwrite=True)

    '''
    Function to retrieve a specific track from those downloaded from MIST webtool.

    NOTE: Tracks from MIST are downloaded from:
          https://waps.cfa.harvard.edu/MIST/interp_tracks.html
          MIST version used: 1.2
          Synthetic Photometry: UBV(RI)c + 2MASS + Kepler + Hipparcos + Gaia (DR2/MAW/EDR3) + Tess

    NOTE: Remove the first lines before the column names on every raw table downloaded from MIST.

    Parameters
    ----------
    mass : int/float, optional
        Enter the mass in M/M_sun of the track you want to retrieve.
        If None as input, all the tracks will be selected.
        Available masses are: 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 1.7, 
        2.0, 2.5, 3, 4, 5, 7, 9, 12, 15, 20, 25, 32, 40, 60, 85, 120.

    av : float, optional
        Enter the extinction (Av) of the isochrone to retrieve. Default is 1.0.

    vini_vcrit : float, optional
        Enter the initial v/v_crit value [0.0/0.4]. Default is 0.4.

    Returns
    -------
    MIST track.
    '''

    if av < 1.0: Av = str(av).replace('.','')
    else: Av = int(av*10)

    vini_vcrit = str(vini_vcrit).replace('.','')

    mass_list = [0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 1.7, 2.0, 2.5, 3, 4, 5, 7, 9, 12, 15, 20, 25, 32, 40, 60, 85, 120]

    if mass != None and mass in mass_list:
        mass = str(float(mass)).replace('.',''); digit = 4-len(mass)
        mass = '0'*digit+mass

        trk_mist = Table.read(modeldir + 'MIST/TRACKS/TRAC_FeH0_Av%s_V%s/TRAC_FeH0_%sMsol_Av%s_V%s.fits' % \
            (Av,vini_vcrit,mass,Av,vini_vcrit), format='fits')

    else:
        trk_mist = Table.read(modeldir + 'MIST/TRACKS/TRAC_FeH0_Av%s_V%s/TRAC_FeH0_08-120_Av%s_V%s.fits' % \
            (Av,vini_vcrit,Av,vini_vcrit))

    return trk_mist


def trackgene(mass=None, vini_vcrit=0.4):

    '''
    Function to retrieve a specific track from Ekstrom et al. (2012).
    
    NOTE: Tracks and isochrones from Geneva are downloaded from:
          https://obswww.unige.ch/Research/evol/tables_grids2011/

    Parameters
    ----------
    mass : int/float, optional
        Enter the mass in M/M_sun of the track you want to retrieve.
        If None as input, all the tracks will be selected.

    vini_vcrit : float, optional
        Enter the initial v/v_crit value [0.0/0.2/0.4]. Default is 0.4.

    Returns
    -------
    Geneva track.

    Notes 'line' column
    -------------------
    1: ZAMS
    2-84: H burning (first part)
    85: minimum of Teff on the MS
    86-109: overall contraction phase before the end of the MS
    110: Turn-off
    111-189: HR diagram crossing and/or pre-He-b core contraction
    190: beginning of He burning
    191-209: He burning (first part)
    210-350: blue loop (if any, maximal extension on point 280)
    351-369: He burning (second part)
    370: core He exhaustion
    371-399: C burning
    400: last model
    '''

    if mass == None:
        print('ERROR: No mass given. Please enter a mass in M/M_sun.')
        return

    if not vini_vcrit in [0.0, 0.2, 0.4]:
        print('WARNING: Geneva tracks are only available for v/vcrit = 0.0, 0.2, and 0.4')
        vini_vcrit = min([0.0, 0.2, 0.4], key=lambda x:abs(x-vini_vcrit))
        print('Choosing %s as the closest value.' % vini_vcrit)

    if vini_vcrit in [0.0, 0.4]:
        mass_list = [0.8, 0.9, 1.0, 1.1, 1.25, 1.35, 1.5, 1.7, 2.0, 2.5, 3, 4, 5, 7, 9, 12, 15, 20, 25, 32, 40, 60, 85, 120]
    elif vini_vcrit == 0.2:
        mass_list = [20, 25, 32, 40, 60, 85, 120]

    vini_vcrit = str(vini_vcrit).replace('0.','')

    if type(mass) is not str and not mass in mass_list:
        print('WARNING: Mass not in list %s' % mass_list)
        mass = min(mass_list, key=lambda x:abs(x-mass))
        print('Choosing %s as the closest value.' % mass)

    # if mass is a round number, it is turned into an integer
    if type(mass) is float and mass.is_integer():
        mass = str(int(mass))
    elif type(mass) is not str and mass in mass_list:
        mass = str(mass).replace('.','p')

    digit = 3-len(mass)
    mass = '0'*digit+mass

    trk_geneva = Table.read(modeldir + 'GENEVA/M%sZ14V%s.dat' % (mass,vini_vcrit), format='ascii', data_start=2, delimiter=' ')

    # FROM Gonzalo
    trk_geneva.rename_columns(['lg(Teff)','lg(L)'],['log_Teff','log_L'])
    #trk_geneva['L'] = (10**trk_geneva['lg(L)'])
    trk_geneva['Teff'] = (10**trk_geneva['log_Teff'])/1e4
    trk_geneva['log_LLsol'] = trk_geneva['log_L'] - np.log10(trk_geneva['mass']) # NOT SURE ABOUT THIS ONE
    trk_geneva['log_g'] = 4*trk_geneva['log_Teff'] + np.log10(trk_geneva['mass']) - trk_geneva['log_L'] - 10.61
    trk_geneva['log_Lspec'] = 4*trk_geneva['log_Teff'] - trk_geneva['log_g'] - 10.61
    trk_geneva['He'] = trk_geneva['4He_surf']/4/trk_geneva['1H_surf']

    # missing to implement a table with everything
    #else:
    #    trk_geneva = Table.read(modeldir + 'MIST/TRACKS/TRAC_FeH0_Av%s_V%s/TRAC_FeH0_08-120_Av%s_V%s.fits' % \
    #        (Av,vini_vcrit,Av,vini_vcrit))

    return trk_geneva


def trackbonn(mass=None, vini_vcrit=0.0):

    '''
    Function to retrieve a specific track from Brott et al. (2011).

    NOTE: Tracks and isochrones are downloaded from:
          https://cdsarc.u-strasbg.fr/viz-bin/qcat?J/A+A/530/A115

    Parameters
    ----------

    mass : int/float, optional
        Enter the mass in M/M_sun of the track you want to retrieve.
        Available masses are: 5, 7, 9, 10, 12, 15, 20, 25, 30, 35, 40, 50, 60.

    vini_vcrit : int, optional
        Enter the initial v/v_crit value. Default is 0.0.
        NOTE: The closest vini will be selected for the given mass.

    Returns
    -------
    BONN track.
    '''

    mass_list = [5, 7, 9, 10, 12, 15, 20, 25, 30, 35, 40, 50, 60]

    while mass not in mass_list:
        print('Mass %s not in list %s' % (str(mass),mass_list))
        mass = int(input('Pick a new mass: '))

    # find all models starting with f+mass
    path = modeldir + 'BONN/'
    models = [f for f in os.listdir(path) if f.startswith('f'+str(mass))]

    # find the vini for each model
    vini_list = [int(m.split('-')[1].split('.mw.fits')[0]) for m in models]
    # calculate the vini/vcrit for each vini
    vcrit_list = [Table.read(modeldir+'BONN/f%s-%s.mw.fits' % (mass,v), format='fits')['Vcrit'][0] for v in vini_list]

    # pick the vini to achieve the nearest v/vcrit
    vini = min(vini_list, key=lambda x:abs(x-vini_vcrit*vcrit_list[vini_list.index(x)]))
    if abs(vini-vini_vcrit*vcrit_list[vini_list.index(vini)]) > 20:
        print('WARNING: The closest vini is more than 20 km/s away from the desired v/vcrit')
    print('Closest v/v_crit is', round(vini/vcrit_list[vini_list.index(vini)],2), 'for vini =', vini, 'km/s')

    trk_brott = Table.read(modeldir+'BONN/f%s-%s.mw.fits' % (mass,vini), format='fits')

    # log(X/H) + 12
    for elem in ['He','C','N','O','Mg','Si']:
        trk_brott[elem] = 10**(trk_brott['eps(%s)' % elem]-12)

    trk_brott.rename_column('logL','log_L')
    trk_brott['log_Teff'] = np.log10(trk_brott['Teff'])
    trk_brott['log_g'] = 4*trk_brott['log_Teff'] + np.log10(trk_brott['Mass']) - trk_brott['log_L'] - 10.61
    trk_brott['log_Lspec'] = 4*trk_brott['log_Teff'] - trk_brott['log_g'] - 10.61

    return trk_brott
