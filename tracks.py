from db import *

def trackbonn11(mass=None, vini_vcrit=0.0):

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
    Brott+11 track.
    '''

    mass_list = [5, 7, 9, 10, 12, 15, 20, 25, 30, 35, 40, 50, 60]

    if type(mass) is not str and not mass in mass_list:
        msg.warn('Mass %s not in list %s' % (str(mass),mass_list))
        mass = min(mass_list, key=lambda x:abs(x-mass))
        msg.info('Choosing %s as the closest value.' % mass)

    mass = str(int(mass)); digit = 3-len(mass)
    mass = '0'*digit+mass

    # find all models starting with f+mass
    path = tracksdir + 'Brott11/'
    tracks = [f for f in os.listdir(path) if f.startswith('M'+str(mass))]

    # find the vini for each model
    vini_list = [m.split('V')[1].split('Av')[0] for m in tracks]
    # calculate the vini/vcrit for each vini
    vcrit_list = [Table.read(tracksdir+'Brott11/M%sZ014V%sAv00.fits' % (mass,v), format='fits')['Vcrit'][0] for v in vini_list]

    # pick the vini to achieve the nearest v/vcrit
    vini_list = [int(i) for i in vini_list]
    vini = min(vini_list, key=lambda x:abs(x-vini_vcrit*vcrit_list[vini_list.index(x)]))
    if abs(vini-vini_vcrit*vcrit_list[vini_list.index(vini)]) > 20:
        msg.warn('The closest vini is more than 20 km/s away from the desired v/vcrit')
    msg.info('Closest v/v_crit is %f for vsini = %f km/s' % (round(vini/vcrit_list[vini_list.index(vini)],2),vini))

    trk_brott = Table.read(tracksdir+'Brott11/M%sZ014V%sAv00.fits' % (mass,f"{int(vini):03d}"), format='fits')

    # log(X/H) + 12
    for elem in ['He','C','N','O','Si']:
        trk_brott[elem] = 10**(trk_brott['eps(%s)' % elem]-12)

    trk_brott.rename_column('logL','log_L')
    trk_brott['log_Teff'] = np.log10(trk_brott['Teff'])
    trk_brott['log_g'] = 4*trk_brott['log_Teff'] + np.log10(trk_brott['Mass']) - trk_brott['log_L'] - 10.61
    trk_brott['log_Lspec'] = 4*trk_brott['log_Teff'] - trk_brott['log_g'] - 10.61

    return trk_brott


def trackmist(mass=None, metallicity=0.014, vini_vcrit=0.4, av=0.0):
    '''
    Function to retrieve a specific track from those downloaded from MIST webtool.

    NOTE: Tracks from MIST are downloaded from:
          https://waps.cfa.harvard.edu/MIST/interp_tracks.html
          MIST version used: 2.5
          List of masses: 9, 12, 15, 20, 25, 30, 35, 40, 45, 60, 85
          Extinction: 0
          Synthetic Photometry: UBV(RI)c + 2MASS + Kepler + Hipparcos + Gaia (DR2/MAW/EDR3) + Tess

    Parameters
    ----------
    mass : int/float, optional
        Enter the mass in M/M_sun of the track you want to retrieve.
        If None as input, all the tracks will be selected.
        Available masses are: 9, 12, 15, 20, 25, 30, 35, 40, 45, 60, 85.

    metallicity : float, optional
        Enter the metallicity value as 0.###. Default is 0.014 (solar)

    vini_vcrit : float, optional
        Enter the initial v/v_crit value. Default is 0.4.
        Available values: 0.0, 0.4.

    av : float, optional
        Enter the extinction (Av) of the isochrone to retrieve. Default is 0.0.
        Available values: 0.0

    Returns
    -------
    MIST track.
    '''

    av_list = [0.0]
    vini_vcrit_list = [0.0, 0.4]
    mass_list = [9, 12, 15, 20, 25, 30, 35, 40, 45, 60, 85]

    # Choose Av
    if not av in av_list:
        msg.warn('Selected MIST Av not within: %s' % av_list)
        av = min(av_list, key=lambda x:abs(x-av))
        msg.info('Choosing %s as the closest value.' % av)
    if av < 1.0:
        print(av)
        Av = str(av).replace('.','')
    else:
        Av = int(av*10)

    # Choose vini_vcrit
    if not vini_vcrit in vini_vcrit_list:
        msg.warn('Selected MIST vini/vcrit not within: %s' % vini_vcrit_list)
        vini_vcrit = min(vini_vcrit_list, key=lambda x:abs(x-vini_vcrit))
        msg.info('Choosing %s as the closest value.' % vini_vcrit)
    vini_vcrit = str(vini_vcrit).replace('.','')

    # Choose mass
    if type(mass) is not str and not mass in mass_list:
        msg.warn('Selected MIST mass not within: %s' % mass_list)
        mass = min(mass_list, key=lambda x:abs(x-mass))
        msg.info('Choosing %s as the closest value.' % mass)
    mass = str(int(mass)); digit = 3-len(mass)
    mass = '0'*digit+mass

    metallicity = str(metallicity)[3:]

    trk_mist = Table.read(tracksdir + 'MIST/M%sZ%sV%sAv%s.fits' % \
        (mass, metallicity, vini_vcrit, Av), format='fits')

    return trk_mist


def trackgene12(mass, vini_vcrit=0.4):
    '''
    Function to retrieve a specific track from Ekstrom et al. (2012).

    NOTE: Tracks and isochrones from Geneva are downloaded from:
          https://obswww.unige.ch/Research/evol/tables_grids2011/

    Parameters
    ----------
    mass : int/float
        Enter the mass in M/M_sun of the track you want to retrieve.

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

    vini_vcrit_list = [0.0, 0.2, 0.4]

    # Choose vini_vcrit
    if not vini_vcrit in vini_vcrit_list:
        msg.warn('Selected Genec vini/vcrit not within: %s' % vini_vcrit_list)
        vini_vcrit = min(vini_vcrit_list, key=lambda x:abs(x-vini_vcrit))
        msg.info('Choosing %s as the closest value.' % vini_vcrit)

    if vini_vcrit in [0.0, 0.4]:
        mass_list = [0.8, 0.9, 1.0, 1.1, 1.25, 1.35, 1.5, 1.7, 2.0, 2.5, 3, 4, 5, 7, 9, 12, 15, 20, 25, 32, 40, 60, 85, 120]
    elif vini_vcrit == 0.2:
        mass_list = [20, 25, 32, 40, 60, 85, 120]
    vini_vcrit = str(vini_vcrit).replace('0.','')

    # Choose mass
    if type(mass) is not str and not mass in mass_list:
        msg.warn('Selected Genec mass not within: %s' % mass_list)
        mass = min(mass_list, key=lambda x:abs(x-mass))
        msg.info('Choosing %s as the closest value.' % mass)

    # if mass is a round number, it is turned into an integer
    if type(mass) is float and mass.is_integer():
        mass = str(int(mass))
    elif type(mass) is not str and mass in mass_list:
        mass = str(mass).replace('.','p')

    digit = 3-len(mass)
    mass = '0'*digit+mass

    trk_geneva = Table.read(tracksdir + 'Ekstrom12/M%sZ14V%s.fits' % (mass,vini_vcrit), format='fits')

    # FROM Gonzalo
    trk_geneva.rename_columns(['lg(Teff)','lg(L)'],['log_Teff','log_L'])
    #trk_geneva['L'] = (10**trk_geneva['lg(L)'])
    trk_geneva['Teff'] = (10**trk_geneva['log_Teff'])/1e4
    trk_geneva['log_LLsol'] = trk_geneva['log_L'] - np.log10(trk_geneva['mass']) # NOT SURE ABOUT THIS ONE
    trk_geneva['log_g'] = 4*trk_geneva['log_Teff'] + np.log10(trk_geneva['mass']) - trk_geneva['log_L'] - 10.61
    trk_geneva['log_Lspec'] = 4*trk_geneva['log_Teff'] - trk_geneva['log_g'] - 10.61
    trk_geneva['He'] = trk_geneva['4He_surf']/4/trk_geneva['1H_surf']

    return trk_geneva


def trackgene26(mass, vini_vcrit=0.4, model='magn', tams='FS'):
    '''
    Function to retrieve a specific track from Sciarini et al. (2026).

    NOTE: Tracks and isochrones from Geneva are downloaded from:
          https://zenodo.org/records/22082876/

    Parameters
    ----------
    mass : int/float
        Enter the mass in M/M_sun of the track you want to retrieve.

    vini_vcrit : float, optional
        Enter the initial v/v_crit value [0.0/0.2/0.4]. Default is 0.4.

    model : str, optional
        Enter 'magn' to retrieve a magnetic model. Alternatively, use
        'hydro' to retrieve a hydrostatic model. Default is 'magn'.

    tams : str, optional
        Enter 'FS' if you want to retrieve a track computed to match the density
        drop of the Full-Sample in de Burgos et al. (2025). Alternatively,
        enter 'SR' to retrieve a track computed to much the density drop of the
        Slow Rotating sample in de Burgos et al. (2025). Default is FS.

    Returns
    -------
    Geneva track.

    Notes 'line' column
    -------------------
    1: ZAMS
    2-84: H burning
    85: minimum of Teff on the MS (main sequence hook)
    86-109: H burning
    110: Turn-off
    111-189: HRD crossing and/or pre-He-burning core contraction
    190-400: last computed model (the models correspond to the limit of the de Burgos+25
            sample they are stopped once the effective temperature drops below 14 kK).
    '''

    vini_vcrit_list = [0.0, 0.05, 0.10, 0.15, 0.25, 0.40, 0.60]
    mass_list = [12, 15, 20, 25, 32, 40]

    # Choose vini_vcrit
    if not vini_vcrit in vini_vcrit_list:
        msg.warn('Selected Genec vini/vcrit not within: %s' % vini_vcrit_list)
        vini_vcrit = min(vini_vcrit_list, key=lambda x:abs(x-vini_vcrit))
        msg.info('Choosing %s as the closest value.' % vini_vcrit)
    vini_vcrit = str(vini_vcrit*10).replace('.','')

    # Choose mass
    if type(mass) is not str and not mass in mass_list:
        msg.warn('Selected Genec mass not within: %s' % mass_list)
        mass = min(mass_list, key=lambda x:abs(x-mass))
        msg.info('Choosing %s as the closest value.' % mass)

    # if mass is a round number, it is turned into an integer
    if type(mass) is float and mass.is_integer():
        mass = str(int(mass))
    elif type(mass) is not str and mass in mass_list:
        mass = str(mass).replace('.','p')

    digit = 3-len(mass)
    mass = '0'*digit+mass

    if model not in ['magn','hydro']:
        msg.error("The model input parameter must be 'magn' or 'hydro'")
    if tams not in ['FS','SR']:
        msg.error("The tams input parameter must be 'FS' or 'SR'")

    trk_geneva = Table.read(tracksdir + 'Sciarini26/M%sZ14V%sAv00_%s_%s.fits' % (mass,vini_vcrit,model,tams), format='fits')

    trk_geneva.rename_columns(['lg(Teff)','lg(L)'],['log_Teff','log_L'])
    #trk_geneva['L'] = (10**trk_geneva['lg(L)'])
    trk_geneva['Teff'] = (10**trk_geneva['log_Teff'])/1e4
    trk_geneva['log_LLsol'] = trk_geneva['log_L'] - np.log10(trk_geneva['mass']) # NOT SURE ABOUT THIS ONE
    trk_geneva['log_g'] = 4*trk_geneva['log_Teff'] + np.log10(trk_geneva['mass']) - trk_geneva['log_L'] - 10.61
    trk_geneva['log_Lspec'] = 4*trk_geneva['log_Teff'] - trk_geneva['log_g'] - 10.61
    trk_geneva['He'] = trk_geneva['4He_surf']/4/trk_geneva['1H_surf']

    return trk_geneva


