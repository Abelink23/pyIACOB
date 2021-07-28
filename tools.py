from db import *

import scipy.constants as cte

from astroquery.nist import Nist
from astroquery.atomic import AtomicLineList

from astropy.time import Time
from astropy.stats import sigma_clip
from astropy.coordinates import EarthLocation


#===============================================================================
# Generate table with Gaia eDR3 data plus corrected parallax zero point offset.
#===============================================================================
def table_zp_edr3(table, ra, dec, search_radius=0.5):

    '''
    Function to query input list of RA/DEC coordinates in Gaia eDR3 and for the
    output sources adds a column with the corrected parallax zero point offset.

    Parameters
    ----------
    table : str
        Name of the input table containing the coordinates.

    ra : str
        Column name of the RA.

    dec : str
        Column name of the DEC.

    search_radius : int/float
        Search radius for Gaia query in arcseconds. Default is 0.5.

    Returns
    -------
    Table containing the output query of Gaia eDR3 including the corrected parallax offset.
    '''

    import sys
    sys.path.append(os.path.expanduser('~')+'/MEGA/PhD/programs/python/edr3_zp')
    import zpt
    zpt.load_tables()
    from astroquery.gaia import Gaia

    sr = round(search_radius/60/60,7)

    table = findtable(table)

    first = True
    for ra_i,dec_i in zip(table[ra],table[dec]):
        job = Gaia.launch_job("select TOP 3 * FROM gaiaedr3.gaia_source "
                    "WHERE 1=CONTAINS(POINT('ICRS',ra,dec), "
                    "CIRCLE('ICRS',%f,%f,%f))" % (ra_i,dec_i,sr)).get_results()
        if len(job) == 0: continue
        elif len(job) > 1: job.sort('phot_g_mean_mag')

        if first == True: table = job; first = False
        else: table.add_row(job[0])

    if first == True:
        print('No objects were found for the input coordinates.')
        return None

    table = table[table['astrometric_params_solved']>3]
    zpvals = zpt.get_zpt(table['phot_g_mean_mag'],table['nu_eff_used_in_astrometry'],\
           table['pseudocolour'],table['ecl_lat'],table['astrometric_params_solved'])
    table.add_column(zpvals,name='zp_offset')

    table.remove_column('designation')

    table.write(maindir+'tables/zp_offset_eDR3.fits',format='fits',overwrite=True)

    return table


#===============================================================================
# Atomic lines and wavelenghts conversions.
#===============================================================================

def atline(lwl, rwl, elements=None, source='ALL'):

    '''
    Function to retrieve spectral lines from the Atomic Line databases.

    Parameters
    ----------
    lwl : float
        Sets the start wavelenght.

    rwl : float
        Sets the end wavelenght.

    elements : str, optional
        Enter a string with the elements associated to the lines.
        If 'OB' it will pick the most characteristic lines of OB stars.

    source : str, optional
        Choose between 'ALL' and 'NIST' databases. Default is 'ALL'

    Returns
    -------
    List with the spectral lines.
    '''

    wavelength_range = (lwl*u.Angstrom,rwl*u.Angstrom)
    if elements == 'OB':
        elements = 'H I, He I-II, O I-IV, C I-IV, Ne I-III, Fe I-III, N I-IV, \
        Si I-IV, Mg I-IV, S I-IV, V I-II, Cr I-II, Ni I-II, Sc I-II, Ti I-II, \
        Ca I-II, Na I-II'

    if source == 'ALL':
        columns = ['spec', 'type', 'conf', 'term', 'angm', 'prob']
        list = AtomicLineList.query_object(wavelength_range,wavelength_type='Air',\
        wavelength_accuracy=20,element_spectrum=elements,output_columns=columns)

    elif source == 'NIST':
        list = Nist.query(lwl*u.Angstrom,rwl*u.Angstrom,wavelength_type='vac+air',\
               linename=elements)

    return list


def atokms(dlamb, lamb0):

    '''
    Function to calculate the velocity in km/s providing delta(lambda) [AA].

    Parameters
    ----------
    dlamb : float
        Enter the delta(lambda) in angstroms.

    lamb0 : float/int
        Enter the central wavelenght for the conversion.

    Returns
    -------
    Corresponding delta in km/s.
    '''

    velocity = float(dlamb)/float(lamb0)*cte.c/1000

    return velocity


def kmstoa(dkms, lamb0):

    '''
    Function to calculate delta(lambda) [AA] providing the equivalent in km/s.

    Parameters
    ----------
    dkms : float
        Enter the delta(km/s)

    lamb0 : float/int
        Enter the central wavelenght for the conversion.

    Returns
    -------
    Corresponding delta in angstroms.
    '''

    angstroms = float(dkms)*float(lamb0)*1000/cte.c

    return angstroms


def vac2air(lamvac):

    '''
    Function to calculate the wavelength on air giving it on vacum. Both in angstroms.

    Parameters
    ----------
    lamvac : float
        Enter the wavelenght in angstroms on vacum.

    Returns
    -------
    Wavelenght on air.
    '''

    s = 1e4/float(lamvac)
    #n = 1 + 0.0000834254 + 0.02406147/(130 - s**2) + 0.00015998/(38.9 - s**2)
    n = 1 + 0.000083 + 0.02406/(130 - s**2) + 0.00016/(38.9 - s**2) # Truncated
    lamair = float(lamvac)/n

    return lamair


#===============================================================================
# COORDINATES AND DISTANCES
#===============================================================================

def sky_dist(radec1, radec2):

    '''
    Function to calculate the distance between two positions in the sky.

    Parameters
    ---------
    radec1 : str
        Enter a string with origin coordinates, e.g. '02:23:32.33 -22:13:00.3'

    radec2 : str
        Enter a string with end coordinates, e.g. '34.222 +12.222'

    Returns
    -------
    Separation in arcsec, and degrees.
    '''

    if any(i in radec1 for i in [':','h']): c1 = SkyCoord(radec1,unit=(u.hourangle,u.deg))
    else: c1 = SkyCoord(radec1,unit=u.deg)

    if any(i in radec2 for i in [':','h']): c2 = SkyCoord(radec2,unit=(u.hourangle,u.deg))
    else: c2 = SkyCoord(radec2,unit=u.deg)

    return c1.separation(c2).arcsec,c1.separation(c2).deg


def pos_ang(radec1, radec2):

    '''
    Function to calculate the possition angle in degrees.

    Parameters
    ---------
    radec1 : str
        Enter a string with origin coordinates, e.g. '02:23:32.33 -22:13:00.3'

    radec2 : str
        Enter a string with end coordinates, e.g. '34.222 +12.222'

    Returns
    -------
    Position angle in degrees.
    '''

    if any(i in radec1 for i in [':','h']): c1 = SkyCoord(radec1,unit=(u.hourangle,u.deg))
    else: c1 = SkyCoord(radec1,unit=u.deg)

    if any(i in radec2 for i in [':','h']): c2 = SkyCoord(radec2,unit=(u.hourangle,u.deg))
    else: c2 = SkyCoord(radec2,unit=u.deg)

    return c2.position_angle(c1).deg


def changecoords(list, infmt, outfmt):

    '''
    Parameters
    ----------
    list : str
        Enter the list of coordinates in .txt/.lst, or coma-separated string format.

    infmt : str
        Enter the input format, either 'hms' for hourangle (##h##m##s +-##d##m##s),
        'hms:' for hourangle with ':' delimiter, 'deg' for degrees, or 'gal' for
        galactic in degrees.

    outfmt : str
        Enter the desired output format between 'hms', 'deg', 'gal'.

    Returns
    -------
    List of coordinates in the specified output format.
    '''

    path = maindir+'lists'

    coords = []

    # To catch wrong int/float inputs:
    if type(list) == float or type(list) == int:
        print('Input cannot be int or float format. \n Exiting...'); return None

    # Lines in a lst/txt file with more information on each line:
    elif '.lst' in list or '.txt' in list: coords = findlist(list)

    # String of lines separated by coma:
    else: coords = list.split(','); coords = [i.strip() for i in coords]

    if infmt == 'hms': coords = SkyCoord(coords,unit=(u.hourangle,u.deg))
    elif format == 'deg': coords = SkyCoord(coords,unit=u.deg)

    if outfmt == 'hms': coords = coords.to_string('hmsdms')
    elif outfmt == 'hms:': coords = [re.sub('h|d|m',':',i.replace('s','')) \
                                    for i in coords.to_string('hmsdms')]
    elif outfmt == 'deg': coords = coords.to_string('decimal')
    elif outfmt == 'gal': coords = coords.galactic.to_string()

    return coords


#===============================================================================
# Others:
#===============================================================================

def exptime(exp, mag, fib=3):

    '''
    IN DEVELOPMENT

    Function to calculate the SNR for a given exposure time and magnitude (Vega).
    It is assumed an average airmass of 1.4 and gray night.

    Parameters
    ----------
    exp : int/float
        Exposure time to estimate the SNR.

    mag : int/float
        Magnitude of the source.

    Returns
    -------
    SNR
    '''

def visplot(time, target, site=None):

    '''
    Function to create visibility plot for observing runs.

    Parameters
    ----------
    time : str
        Date and time of observation (e.g. '2018-01-02 19:00').

    target : str
        Name of the target (e.g. 'HD 189733').

    Returns
    -------
    SNR
    '''

    import matplotlib.pyplot as plt

    from astroplan import FixedTarget, Observer
    from astroplan.plots import plot_airmass

    time = Time(time)
    target = FixedTarget.from_name(target)

    if site == None:
        print(EarthLocation.get_site_names())
        site = input('Plese use one of the above: ')

    apo = Observer.at_site(site)

    plot_airmass(target, apo, time, brightness_shading=True, altitude_yaxis=True)
    plt.show(block=False)


def outliers(data, iter=3, siglo=2.5, sigup=2.5):

    '''
    Function to perform a sigma clipping to a list of values.

    Parameters
    ----------
    data : list/array
        List or array-like contaning the values to clip.

    iter : int, optional
        Enter number of iterations performed by the sigma clipping. Default is 3.

    siglo : int/float
        Enter lower sigma clipping value.

    sigup : int/float
        Enter upper sigma clipping value.

    Returns
    -------
    Masked array with in values, out values, and bounds values, after clipping.
    '''

    if type(data) == list: data = np.asarray(data)

    values = sigma_clip(data,sigma_lower=siglo,sigma_upper=sigup,maxiters=iter,
            masked=False,axis=-1,cenfunc='mean',return_bounds=True)
    bounds = values[1:]; values = values[0]
    mask_in = ~np.isnan(values); mask_out = np.isnan(values)

    return data[mask_in],data[mask_out],bounds


def rv_corr(spectrum, observatory, correction):

    '''
    Function to calculate the heliocentric/barycentric radial velocity correction in km/s.

    Parameters
    ----------
    spectrum : str
        See help for db.findstar

    observatory : str
        Enter the name of the observatory.

    correction : str
        Enter the correction between barycentric or heliocentric

    Returns
    -------
    Value of the radial velocity correction.
    '''

    fits_dir = findstar(spectrum)

    if len(fits_dir) > 1:
        print('Error: More than one spectrum selected.\nExitting...')
        return None

    '''Extract the key values from the FITS'''
    hdu = fits.open(fits_dir)    # Open the fits image file
    hdu0 = hdu[0]            # Load the header list of primary header
    header0 = hdu0.header    # Read the values of the headers

    if observatory == 'INT':
        OBJECT = header0['OBJECT'].replace('.','_').split('_')[0] # Retrieve the target name
        DATE_OBS = header0['DATE-OBS']+'T'+header0['UTSTART']     # Retrieve the datetime of observation
        RA = header0['RA']                                        # Retrieve Right Ascension in hh:mm:ss.sss
        DEC = header0['DEC']                                      # Retrieve Declination in +-dd:mm:ss.ss
        EQUINOX = header0['EQUINOX']                              # Retrieve equinox
        LATITUDE = header0['LATITUDE']                            # Retrieve observatory latitude in deg
        LONGITUD = header0['LONGITUD']                            # Retrieve observatory longitude in deg
        HEIGHT = header0['HEIGHT']                                # Retrieve observatory altitude

        COORDINATES = SkyCoord(RA.strip()+' '+DEC.strip(),unit=(u.hour,u.deg),frame='icrs')
        LOC = EarthLocation.from_geodetic(float(LONGITUD),float(LATITUDE),float(HEIGHT))

    elif observatory == 'GTC':
        OBJECT = header0['OBJECT'].replace('.','_').split('_')[0] # Retrieve the target name
        DATE_OBS = header0['DATE-OBS']                            # Retrieve the datetime of observation
        RA = header0['RADEG']                                     # Retrieve Right Ascension in degrees
        DEC = header0['DECDEG']                                   # Retrieve Declination in degrees
        EQUINOX = header0['EQUINOX']                              # Retrieve equinox
        LATITUDE = header0['LATITUDE']                            # Retrieve observatory latitude
        LONGITUD = header0['LONGITUD']                            # Retrieve observatory longitude
        HEIGHT = header0['HEIGHT']                                # Retrieve observatory altitude

        COORDINATES = SkyCoord(RA*u.deg,DEC*u.deg,frame='icrs')
        LOC = EarthLocation.from_geodetic(float(LONGITUD),float(LATITUDE),float(HEIGHT))

    elif observatory == 'NOT':
        OBJECT = header0['OBJECT'].replace('.','_').split('_')[0] # Retrieve the target name
        DATE_OBS = header0['DATE-AVG']                            # Retrieve the datetime of observation
        RA = header0['RA']                                        # Retrieve Right Ascension in degrees
        DEC = header0['DEC']                                      # Retrieve Declination in degrees
        EQUINOX = header0['EQUINOX']                              # Retrieve equinox
        LATITUDE = header0['OBSGEO-Y']                            # Retrieve observatory latitude
        LONGITUD = header0['OBSGEO-X']                            # Retrieve observatory longitude
        HEIGHT = header0['OBSGEO-Z']                              # Retrieve observatory altitude

        COORDINATES = SkyCoord(RA*u.deg,DEC*u.deg,frame='icrs')
        LOC = EarthLocation.from_geocentric(float(LONGITUD),float(LATITUDE),float(HEIGHT),unit='m')

    elif observatory == 'MERCATOR':
        OBJECT = header0['OBJECT'].replace('.','_').split('_')[0] # Retrieve the target name
        DATE_OBS = header0['DATE-AVG']                            # Retrieve the datetime of observation
        RA = header0['RA']                                        # Retrieve Right Ascension in degrees
        DEC = header0['DEC']                                      # Retrieve Declination in degrees
        EQUINOX = header0['EQUINOX']                              # Retrieve equinox
        LATITUDE = header0['OBSGEO-Y']                            # Retrieve observatory latitude
        LONGITUD = header0['OBSGEO-X']                            # Retrieve observatory longitude
        HEIGHT = header0['OBSGEO-Z']                              # Retrieve observatory altitude

        COORDINATES = SkyCoord(float(RA)*u.deg,float(DEC)*u.deg,frame='icrs')
        LOC = EarthLocation.from_geocentric(float(LONGITUD),float(LATITUDE),float(HEIGHT),unit='m')

    elif observatory == 'LCO100': # NOTE THAT IT IS NOT GETTING THE AVERAGE DATETIME OF THE EXPOSURE
        OBJECT = header0['OBJECT'].replace('.','_').split('_')[0] # Retrieve the target name
        DATE_OBS = header0['DATE-OBS'] + 'T' + header0['UTSTART'] # Retrieve the datetime of observation
        RAhms = header0['RA']                                     # Retrieve Right Ascension in hms
        DECdms = header0['DEC']                                   # Retrieve Declination in dms
        RA, DEC = hms_deg(RAhms + ' ' + DECdms).split(' ')
        EQUINOX = header0['EQUINOX']                              # Retrieve equinox

        COORDINATES = SkyCoord(float(RA)*u.deg,float(DEC)*u.deg,frame='icrs')
        LOC = EarthLocation.of_site('Las Campanas Observatory')

    elif observatory == 'MPI2.2': # NOTE THAT IT IS NOT GETTING THE AVERAGE DATETIME OF THE EXPOSURE
        OBJECT = header0['OBJECT'].replace('.','_').split('_')[0] # Retrieve the target name
        DATE_OBS = header0['ARCFILE'].split('.')[-3]              # Retrieve the datetime of observation
        RA = header0['RA']                                        # Retrieve Right Ascension in degrees
        DEC = header0['DEC']                                      # Retrieve Declination in degrees
        EQUINOX = header0['EQUINOX']                              # Retrieve equinox

        COORDINATES = SkyCoord(float(RA)*u.deg,float(DEC)*u.deg,frame='icrs')
        LOC = EarthLocation.of_site('La Silla Observatory')

    elif observatory == 'MagellanII': # NOTE THAT IT IS NOT GETTING THE AVERAGE DATETIME OF THE EXPOSURE
        OBJECT = header0['OBJECT'].replace('.','_').split('_')[0] # Retrieve the target name
        DATE_OBS = header0['UT-DATE'] + 'T' + header0['UT-START'] # Retrieve the datetime of observation
        RA = header0['RA-D']                                      # Retrieve Right Ascension in degrees
        DEC = header0['DEC-D']                                    # Retrieve Declination in degrees
        EQUINOX = header0['EQUINOX']                              # Retrieve equinox

        COORDINATES = SkyCoord(float(RA)*u.deg,float(DEC)*u.deg,frame='icrs')
        LOC = EarthLocation.of_site('Las Campanas Observatory')

    rv_correction = COORDINATES.radial_velocity_correction(kind=correction, \
                    obstime=Time(DATE_OBS),location=LOC).to('km/s')

    return rv_correction.value
