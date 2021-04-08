import re
import time
import os.path
import platform
import progressbar as pb

import warnings; warnings.filterwarnings("ignore")

import numpy as np

import pandas as pd # Substitute at some point

import astropy.units as u
from astropy.io import fits
from astropy.table import Table, vstack, hstack
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
Simbad.add_votable_fields('flux(B)','flux(V)','sptype')

def mainpath(path=None):
    '''
    Function to set the main directory (working directory) where the programs,
    lists, tables, etc. are.

    Parameters
    ----------
    path : str, optional
        If 'default' or 'def', it will choose the default path (see below).

    Returns: Selected main directory path.
    '''

    if path in ['def','default']:
        if platform.system() == 'Darwin':
            defmainpath = '/Users/abelink/MEGA/PhD/'
        elif platform.uname().node == 'msi':
            defmainpath = '/home/abelink/MEGA/PhD/'
        elif platform.uname().node == 'dama.dyn.iac.es':
            defmainpath = '/net/nas/proyectos/hots/adeburgos/'

        mainpath = defmainpath

    elif path == None:
        mainpath = input('Working directory path (default is %s) : ' % defmainpath)
        if mainpath == '': mainpath = defmainpath

    else: mainpath = path

    return mainpath

maindir = mainpath('def')

def datapath(path=None):
    '''
    Function to set the directory where the data (fits) are.

    Parameters
    ----------
    path : str, optional
        If 'default' or 'def', it will choose the default path (see below).

    Returns: Selected data directory path.
    '''

    if path in ['def','default']:
        if platform.system() == 'Darwin':
            defdatapath = '/Users/abelink/Documents/DB/'
        elif platform.uname().node == 'msi':
            defdatapath = '/media/abelink/Orange/PhData/DB/'
        elif platform.uname().node == 'dama.dyn.iac.es':
            defdatapath= '/net/nas/proyectos/hots/masblue/obs_iac/spec_opt/IACOB_DB/'

        datapath = defdatapath

    elif path == None:
        datapath = input("Data directory path (default is %s) : " % defdatapath)
        if datapath == '': datapath = defdatapath

    else: datapath = path

    return datapath

datadir = datapath('def')

def search(myfile,path):
    '''
    Function to search a file within a directory.

    Parameters
    ----------
    myfile : str
        Name of the file to search.

    path : str
        Path where to search for the file.

    Returns: Path to the searched file.
    '''

    f_dir = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file == myfile: f_dir = os.path.join(root,file)
    if f_dir == []:
        print('File %s not found.\n' % myfile)
        return None

    return f_dir


def findstar(spectra=None,SNR=None):
    '''
    Function to get the paths of the searched spectra allowing to limitate the
    results by a minimum SNR if provided in the header, or the one with best SNR.

    Parameters
    ----------
    spectra : str, optional
        Enter the input spectra, either name(s) of the star(s), the fits files
        separated by coma, a .txt/.lst file containing the filenames, or '*'
        if you want to select all the fits files inside the working folder.

    SNR : str/int, optional
        If 'best' as input, it finds only the best SNR spectrum for each star.
        If 'bestMF' same as 'best' but prioritizing spectra from HERMES/FEROS.
        If specified, it returns all the spectra above the chosen SNR.

    Returns: Paths to the files found.
    '''

    if spectra == None:
        while spectra == '':
            print('No file/files were selected.\n')
            spectra = input('Enter the input spectra (name(s), *, *.fits, *.txt/.lst): ')

    if spectra == '*': list_spectra = ['*']

    elif '.lst' in spectra or '.txt' in spectra:
        print('\nSearching file in %s... \n' % (maindir))

        list_dir = search(spectra,maindir); list_spectra = []
        with open(list_dir,'r') as spectra: list_spectra = spectra.read().splitlines()
        if len(list_spectra) == 0:
            print('No spectra in list found.\nExiting...'); return None
        list_spectra = [spectrum.split()[0] for spectrum in list_spectra \
                       if not spectrum.startswith('#') and not spectrum == '']

    elif 'fits' in spectra: list_spectra = spectra.split(',')

    else: list_spectra = spectra.replace(' ','').split(',') # This catches 'HDXXXX, HDYYYY'

    dir_spectra = []
    for spectrum in list_spectra:
        match = 0
        for root, dirs, files in os.walk(datadir):
            for file in files:
                if file.startswith('._'): continue
                elif spectrum == '*' and file.endswith('.fits'):
                    dir_spectra.append(os.path.join(root,file))
                    match = 1
                elif file == spectrum:
                    dir_spectra.append(os.path.join(root,file))
                    match = 1
                elif spectrum + '_' in file and file.endswith('.fits'):
                    dir_spectra.append(os.path.join(root,file))
                    match = 1
        if not match == 1: print('File %s not found.\n' % (spectrum))

    if len(dir_spectra) == 0: quit()

    # Spectra selection based on selected SNR.
    if SNR == 'best': dir_spectra = snr(dir_spectra)
    elif SNR == 'bestMF': dir_spectra = snr(dir_spectra,get_MF=True)
    elif type(SNR) == int: dir_spectra = snr(dir_spectra,snrcut=SNR)

    # Order all spectra from a single target by date.
    if len(list_spectra) == 1:
        dir_spectra = sorted(dir_spectra, key = lambda path: \
        path.split('/')[-1].split('_')[1] + path.split('/')[-1].split('_')[2])

    return dir_spectra


def searchlines(line,tol=1):
    line = float(line.replace('?',''))
    linesdb = findlines('synt_lines_OB.lst')
    indexes = [linesdb[0].index(i) for i in linesdb[0] if i > line-tol and i < line+tol]
    print('Nearest lines are:')
    [print(linesdb[0][idx],linesdb[1][idx],linesdb[2][idx]) for idx in indexes]
    line = input('Choose a line from the list: ')

    return line


def findlines(list):
    '''
    Function to extract atomic lines from a list containing their information or
    provide output format for a given wavelenght.

    Parameters
    ----------
    list : str/float
        Enter the list of lines to fit in .txt/.lst, or coma-separated string
        with wavelenghts, or single float/int wavelenght.

    Returns: List of wavelenghts, element names, and loggf.
    '''

    path = maindir+'lists/lines'

    lines = []

    # Single input in float or int format without quotation marks:
    if type(list) == float or type(list) == int:
        lines.append(float(list)); elements = loggf = [None]

    # Lines in a lst/txt file with more information on each line:
    elif '.lst' in list or '.txt' in list:
        list_dir = search(list,path)
        with open(list_dir,'r') as file_lines: lines = file_lines.read().splitlines()
        rows = [line.split(',') for line in lines if not line.startswith('#') and not line == '']
        lines = [float(line[0]) for line in rows]
        try:
            elements = [line[1].strip() for line in rows]
            loggf = [float(line[2].strip()) for line in rows]
        except: elements = loggf = [None]*len(lines)

    # String of lines separated by coma:
    else:
        lines = list.split(','); elements = loggf = [None]*len(lines)

        lines_f = []
        for n in range(len(lines)):
            if '?' in lines[n]: lines_f.append(searchlines(lines[n]))
            else: lines_f.append(lines[n])
        lines = [float(line) for line in lines_f]

    return lines,elements,loggf


def findlist(list):
    '''
    Function to extract items in a txt/lst file.

    Parameters
    ----------
    list : str
        Enter the list of items in .txt/.lst.

    Returns: Items contained in the input list.
    '''

    path = maindir+'lists'

    # To catch wrong int/float inputs:
    if type(list) == float or type(list) == int:
        print('Input cannot be int or float format. \n Exiting...'); return None

    # Lines in a lst/txt file with more information on each line:
    elif '.lst' in list or '.txt' in list:
        list_dir = search(list,path)
        with open(list_dir,'r') as file_list: list = file_list.read().splitlines()

        items = [row.split(',')[0] for row in list if not row.startswith('#') and not row == '']

    return items


def findtable(table,path=None):
    '''
    Function to get the data from a FITS-format table.

    Parameters
    ----------
    table : str
        Enter the fits table containing the data.

    path : str,optional
        Path where to search for the file.

    Returns: Data in table, in table format.
    '''

    if path == None: path = maindir+'tables'

    table_dir = search(table,path)

    if '.fits' in table:
        #try:
        #    with fits.open(table_dir,mode='readonly') as hdu_list:
        #        data = hdu_list[1].data
        #except: data = Table.read(table_dir,format='fits')

        data = Table.read(table_dir,format='fits')

    else: data = Table.read(table_dir,format='ascii',delimiter=' ')

    return data


def snr(spectra,snrcut=None,get_MF=None):
    '''
    Function to provide the spectrum with best signal to noise ratio, or all the
    spectra above a given value of signal to noise ratio, taken from the header.

    Parameters
    ----------
    spectra : list
        List of paths of the spectra to filter.

    snrcut : int, optional
        If established, it returns all the spectra above the chosen SNR.

    get_MF : Boolean, optional
        If True, it returns available spectra from either Mercator or Feros with
        SNR within 15% less than the best SNR spectra taken with FIES.

    Returns: Paths to filtered spectrum/spectra.
    '''

    names_stars = []
    for spectrum in spectra:
        name_star = spectrum.split('/')[-1].split('_')[0]
        if name_star not in names_stars: names_stars.append(name_star)

    best_spectra = []
    for star in names_stars:
        SNR_best = 0; SNR_best_MF = 0
        for spectrum in spectra:
            filename = spectrum.split('/')[-1].split('_')
            name_star = filename[0]
            date = int(filename[1]+filename[2])
            instr = filename[3]
            if star != name_star: continue
            else:
                # Retrieve the key values fron the fits header
                hdu = fits.open(spectrum)# Open the fits image file
                hdu0 = hdu[0]            # Load the header list of primary header
                header0 = hdu0.header    # Read the values of the headers
                SNR = float(header0['I-SNR'])  # Estimated Signal to Noise Ratio

                if snrcut == None:
                    # Date is used for spectra with same SNR choosing the newest one.
                    # Instr is used for when get_MF is enabled.
                    if SNR > SNR_best or (SNR == SNR_best and date > best_spec_date):
                        SNR_best = SNR; best_spec = spectrum
                        best_spec_date = date; best_spec_inst = instr

                    if get_MF == True and instr in ['M','F'] and SNR > SNR_best_MF:
                        SNR_best_MF = SNR; best_spec_MF = spectrum

                elif SNR > int(snrcut): best_spectra.append(spectrum)

        if snrcut == None:
            if get_MF == True and best_spec_inst=='N' and SNR_best-SNR_best_MF<.15*SNR_best:
                 best_spec = best_spec_MF

            best_spectra.append(best_spec)

        elif len(best_spectra) == 0:
            print('No spectra found with SNR higher than %s.' % (snrcut))

    return best_spectra


def gen_fits(list, db, coords=None, limdist=None, spt=None, lc=None, snrcut=None,
    spccode=None, bmag=None, vmag=None, gaia=None, skip=None):
    '''
    Function to generate a FITS table with information about sources coming from
    IACOB/FEROS database, a list of names or coordinates, allowing to limitate
    the results by B/V magnitudes, SpT, LC, distance, and also providing Gaia data.

    Parameters
    ----------
    list : str
        Enter the input list, either name(s)/FITS of the source(s) separated by coma,
        or a .txt/.lst file containing the source names or coordinates.
        (Coordinates must be provided as h:m:s +-d:m:s or d +-d without comas)
        (If "db" is "IACOB", type '*' to select all the availabla FITS)

    db : str
        Enter the input database: IACOB/Simbad

    coords : str, optional
        Enter 'header' to take the coordinates from header. Otherwise it
        takes the coordinates from Simbad, quering the fits filename.

    limdist : list, optional
        Enter the RADEC coordinates [hms] or [deg] of the origin point from
        where to measure the distance and the distance [deg] where to find stars.
        e.g. ['12:30:30.2 +40:20:10.3', 3] or ['35.2368 +57.6622',4.5]

    spt : str, optional
        Enter a desired spectral types to search, separated by coma e.g. 'O,B1'.

    lc : str, optional
        Enter a desired luminosity classes to search, separated by coma e.g. 'I,V'.

    snrcut : int, optional
        If specified, it returns all the spectra above the chosen SNR.

    spccode : str, optional
        If 'y' or 'yes', it will create separate columns with SpT and LC numbers.

    bmag : str, optional
        Enter a desired bmag to cut the results e.g. '<6'.

    vmag : str, optional
        Enter a desired vmag to cut the results e.g. '>8.5'.

    gaia : str, optional
        If 'y' or 'yes', it will create separate columns with Gaia data.

    skip : str, optional
        Enter a coma separated list of targets to exclude in the table.

    Returns: None (Saves the generated table in the */table/ folder)
    '''

    if db == 'IACOB':
        lst_sources_all = findstar(spectra=list,SNR=snrcut)
        lst_sources_f = snr(lst_sources_all)
        #lst_sources_f.sort()
        type_list = 'names'

    elif db == 'Simbad':
        lst_sources_f = findlist(list)
        type_list = input('The list contains names or coordinates? [names/coords]: ')
        if not type_list in ['names','coords']:
             print('Input answer is not valid. Exiting... \n'); return None

    else: print('Database not recognised. Exiting... \n'); return None

    # Tuning the format of the input variables
    if limdist != None:
        RADEC = limdist[0]
        dist = float(limdist[1])

    if spt != None:
        spt = re.split(' |,',spt)

    if lc != None:
        lc = re.split(' |,',lc)


    '''=============================== Queries =============================='''
    columns = ['BPmag','e_BPmag','+Gmag','e_Gmag','RPmag','e_RPmag',\
               'pmRA','e_pmRA','pmDE','e_pmDE','Plx','e_Plx']

    if gaia in ['y','yes']:
        ruwe = input('Calculate RUWE? [y/n]: ')

        if ruwe == 'y':
            columns.extend(['astrometric_n_good_obs_al','astrometric_chi2_al'])
            table_u0 = pd.read_csv('/home/abelink/PhD/tables/Gaia/table_u0_g_col.txt',
                       delimiter=',',header=1,names=['g_mag','bp_rp','u0'])

        offset = input('Apply +0.03 mas offset to parallax? [y/n]: ')

    v = Vizier(columns=columns); v.ROW_LIMIT = 1


    '''============================ Progress Bar ============================'''
    bar = pb.ProgressBar(maxval=len(lst_sources_f),
                         widgets=[pb.Bar('=','[',']'),' ',pb.Percentage()])
    bar.start()


    '''=============================== Sources =============================='''
    delete = []; first = True
    for source,i in zip(lst_sources_f,range(len(lst_sources_f))):


        '''=========== Retrieve the key values fron the fits header ========='''
        OBJRA = OBJDEC = None
        if db == 'IACOB':
            hdu = fits.open(source)  # Open the fits image file
            hdu0 = hdu.verify('fix') # Fix header keywords
            hdu0 = hdu[0]            # Load the header list of primary header
            header0 = hdu0.header    # Read the values of the headers
            if '_M_' in source:
                try: OBJRA = header0['OBJ_RA']; OBJDEC = header0['OBJ_DEC']
                except: None # Only new fits include it
            elif '_N_' in source:
                OBJRA = header0['OBJRA']*360/24; OBJDEC = header0['OBJDEC']
            elif '_F_' in source:
                OBJRA = header0['RA']; OBJDEC = header0['DEC']

            SNR_best = float(header0['I-SNR'])

            source = source.split('/')[-1].split('_')[0]


        '''========================= Skip bad sources ======================='''
        if db == 'IACOB':
            if any(bad in source[:3] for bad in ['DO2']): continue

        if skip != None and any(bad in source for bad in skip.split(',')): continue


        '''============= Simbad query by object name/coordinates ============'''
        if type_list == 'names':

            simbad = SB(source,OBJRA,OBJDEC)

            if simbad == None: continue

            source_f = source

        elif type_list == 'coords':

            if ':' in source: c = SkyCoord(source,unit=(u.hourangle,u.deg))
            else: c = SkyCoord(source,unit=u.deg)
            try: simbad = Simbad.query_region(c,radius=2*u.arcsec)
            except: time.sleep(3); simbad = Simbad.query_region(c,radius=2*u.arcsec)

            try:
                if len(simbad) > 1: simbad = Table(simbad[0])
            except: print('Source %s not found in Simbad' % (source)); continue

            source_f = simbad['MAIN_ID'][0].decode()


        '''======================= Get the coordinates ======================'''
        if coords == 'header':
            RADEC_0 = SkyCoord(ra=header0['RA'],dec=header0['DEC'],unit=(u.deg))

        else:
            RADEC_0 = SkyCoord(ra=simbad['RA'],dec=simbad['DEC'],unit=(u.hourangle,u.deg))[0]

        RAdeg = RADEC_0.ra.deg; DECdeg = RADEC_0.dec.deg

        RADEC_0 = re.sub('[h,d,m]',':',RADEC_0.to_string('hmsdms')).replace('s','')
        RAhms = RADEC_0.split()[0]; DECdms = RADEC_0.split()[1]


        '''======================= Limit by distance ========================'''
        if limdist != None:
            c1 = SkyCoord(RADEC_0,unit=(u.hourangle,u.deg))
            if any(i in RADEC for i in [':','h']):
                  c2 = SkyCoord(RADEC,unit=(u.hourangle,u.deg))
            else: c2 = SkyCoord(RADEC,unit=u.deg)

            if c1.separation(c2).deg > dist: delete.append(source); continue


        '''===================== Get the spectral class ====================='''
        if db == 'IACOB':
            SpC_0 = header0['I-SPC']
            try: SpC_ref = header0['I-SPCREF']
            except: SpC_ref = '-'
            if (not type(SpC_0) == str or SpC_0.strip() == '-'):
                SpC_0 = simbad['SP_TYPE'][0]; SpC_ref = 'SIMBAD'
        else:
            SpC_0 = simbad['SP_TYPE'][0]; SpC_ref = 'SIMBAD'


        '''======================= Limit by SpT or LC ======================='''
        spt_0 = []; lc_0 = []
        if spt != None or lc != None:
            if 'I' in SpC_0:
                spt_0 = SpC_0[:SpC_0.index('I')]
                lc_0 = SpC_0[SpC_0.index('I'):]
            elif 'V' in SpC_0:
                spt_0 = SpC_0[:SpC_0.index('V')]
                lc_0 = SpC_0[SpC_0.index('V'):]
            elif len(re.split(':|pe',SpC_0.strip())[0]) <= 4:
                spt_0 = SpC_0.strip()

            if spt != None:
                if any(j in spt_0 for j in spt) == False:
                    delete.append(source); continue

            match = False
            if lc != None and lc_0 != []:
                if 'IV' in lc and 'IV' in lc_0: match = True
                lc_0 = lc_0.replace('IV','')
                if 'V' in lc and 'V' in lc_0: match = True
                if 'III' in lc and 'III' in lc_0: match = True
                lc_0 = lc_0.replace('III','')
                if 'II' in lc and 'II' in lc_0: match = True
                lc_0 = lc_0.replace('II','')
                if 'I' in lc and 'I' in lc_0: match = True
                if match == False: delete.append(source); continue


        '''=================== Get the spectral class code =================='''
        if spccode in ['y','yes'] and SpC_0 != '':
            spc_c = SpC_0.split('+')[0].split('/')[0].replace(':','')

            if (spc_c in ['~']) or (spc_c.startswith(('WC','WN','WR','C')) == True):
                spt_c = float('NaN'); lc_c = float('NaN')
            else:
                spt_code = {'O': 1, 'B': 2, 'A': 3, 'F': 4, 'G': 5, 'K': 6, 'M': 7}
                total_spt_c = []

                spt_c = re.findall('[O,B,A,F,G,K,M]',spc_c); num = len(spt_c)

                for spt_c_i in spt_c: total_spt_c.append(spt_code[spt_c_i])

                spt_c = re.findall('[0-9.]+',spc_c)

                for spt_c_i in spt_c: total_spt_c.append(float(spt_c_i)/10)

                try:
                    # This is for B1-2 in which there is B(1) and 1-2(2)
                    if num < len(total_spt_c)-num:
                        spt_c = (sum(total_spt_c[:num])/num + \
                                 sum(total_spt_c[num:])/len(total_spt_c[num:]))
                    else: spt_c = sum(total_spt_c)/num
                except: print(source,total_spt_c,num)

                lc_code = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5}
                total_lc_c = []

                lc_c = re.findall('[I,V]+',spc_c)

                for lc_c_i in lc_c: total_lc_c.append(lc_code[lc_c_i])

                lc_c = np.asarray(total_lc_c).mean()

        elif spccode in ['y','yes'] and SpC_0 == '': spt_c = lc_c = np.nan


        '''======================= Limit by magnitude ======================='''
        bmag_0 = simbad['FLUX_B']
        vmag_0 = simbad['FLUX_V']

        if bmag != None and str(bmag_0) != '--':
            bmag_0 = round(bmag_0, 2)
            if bmag[0] == '<':
                if bmag_0 > float(bmag[1:]): delete.append(source); continue
            elif bmag[0] == '>':
                if bmag_0 < float(bmag[1:]): delete.append(source); continue

        if vmag != None and str(vmag_0) != '--':
            vmag_0 = round(vmag_0, 2)
            if vmag[0] == '<':
                if vmag_0 > float(vmag[1:]): delete.append(source); continue
            elif vmag[0] == '>':
                if vmag_0 < float(vmag[1:]): delete.append(source); continue


        '''======================= Get Gaia DR2 data ========================'''
        if gaia in ['y','yes']:
            gFlag = True; catalog = "I/345/gaia2"

            if vmag_0 != '' and float(vmag_0) < 5:
                print('\nWARNING: '+str(source_f)+' is too bright (Vmag = %r)' %vmag_0)

            try:
                gaiaq = v.query_object(source_f,catalog=catalog,radius=2*u.arcsec)[0]

                # Correct Gaia photometry:
                if str(gaiaq['Gmag'][0]) != '--' or str(gaiaq['e_Gmag'][0]) != '--':
                    if gaiaq['e_Gmag'] < 8e-3:  gaiaq['e_Gmag'] = 8e-3
                    if (6 <= gaiaq['Gmag'] <= 16):
                        gaiaq['Gmag'] = gaiaq['Gmag'] - 0.0032*(gaiaq['Gmag']-6)
                    if gaiaq['Gmag'] < 6:
                        gaiaq['Gmag'] = gaiaq['Gmag'] + 0.0271*(6-gaiaq['Gmag'])
                    gaiaq['Gmag'].unit = u.mag
                else: gFlag = False; print(str(source_f) + ' has missing Gaia G photometry.')

                if str(gaiaq['e_BPmag'][0]) != '--':
                    if gaiaq['e_BPmag'] < 9e-3:  gaiaq['e_BPmag'] = 9e-3
                else: gFlag = False; print(str(source_f) + ' has missing Gaia B photometry.')
                if str(gaiaq['e_RPmag'][0]) != '--':
                    if gaiaq['e_RPmag'] < 10e-3: gaiaq['e_RPmag'] = 10e-3
                else: gFlag = False; print(str(source_f) + ' has missing Gaia R photometry.')

                # Correct Gaia astrometry:
                try:
                    if offset == 'y': gaiaq['Plx'] = gaiaq['Plx'] + 0.03
                except: print(str(source_f) + ' has missing Gaia parallax.')

            except:
                gaiaq = []; print('\n' + str(source_f) + ' could not be queried in Gaia.')
                #catalog = "/I/311/hip2"


        '''========================= Calculate RUWE ========================='''
        if gaia in ['y','yes'] and ruwe == 'y' and not gaiaq == []:

            if gFlag == True:
                UWE = np.sqrt(float(gaiaq['chi2AL'])/(float(gaiaq['NgAL']) - 5))

                diff = abs(gaiaq['Gmag'][0] - table_u0['g_mag']) + \
                       abs(gaiaq['BPmag'][0] - gaiaq['RPmag'][0] - table_u0['bp_rp'])
                diff = diff.tolist()

                RUWE = round(UWE/table_u0['u0'][diff.index(min(diff))],4)

            gaiaq.remove_columns(['NgAL','chi2AL'])


        '''============= Counting FIES / HERMES / FEROS spectra ============='''
        if db == 'IACOB':
            FIES = 0; HERMES = 0; FEROS = 0
            for source_ in lst_sources_all:
                if source_f + '_' in source_ and source_.endswith('.fits'):
                    if '_N_' in source_: FIES = FIES + 1
                    elif '_M_' in source_: HERMES = HERMES + 1
                    elif '_F_' in source_: FEROS = FEROS + 1


        '''=========================== Export row ==========================='''
        output = Table([[source_f],[RAhms],[DECdms],[RAdeg],[DECdeg],[SpC_0],[SpC_ref]],\
         names = ('Name','RA_J2000','DEC_J2000','RAdeg_J2000','DECdeg_J2000','SpC','SpC_ref'))

        #output['RA_J2000'].unit = u.hourangle; output['DEC_J2000'].unit = u.deg
        output['RAdeg_J2000'].unit = u.deg; output['DECdeg_J2000'].unit = u.deg

        if spccode in ['y','yes']:
            spc_code = Table([[spt_c],[lc_c]],names=('SpT_code','LC_code'))
            output = hstack([output,spc_code])

        magnitudes = Table([[bmag_0],[vmag_0]],names=('mag_B','mag_V'))
        magnitudes['mag_B'].unit = magnitudes['mag_V'].unit = u.mag
        output = hstack([output,magnitudes])

        if gaia in ['y','yes'] and not gaiaq == []:
            output = hstack([output,gaiaq])

            if ruwe == 'y' and gFlag == True:
                ruwe_table = Table([[RUWE]],names=['RUWE'])
                output = hstack([output,ruwe_table])

        if db == 'IACOB':
            nspect = Table([[FIES],[HERMES],[FEROS]],names=('FIES','HERMES','FEROS'))
            output = hstack([output,nspect])

            SNR_best = Table([[SNR_best]],names=('SNR_best',))
            output = hstack([output,SNR_best])

        if first == True: table = output; first = False
        else: table = vstack([table,output])

        bar.update(i)#; print(' ')
        time.sleep(0.1)

    bar.finish()


    '''============================ Export table ============================'''
    try:
        hdu_f = fits.BinTableHDU(data=table.filled(np.nan))
        hdu_f.writeto(maindir+'tables/tablestars.fits',overwrite=True)
    except: print('Table is empty, no sources were found.'); return None

    lst_sources_f = [source for source in lst_sources_f if not source in delete]

    return None


def SB(name=None,ra=None,dec=None):
    '''
    Function to query an object in Simbad database.

    Parameters
    ----------
    name : str, optional
        Enter the name of the source to query.

    ra : float, optional
        Enter the right ascension of the source, in degrees.

    dec : float, optional
        Enter the declination of the source, in degrees.

    Returns: Queried object in Table format.
    '''

    if name != None:
        try: simbad = Simbad.query_object(name)
        except: time.sleep(3); simbad = Simbad.query_object(name)
    else: name = 'empty'; simbad = None

    while type(simbad) == type(None):
        print('Provide alternative name for %s in Simbad.' % name)
        print('In some cases try replacing "HD" by "HD ".')
        print('Type "sky" to query 5" around the central coordinates (if given).')
        print('Hit return to skip this source.')
        check = input('Name: ')

        if check == '': print('Skipping source: %s\n' % name); break

        elif check == 'sky' and ra != None and dec != None:
            simbad = Simbad.query_region(SkyCoord(ra,dec,unit='deg'),radius='5s')
            if type(simbad) == type(None): print('No objects found.')

        else: simbad = Simbad.query_object(check)

        if type(simbad) != type(None) and len(simbad) > 1:
            print('More than one Simbad result, choosing the brigtest source...')
            simbad.sort('FLUX_V'); simbad = Table(simbad[0])

    return simbad


def zp_edr3(ra,dec,radius=1):
    '''
    Function to obtain avaliable parallax zero-point offsets.
    Input is a list of coordinates in degrees which are used for the Gaia query.

    Parameters
    ----------
    ra : float,list
        Enter the right ascension of the source, in degrees.

    dec : float,list
        Enter the declination of the source, in degrees.

    radius : int,float
        Enter the search radius in arcsec for the Gaia query. Default is 1.

    Returns: Table with the queried sources including the zero-point offset.

    '''

    from astroquery.gaia import Gaia

    radius = radius/60/60

    first = True
    for ra,dec in zip(ra,dec):
        job = Gaia.launch_job("select TOP 1 * FROM gaiaedr3.gaia_source "
                    "WHERE 1=CONTAINS(POINT('ICRS',ra,dec), "
                    "CIRCLE('ICRS',%f,%f,%f))" % (ra,dec,radius)).get_results()
        if len(job) == 0: continue
        elif len(job) > 1: job.sort('phot_g_mean_mag')

        if first == True: table = job; first = False
        else: table.add_row(job[0])

    table = table[table['astrometric_params_solved']>3]
    zpvals = zpt.get_zpt(table['phot_g_mean_mag'],table['nu_eff_used_in_astrometry'],\
           table['pseudocolour'],table['ecl_lat'],table['astrometric_params_solved'])
    table.add_column(zpvals,name='zp_offset')

    table.remove_column('designation')

    return table


def checknames():
    '''
    Function to detect errors in the filenames/headers.

    Parameters
    ----------
    (empty)

    Returns: None (Generate output files with the found issues)
    '''

    dir_spectra = findstar()

    errorslist = open(maindir+'lists/ErrorNames.txt', 'w')

    bar = pb.ProgressBar(maxval=len(dir_spectra),
                        widgets=[pb.Bar('=','[',']'),' ',pb.Percentage()])
    bar.start()

    date_0_long = []; errors = 0; i = 0
    for spectrum in dir_spectra:

        bar.update(i + dir_spectra.index(spectrum))
        time.sleep(0.1)

        # Retrieve the key values fron the fits header
        hdu = fits.open(spectrum)# Open the fits image file
        hdu0 = hdu[0]            # Load the header list of primary header
        header0 = hdu0.header    # Read the values of the headers

        filename = spectrum.split('/')[-1]
        namestar = spectrum.split('/')[-1].split('_')[0]

        simbad = SB(namestar)

        name_0 = header0['OBJECT'].strip().replace(' ','')

        date_filename = spectrum.split('/')[-1].split('_')[1] # date from the filename
        if '_N_' in filename or '_M_' in filename:
            date_0 = header0['DATE-AVG'] #; print date_0
            date_0_short = str(date_0[0:4]) + str(date_0[5:7]) + str(date_0[8:10])
        elif '_F_' in filename:
            date_0 = header0['ARCFILE']
            date_0_short = str(date_0[6:10]) + str(date_0[11:13]) + str(date_0[14:16])

        RA_0 = header0['RA']   # 'RA'
        DEC_0 = header0['DEC'] # 'DEC'
        RADEC_0 = str(RA_0) + ' ' + str(DEC_0)

        try:
            RADEC_sim = str(simbad['RA'][0]).replace(' ',':') + ' ' + \
                        str(simbad['DEC'][0]).replace(' ',':')
        except:
            errorslist.write('Problem quering ' + filename + ' in Simbad' + '\n')
            errors = errors + 1; continue

        c1 = SkyCoord(RADEC_0,unit=u.deg); c2 = SkyCoord(RADEC,unit=(u.hourangle,u.deg))
        difcoord = round(c1.separation(c2).arcsec,3)

        if difcoord > 90:
            errorslist.write('Distance from Simbad query is ' + str(difcoord) + \
                            ' arcsec for spectrum ' + filename + '\n')
            errors = errors + 1

        # Catch files with wrong object name comparing filename and header
        if namestar != name_0.upper():
            errorslist.write(filename + ' has a wrong object file name!' + '\n')
            errors = errors + 1

        # Catch two+ files with same date in header but different filenames
        # e.g. HD111111 and HDE111111
        if date_0 not in date_0_long:
            date_0_long.append(date_0)
        else:
            errorslist.write(filename + \
            ' is a duplicated file with same full date but different filename!' + '\n')
            errors = errors + 1

        # Catch wrong dates in the filenames when comparing with the date in the header
        if not date_filename == date_0_short:
            errorslist.write(filename + ' is a repeated file with wrong name!' + '\n')
            errors = errors + 1

    bar.finish()

    errorslist.close()

    if errors == 0: print('\nEverything is OK...')
    else: print('\nSome errors found, check ErrorNames.txt')
