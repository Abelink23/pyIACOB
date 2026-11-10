# Python Standard Library packages:
import re
import time
from pathlib import Path

# Other main packages
import numpy as np
import progressbar as pb

# Astro-packages
import astropy.units as u
from astropy.io import ascii, fits
from astropy.table import Table, join, vstack, hstack
from astropy.coordinates import SkyCoord

from color_msg import msg

msg = msg()
msg.bold_under('c', '\nWelcome to the pyIACOB package! - v1.24')
msg.c('developed by Abel de Burgos et al. (2026/11/09)\n')


def load_default_paths(dir_paths_file='paths.txt'):
    '''
    Function to load the default paths for the pyIACOB package from a text file.

    Parameters
    ----------
    dir_paths_file : str | Path, optional
        Name of the text file containing the paths. It must be in the same path
        where the package is being run. If not found, it will ask for the full
        path to the file.

    Returns
    -------
    tuple[Path, Path, Path, Path, Path, Path]
        Paths for the main directories (main, data, ib, maui, tracks, tess).
    '''
    paths_file = Path(dir_paths_file)
    while not paths_file.is_file():
        msg.error(f"File '{paths_file}' not found...")
        paths_file = Path(input("Please provide the full path to the file: ").strip())

    with paths_file.open('r') as f:
        dirs = {
            key.strip(): Path(val.strip()).expanduser()
            for line in f
            if line.strip() and not line.startswith('#') and '=' in line
            for key, val in [line.split('=', 1)]
            }

    return dirs['main'], dirs['data'], dirs['maui']


maindir, datadir, mauidir = load_default_paths('paths.txt')


def search(myfile, path, silent_search=False):
    '''
    Function to search a file within a directory.

    Parameters
    ----------
    myfile : str
        Name of the file to search.
    path : str | Path
        Path where to search for the file.
    silent_search : bool, optional
        If True, suppresses warning and error messages.

    Returns
    -------
    Path | None
        Path to the searched file, or None if not found.
    '''

    search_path = Path(path).expanduser()
    if not search_path.is_dir():
        if not silent_search:
            msg.error(f"Directory '{search_path}' not found.\n")
        return None

    matches = [p for p in search_path.rglob(str(myfile)) if p.is_file() and not p.name.startswith('._')]

    if not matches:
        if not silent_search:
            msg.error(f"File '{myfile}' not found.\n")
        return None

    if len(matches) > 1 and not silent_search:
        msg.warn(f"More than one file found, selecting the first one: {matches[0]}")

    return matches[0]


def findstar(spectra=None, snr=0, path=None):
    '''
    Function to get the paths of the searched spectra allowing to limitate the
    results by a minimum SNR if provided in the header, or the one with best SNR.

    Parameters
    ----------
    spectra : str, optional
        Enter the input spectra, either name(s) of the star(s), the fits files
        separated by coma, a .txt/.lst file containing the (file)names, or '*'
        if you want to select all the fits files inside the working folder.

    snr : str/decimal, optional
        If 'best' as input, it finds only the best SNR spectrum for each star.
        If 'bestHF' same as 'best' but prioritizing spectra from HERMES/FEROS.
        If specified, it returns all the spectra above the chosen SNR.
        If not specified it will search for all the available spectra with SNR > 0.

    path : str, optional
        Path where to search for the spectra. If not specified, it will search
        in the default data directory.

    Returns
    -------
    list[Path]
        Paths to the files found.
    '''

    if spectra is None:
        while not spectra:
            print('No file/files were selected.\n')
            spectra = input('Enter the input spectra (name(s), *, *.fits, *.txt/.lst): ').strip()

    if spectra == '*':
        list_spectra = ['*']
    elif '.lst' in spectra or '.txt' in spectra:
        print(f"\nSearching file in {maindir}... \n")
        list_file = search(spectra, maindir)
        if list_file is None:
            return []

        with open(list_file, 'r') as f:
            raw_lines = f.read().splitlines()

        if not raw_lines:
            msg.error('No spectra in list found.\nExiting...')
            return []

        list_spectra = [line.split()[0] for line in raw_lines if not line.startswith('#') and line.strip()]
    elif 'fits' in spectra:
        list_spectra = [s.strip() for s in spectra.split(',')]
    else:
        list_spectra = spectra.replace(' ', '').split(',')

    search_path = Path(path if path else datadir).expanduser()
    all_files = [p for p in search_path.rglob('*') if p.is_file() and not p.name.startswith('._')]

    dir_spectra = []
    for spectrum in list_spectra:
        match = False
        for p in all_files:
            file_name = p.name
            if spectrum == '*' and p.suffix == '.fits':
                dir_spectra.append(p)
                match = True
            elif file_name == spectrum:
                dir_spectra.append(p)
                match = True
            elif f"{spectrum}_" in file_name and p.suffix == '.fits':
                dir_spectra.append(p)
                match = True

        if not match:
            msg.error(f"File/source {spectrum} not found.\n")

    if not dir_spectra:
        return []

    # Spectra selection based on selected SNR
    if any('ascii' in str(s) for s in list_spectra) and snr is not None:
        msg.warn('SNR for ascii files not available, ignoring SNR keyword.')
    else:
        if snr == 'best' or snr is None:
            dir_spectra = spec_snr(dir_spectra)
        elif snr == 'bestHF':
            dir_spectra = spec_snr(dir_spectra, get_HF=True)
        elif isinstance(snr, (int, float)) or (isinstance(snr, str) and snr.isnumeric()):
            dir_spectra = spec_snr(dir_spectra, snrcut=float(snr))
        else:
            msg.error('SNR keyword not provided properly.\nExiting...')
            return []

    # Order all spectra from a single target by date
    if len(list_spectra) == 1:
        dir_spectra = sorted(dir_spectra, key=lambda p: Path(p).name.split('_')[1] + Path(p).name.split('_')[2])

    return dir_spectra


def searchlines(line, tol=1):
    '''
    IN DEVELOPMENT - this is to easily find spectral lines in a future included list of
    spectral lines developed within the IACOB team.

    Parameters
    ----------
    line : float
        Approximate wavelength of the line to search

    tol : int/float, optional
        Tolerance in the search of the line

    Returns
    -------
    List of lines found for the input line.
    '''

    line = float(line.replace('?', ''))
    linesdb = findlines('synt_lines_OB.lst')
    indexes = [linesdb[0].index(i) for i in linesdb[0] if line - tol < i < line + tol]
    print('Nearest lines are:')
    for idx in indexes:
        print(linesdb[0][idx], linesdb[1][idx], linesdb[2][idx])
    return input('Choose a line from the list: ')


def findlines(list_lines):
    '''
    Function to extract atomic lines from a list containing their information or
    provide output format for a given wavelength.

    Parameters
    ----------
    list_lines : str/float
        Enter the list of lines to fit in .txt/.lst, or coma-separated string
        with wavelengths, or single float/int wavelength.

    Returns
    -------
    List of wavelengths, element names, and loggf.
    '''

    path = maindir / 'lists' / 'lines'
    elements = loggf = []

    if isinstance(list_lines, (float, int)):
        return [float(list_lines)], [None], [None]

    if isinstance(list_lines, str) and ('.lst' in list_lines or '.txt' in list_lines):
        list_file = search(list_lines, path)
        if list_file is None:
            return [], [], []

        with open(list_file, 'r') as f:
            raw_lines = f.read().splitlines()

        rows = [line.split(',') for line in raw_lines if not line.startswith('#') and line.strip()]
        lines = [float(row[0]) for row in rows]

        try:
            elements = [row[1].strip() for row in rows]
        except IndexError:
            elements = [None] * len(lines)

        try:
            loggf = [float(row[2].strip()) for row in rows]
        except IndexError:
            loggf = [None] * len(lines)

        return lines, elements, loggf

    # Comma-separated string
    raw_entries = list_lines.split(',')
    elements = loggf = [None] * len(raw_entries)
    lines_f = [searchlines(entry) if '?' in entry else entry for entry in raw_entries]
    lines = [float(line) for line in lines_f]

    return lines, elements, loggf


def findlist(list_name, path=None):
    '''
    Function to extract items in a txt/lst file.

    Parameters
    ----------
    list_name : str
        Enter the list of items in .txt/.lst.

    path : str, optional
        Path where to search for the file.

    Returns
    -------
    Items contained in the input list.
    '''

    search_path = Path(path) if path is not None else maindir / 'lists'

    if isinstance(list_name, (float, int)):
        msg.error('Input cannot be int or float format. \n Exiting...')
        return None

    if isinstance(list_name, str) and ('.lst' in list_name or '.txt' in list_name):
        list_file = search(list_name, search_path)
        if list_file is None:
            return []

        with open(list_file, 'r') as f:
            raw_lines = f.read().splitlines()

        return [row.split(',')[0].strip() for row in raw_lines if not row.startswith('#') and row.strip()]

    return list_name.replace(' ', '').split(',')


def findtable(table, path=None, format=None, delimiter=',', header_start=0,
    fits_strip_end=True, fix_missing=False, guess=False, quicklook=False,
    silent_search=False):
    '''
    Function to get the data from a table with different formats.
    See https://docs.astropy.org/en/stable/io/ascii/read.html for more details
    of 'format', 'delimiter' and 'header_start' keywords.

    Parameters
    ----------
    table : str
        Enter the fits table containing the data.

    path : str, optional
        Path where to search for the file.

    format : str, optional
        Format of the input table. See ascii.read() and Table.read() documentation.
        Default is None.

    delimiter : str, optional
        The string used to separate values. Default is a comma ','.

    header_start : int, optional
        Index of the row which contains the header. Default is 0.

    fits_strip_end : boolean, optional
        If 'True' it strips all strings within the data of a fits table.

    fix_missing : boolean, optional
        If 'True' it will replace missing strings by '' abd missing floats by float 'nan'.

    guess : boolean, optional
        If 'True' it will try to guess the format of the table.

    quicklook : boolean, optional
        If 'True' it will show the table in the web browser. Default is 'False'.

    silent_search : boolean, optional
        If 'True' it will not print any message if the file is not found of more
        than one file is found. Default is 'False'.

    Note: If the table has a commented header, use 'commented_header' as format.

    Returns
    -------
    Data in table, in table format.
    '''

    search_path = Path(path) if path is not None else maindir / 'tables'
    table_file = search(table, search_path, silent_search=silent_search)

    if table_file is None:
        return None

    if '.fits' in table_file.suffix or format == 'fits':
        data = Table.read(table_file, format='fits')
        if fits_strip_end:
            to_strip = [col for col in data.colnames if data.dtype[col].type in (np.str_, np.bytes_)]
            for col in to_strip:
                data[col] = [val.strip() if not np.ma.is_masked(val) else val for val in data[col]]
    else:
        data = ascii.read(table_file, header_start=header_start, format=format, delimiter=delimiter, guess=guess)

    if fix_missing:
        for col in [c for c in data.colnames if data.dtype[c].type in (np.str_, np.bytes_)]:
            data[col].fill_value = ''
        for col in [c for c in data.colnames if data.dtype[c].type == np.float64]:
            data[col].fill_value = np.nan
        data = data.filled()

    if quicklook:
        data.show_in_browser(jsviewer=True)

    return data


def xmatch_table(table1, table2, match_col=None, output_cols=None,
                 output_name='xmatch_table', format='fits'):
    '''
    Function to generate xmatched tables with selected output columns and format.

    Parameters
    ----------
    table1 : str
        Name of the main table.

    table2 : str
        Name of the second table to cross-match with the first.

    match_col : str, optional
        Name of the column to be used as anchor for the cross-match.
        If empty, the user is asked among the colums that are in both tables.

    output_cols : str, optional
        Coma separated string with the list of column names to be included in
        the output file. If empty, all the columns are selected.

    output_name : str, optional
        Name of the output file. Default is 'xmatch_table'.
        Note: Do not include the extension.

    format : str, optional
        Enter the output format for the table: 'fits' (default), 'ascii' or 'csv'.

    Returns
    -------
    Nothing but the output table with cross-match is generated in the
    same path where the two input tables are located.
    '''

    t1 = findtable(table1)
    t2 = findtable(table2)
    if t1 is None or t2 is None:
        return None

    # Cross-match column:
    same_cols = [col for col in t1.colnames if col in t2.colnames]
    if match_col is None:
        if not same_cols:
            print('No columns found with the same name in both tables. Exiting... \n')
            return None
        msg.info(f"Columns {', '.join(same_cols)} with the same name found in both tables.")
        match_col = input('Select one of the above column names for the cross-match: ').strip()


    # Remove extra columns:
    if len(same_cols) > 1:
        msg.info('More than one column has the same name.')
        msg.info('Output table will contain *_1, *_2 for them.')
        rm = input('Remove repeated columns from second table? [y/n]: ').strip().lower()
        if rm == 'y':
            t2.remove_columns([c for c in t2.colnames if c in t1.colnames and c != match_col])

    # Cross-match:
    if t1[match_col].dtype.char == 'S':
        t1[match_col] = [val.strip() for val in t1[match_col]]
        t2[match_col] = [val.strip() for val in t2[match_col]]

    matched = join(t1, t2, keys=match_col)
    if len(matched) == 0:
        msg.error('No matches for the joined tables. Exiting... \n')
        return None

    # Output columns:
    if output_cols is not None:
        cols_to_keep = [c.strip() for c in output_cols.split(',')]
        matched.remove_columns([c for c in matched.colnames if c not in cols_to_keep])

    # Saving the file:
    out_dir = maindir / 'tables'
    out_dir.mkdir(parents=True, exist_ok=True)
    full_path = out_dir / f"{output_name}.{format}"

    out_format = format + '.fixed_width_two_line' if format == 'ascii' else format
    matched.write(full_path, format=out_format, overwrite=True)


def spec_snr(spectra, snrcut=None, get_HF=None):
    '''
    Function to provide the spectrum with best signal to noise ratio, or all the
    spectra above a given value of signal to noise ratio, taken from the header.

    Parameters
    ----------
    spectra : list
        List of paths of the spectra to filter.

    snrcut : int, optional
        If established, it returns all the spectra above the chosen SNR.

    get_HF : Boolean, optional
        If True, it returns available spectra from either Mercator or Feros with
        SNR within 15% less than the best SNR spectra taken with FIES.

    Returns
    -------
    Paths to filtered spectrum/spectra
    '''

    names_stars = []
    for spectrum in spectra:
        # Check that the spectrum is a valid file
        if not Path(spectrum).is_file():
            msg.error(f"File '{spectrum}' not found. Skipping...")
            continue

        id_star = Path(spectrum).name.split('_')[0]
        if id_star not in names_stars:
            names_stars.append(id_star)

    best_spectra = []
    for star in names_stars:
        SNR_best = 0
        SNR_best_HF = 0
        best_spec = None
        best_spec_date = 0
        best_spec_inst = ''
        best_spec_MF = None

        for spectrum in spectra:
            filename = Path(spectrum).name.split('_')
            id_star = filename[0]
            date = int(filename[1] + filename[2])
            instr = filename[3]

            if star != id_star:
                continue

            with fits.open(Path(spectrum)) as hdu:
                # Read the values of the headers
                header = hdu[0].header
                SNR = float(header['I-SNR'])

            if snrcut is None:
                # Date is used for spectra with same SNR choosing the newest one.
                # Instr is used for when get_HF is enabled.
                if SNR > SNR_best or (SNR == SNR_best and SNR != 0 and date > best_spec_date):
                    SNR_best = SNR
                    best_spec = spectrum
                    best_spec_date = date
                    best_spec_inst = instr

                if get_HF is True and instr in ('M', 'F') and SNR > SNR_best_HF:
                    SNR_best_HF = SNR
                    best_spec_MF = spectrum
            elif SNR > float(snrcut):
                best_spectra.append(spectrum)

        if snrcut is None and best_spec is not None:
            if get_HF is True and best_spec_inst == 'N' and (
                0.80 * SNR_best < SNR_best_HF or (0.65 * SNR_best < SNR_best_HF and SNR_best_HF > 100)
            ):
                best_spec = best_spec_MF
            best_spectra.append(best_spec)

    if snrcut is not None and not best_spectra:
        msg.error(f"No spectra found with SNR higher than {snrcut}.")

    return best_spectra


def table_db(list_name, db, coords=None, limdist=None, lim_lb=None, spt=None, lc=None,
             snrcut=None, spccode=False, bmag=None, vmag=None, gaia=False, radius=1, skip=None):
    '''
    Function to generate a FITS table with information about sources coming from
    IACOB/FEROS database, a list of names or coordinates, allowing to limit
    the results by B/V magnitudes, SpT, LC, distance, and also providing Gaia data.

    Parameters
    ----------
    list_name : str
        Enter the input list, either name(s)/FITS of the source(s) separated by coma,
        or a .txt/.lst file containing the source names or coordinates.
        (Coordinates must be provided as h:m:s +-d:m:s or d +-d without comas)
        (If "db" is "IACOB", type '*' to select all the available FITS)

    db : str
        Enter the input database: IACOB/Simbad. If list contains coordinates,
        Simbad option must be used.

    coords : str, optional
        Enter 'header' to take the coordinates from header. Otherwise it
        takes the coordinates from Simbad, quering the fits filename.

    limdist : list, optional
        Enter the RADEC coordinates [hms] or [deg] of the origin point from
        where to measure the distance and the distance [deg] where to find stars.
        e.g. ['12:30:30.2 +40:20:10.3', 3] or ['35.2368 +57.6622',4.5]

    lim_lb : list, optional
        Enter the galactic l(min,max), b(min,max) coordinates in degrees if you want to
        filter the sources which lie outside the input boundaries. e.g. [100,150,-15,15]

    spt : str, optional
        Enter a desired spectral types to search, separated by coma e.g. 'O,B1'.

    lc : str, optional
        Enter a desired luminosity classes to search, separated by coma e.g. 'I,V'.

    snrcut : int, optional
        If specified, it returns all the spectra above the chosen SNR.

    spccode : boolean, optional
        If True, it will create separate columns with SpT and LC numbers.
        Default is False.

    bmag : str, optional
        Enter a desired bmag to cut the results e.g. '<6'.

    vmag : str, optional
        Enter a desired vmag to cut the results e.g. '>8.5'.

    gaia : str, optional
        If 'DR2'/'DR3', it will create separate columns with Gaia DR2/DR3 data.
        Default is False.

    radius : int,float
        Enter the search radius in arcsec for the Gaia query. Default is 1.

    skip : str, optional
        Enter a coma separated list of targets to exclude in the table.

    Returns
    -------
    Nothing, but saves the generated table in the */table/ folder.
    '''

    if db == 'IACOB':
        lst_sources_all = findstar(spectra=list_name, snr=snrcut)
        # For each source, the best available SNR is picked
        lst_sources_f = spec_snr(lst_sources_all, get_HF=True)
        type_list = 'names'
    elif db == 'Simbad':
        lst_sources_f = findlist(list_name)
        type_list = input('The list contains names or coordinates? [names/coords]: ').strip()
        if type_list not in ('names', 'coords'):
            msg.error('Input answer is not valid. Exiting... \n')
            return None
    else:
        msg.error('Database not recognised. Exiting... \n')
        return None

    # Tuning the format of the input variables
    if limdist is not None:
        RADEC = limdist[0]
        dist = float(limdist[1])

    if spt is not None:
        spt = re.split(' |,', spt)
    if lc is not None:
        lc = re.split(' |,', lc)

    #=============================== Query Gaia ================================
    if gaia in ('DR2', 'DR3'):
        # Vizier - Only used to query in Gaia
        from astroquery.vizier import Vizier
        if gaia == 'DR2':
            gaia_columns = ['BPmag', 'e_BPmag', '+Gmag', 'e_Gmag', 'RPmag', 'e_RPmag',
                            'pmRA', 'e_pmRA', 'pmDE', 'e_pmDE', 'Plx', 'e_Plx',
                            'astrometric_n_good_obs_al', 'astrometric_chi2_al']
            table_u0 = findtable('table_u0_g_col.txt', delimiter=',')
            offset = input('Apply +0.03 mas offset to parallax? [y/n]: ').strip().lower()
            v = Vizier(columns=gaia_columns)
        else:
            v = Vizier()
        v.ROW_LIMIT = 1

    # Progress Bar
    bar = pb.ProgressBar(maxval=len(lst_sources_f),
        widgets=[pb.Bar('=', '[', ']'), ' ', pb.Percentage(), ' ', pb.ETA()],
        redirect_stdout=True, redirect_stderr=True)
    bar.start()

    #================================= Sources =================================
    table = None
    for i, source in enumerate(lst_sources_f):
        row = Table()
        #============= Retrieve the key values from the fits header ============
        OBJRA = OBJDEC = None

        if db == 'IACOB':
            source_path = Path(source)
            filename = source_path.name
            id_star = filename.split('_')[0]

            with fits.open(source_path) as hdu:
                hdu.verify('fix')       # Fix header keywords
                header = hdu[0].header  # Read the values of the headers

            row['ID'] = [id_star]
            # Official name (not implemented yet)
            row['Name'] = [header.get('I-Name', id_star)]

            # Gather the coordinates:
            if '_M_' in filename:
                OBJRA = header.get('OBJ_RA')
                OBJDEC = header.get('OBJ_DEC')
            elif '_N_' in filename:
                OBJRA = header['OBJRA'] * 360 / 24
                OBJDEC = header['OBJDEC']
            elif '_F_' in filename:
                OBJRA = header['RA']
                OBJDEC = header['DEC']

            source = row['ID'][0]

        elif db == 'Simbad':
            row['ID'] = [source]
            source = row['ID'][0]

        #=========================== Skip bad sources ==========================
        if db == 'IACOB' and any(bad in source[:3] for bad in ['DO2']):
            continue

        if skip is not None and any(bad in source for bad in skip.split(',')):
            continue

        #=============== Simbad query by object name/coordinates ===============
        if type_list == 'names':
            simbad = query_Simbad(source, OBJRA, OBJDEC)
            if simbad is None or len(simbad) == 0:
                continue
        else:
            c = SkyCoord(source, unit=(u.hour, u.deg) if ':' in source else u.deg)
            try: # The actual query
                simbad = query_Simbad(c.ra.deg, c.dec.deg, radius=radius * u.arcsec)
            except Exception: # Retry after 2 seconds if the query fails (e.g. timeout)
                time.sleep(2)
                msg.info(f'Trying using a larger radius of {2 * radius} arcsec')
                simbad = query_Simbad(c.ra.deg, c.dec.deg, radius=2 * radius * u.arcsec)

            if simbad is None or len(simbad) == 0:
                msg.warn(f"Source {source} not found in Simbad")
                continue
            if len(simbad) > 1:
                simbad = Table(simbad[0])

            source = simbad['MAIN_ID'][0]

            row['ID'] = row['Name'] = [source]

        #========================= Get the coordinates =========================
        if coords == 'header':
            RADEC_0 = SkyCoord(ra=header['RA'], dec=header['DEC'], unit=u.deg)
        else:
            RADEC_0 = SkyCoord(ra=simbad['ra'], dec=simbad['dec'], unit=(u.hour, u.deg))[0]

        row['RA_J2000'] = RADEC_0.ra.to_string(unit=u.hour, sep=':', pad=True, alwayssign=True)
        row['DEC_J2000'] = RADEC_0.dec.to_string(unit=u.deg, sep=':', pad=True, alwayssign=True)
        row['RAdeg_J2000'] = round(RADEC_0.ra.deg, 6)
        row['RAdeg_J2000'].unit = u.deg
        row['DECdeg_J2000'] = round(RADEC_0.dec.deg, 6)
        row['DECdeg_J2000'].unit = u.deg

        #========================= Limit by distance ===========================
        if limdist is not None:
            RADEC_f = SkyCoord(RADEC, unit=(u.hour, u.deg) if any(j in RADEC for j in (':', 'h')) else u.deg)
            if RADEC_0.separation(RADEC_f).deg > dist:
                continue

        #============================ Limit by l,b =============================
        if lim_lb is not None:
            l_0 = RADEC_0.galactic.l.deg
            b_0 = RADEC_0.galactic.b.deg
            bounds = [float(val) for val in lim_lb]
            if not (bounds[0] < l_0 < bounds[1] and bounds[2] < b_0 < bounds[3]):
                continue

        #===================== Get the photometric magnitudes ==================
        # Photometry
        for mag_band in ('U', 'B', 'V'):
            val = simbad[mag_band][0]
            val = np.nan if str(val) == '--' else float(val)
            row[f'mag_{mag_band}'] = round(val, 2)
            row[f'mag_{mag_band}'].unit = u.mag

        if bmag is not None and not np.isnan(row['mag_B'][0]):
            b_val = row['mag_B'][0]
            if bmag.startswith('<') and float(bmag[1:]) < b_val:
                continue
            if bmag.startswith('>') and float(bmag[1:]) > b_val:
                continue

        if vmag is not None and not np.isnan(row['mag_V'][0]):
            v_val = row['mag_V'][0]
            if vmag.startswith('<') and float(vmag[1:]) < v_val:
                continue
            if vmag.startswith('>') and float(vmag[1:]) > v_val:
                continue

        #======================= Get the spectral class ========================
        if db == 'IACOB':
            SpC_0 = header.get('I-SPC', '-')
            if not isinstance(SpC_0, str) or SpC_0.strip() == '-':
                row['SpC'] = SpC_0 = simbad['sp_type'][0].replace(' ', '')
                row['SpC_ref'] = 'SIMBAD'
            else:
                row['SpC'] = SpC_0.replace(' ', '')
                row['SpC_ref'] = header.get('I-SPCREF', '???')
        else:
            row['SpC'] = SpC_0 = simbad['sp_type'][0].replace(' ', '')
            row['SpC_ref'] = 'SIMBAD'

        #========================= Limit by SpT or LC ==========================
        if spt is not None or lc is not None:
            if 'I' in SpC_0:
                spt_0 = SpC_0[:SpC_0.index('I')]
                lc_0 = SpC_0[SpC_0.index('I'):]
            elif 'V' in SpC_0:
                spt_0 = SpC_0[:SpC_0.index('V')]
                lc_0 = SpC_0[SpC_0.index('V'):]
            elif len(re.split(r':|p|n|f|\(|\)', SpC_0.strip())[0]) <= 4:
                spt_0 = SpC_0.strip()
                lc_0 = ''
            else:
                spt_0 = SpC_0
                lc_0 = ''

            if spt is not None and not any(s in spt_0 for s in spt):
                continue

            if lc is not None and lc_0:
                match_lc = False
                for token in ('IV', 'V', 'III', 'II', 'I'):
                    if token in lc and token in lc_0:
                        match_lc = True
                    lc_0 = lc_0.replace(token, '')
                if not match_lc:
                    continue

        #===================== Get the spectral class code =====================
        if spccode:
            row['SpT'], row['LC'] = spc_code(SpC_0) if SpC_0 else (np.nan, np.nan)

        #============== Add extra columns from header or database ==============
        # Counting FIES / HERMES / FEROS spectra
        if db == 'IACOB':
            FIES = HERMES = FEROS = 0
            for item in lst_sources_all:
                item_name = Path(item).name
                if item_name.startswith(f"{source}_") and item_name.endswith('.fits'):
                    if '_N_' in item_name:
                        FIES += 1
                    elif '_M_' in item_name:
                        HERMES += 1
                    elif '_F_' in item_name:
                        FEROS += 1

            row['FIES'] = FIES
            row['HERMES'] = HERMES
            row['FEROS'] = FEROS
            row['Ref_file'] = filename
            row['SNR_best'] = int(header['I-SNR'])
            row['SB'] = header.get('I-SB', '----')
            row['Comments'] = header.get('I-commen', '----')

        #======================== Get Gaia DR2/3 data ==========================
        if gaia in ('DR2', 'DR3'):
            catalog = "I/345/gaia2" if gaia == 'DR2' else "I/355/gaiadr3"
            try:
                gaiaq = v.query_object(source, catalog=catalog, radius=radius * u.arcsec)[0]
                if gaia == 'DR2':
                    if str(gaiaq['Gmag'][0]) != '--':
                        gaiaq['e_Gmag'] = max(gaiaq['e_Gmag'], 8e-3)
                        if 6 <= gaiaq['Gmag'] <= 16:
                            gaiaq['Gmag'] -= 0.0032 * (gaiaq['Gmag'] - 6)
                        elif gaiaq['Gmag'] < 6:
                            gaiaq['Gmag'] += 0.0271 * (6 - gaiaq['Gmag'])
                        gaiaq['Gmag'].unit = u.mag

                    if offset == 'y' and str(gaiaq['Plx'][0]) != '--':
                        gaiaq['Plx'] = round(gaiaq['Plx'] + 0.03, 6)

                    UWE = np.sqrt(float(gaiaq['chi2AL']) / (float(gaiaq['NgAL']) - 5))
                    diff = (abs(gaiaq['Gmag'][0] - table_u0['g_mag']) +
                            abs(gaiaq['BPmag'][0] - gaiaq['RPmag'][0] - table_u0['bp_rp'])).tolist()
                    gaiaq['RUWE'] = round(UWE / table_u0['u0'][diff.index(min(diff))], 4)
                    gaiaq.remove_columns(['NgAL', 'chi2AL'])
                elif gaia == 'DR3':
                    gaiaq = gaiaq[[c for c in gaiaq.columns if not any(x in c for x in ['dr2', 'scan_', '_corr', '_obs_', 'ipd_', 'transits'])]]

                row = hstack([row, gaiaq])
            except Exception as e:
                msg.warn(f"{source} could not be queried in Gaia due to: {e}.")

        table = row if table is None else vstack([table, row])
        bar.update(i)

    bar.finish()

    if table is None or len(table) == 0:
        msg.warn('Table is empty, no sources were found.')
        return None

    msg.info(f"Table has been successfully created with length {len(table)}")

    table['mag_U'] = table['mag_U'].astype('float32')
    table['mag_B'] = table['mag_B'].astype('float32')
    table['mag_V'] = table['mag_V'].astype('float32')
    if 'SB' in table.colnames:
        table['SB'] = table['SB'].astype('<U8')
    if 'Comments' in table.colnames:
        table['Comments'] = table['Comments'].astype('<U32')

    out_file = maindir / 'tables' / 'tablestars.fits'
    out_file.parent.mkdir(parents=True, exist_ok=True)
    table.write(out_file, format='fits', overwrite=True)
    return None


def spc_code(spc, sub_class_I=False):
    '''
    Function to translate input SpC (Spectral Class) into SpT and LC codes.

    Parameters
    ----------
    spc : str
        Enter the spectral class to turn into SpT and LC codes.

    Returns
    -------
    SpT and LC codes.
    '''

    spt_dic = {'O': 1, 'B': 2, 'A': 3, 'F': 4, 'G': 5, 'K': 6, 'M': 7}
    if not sub_class_I:
        lc_dic = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5}
    else:
        lc_dic = {'Iab': 1.2, 'Ia+': 0.9, 'Ia': 1.1, 'Ib': 1.3, 'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5}

    spc = spc.replace('var', '').replace('Mn', '').replace('Ve', 'V').replace('Va', 'V')
    spc = spc.replace('/ab', '/Iab').replace('/a', '/Ia').replace('/b', '/Ib')
    spc = spc.split('+')[0].replace('/', '-').replace(':', '')

    if spc in ['~'] or spc.startswith(('WC', 'WN', 'WR', 'C')):
        return np.nan, np.nan

    # Spectral type
    spt_c = re.findall('[O,B,A,F,G,K,M]', spc)
    num = len(spt_c)
    spt_c_lst = [spt_dic[i] for i in spt_c]

    spt_digits = re.findall(r'[0-9.]+', spc)
    spt_c_lst += [float(i) / 10 for i in spt_digits]

    try:
        # This is for B1-2 in which there is B(1) and 1-2(2)
        if num < len(spt_c_lst) - num:
            spt_val = (sum(spt_c_lst[:num]) / num + sum(spt_c_lst[num:]) / len(spt_c_lst[num:]))
        else:
            spt_val = sum(spt_c_lst) / num
    except Exception:
        msg.warn(f"SpC separation into SpT and LC code failed for: {spc}, {spt_c_lst}, {num}")
        spt_val = np.nan

    # Luminosity class
    pattern = '[Iab, Ia+, Ia, Ib, I, V]+' if sub_class_I else '[I,V]+'
    lc_matches = re.findall(pattern, spc)
    lc_c_lst = [lc_dic[i] for i in lc_matches if i in lc_dic]

    lc_val = np.asarray(lc_c_lst).mean() if lc_c_lst else np.nan

    return round(spt_val, 2), round(lc_val, 2)


def query_Simbad(name=None, ra=None, dec=None, radius='5s', otypes=False, other_names=False):
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

    radius : str, optional
        Enter a string with the radius for the sky search.

    otypes : boolean, optional
        If True, the object types will be also queried. Default is False.
        Note that sometimes this issues some errors.

    other_names : boolean, optional
        If True, it asks Simbad for all the other names of the source.
        Default is False.

    Returns
    -------
    Queried object in Table format.
    '''

    # if SIMBAD module is not imported, import it
    if 'Simbad' not in globals():
        from astroquery.simbad import Simbad
        Simbad.server = 'simbad.harvard.edu'
        Simbad.add_votable_fields("U", "B", "V", "sp_type")
        Simbad.columns_in_output = [i for i in Simbad.columns_in_output if 'coo' not in i.name]

    if otypes and 'otypes' not in Simbad.columns_in_output:
        Simbad.add_votable_fields('otypes')
        Simbad.columns_in_output = [i for i in Simbad.columns_in_output if i.name not in ['origin', 'otype_txt']]

    if name is not None:
        if other_names:
            simbad = Simbad.query_objectids(name)
            if simbad and len(simbad) > 0:
                for entry in simbad:
                    print(entry['ID'])
            else:
                print(f"No objects found for {name} in Simbad.")
            return None

        # For some reason sometimes adding a whitespace fix some querying issues
        simbad = None
        for name_variant in (name, name + ' ', ' ' + name):
            simbad = Simbad.query_object(name_variant)
            if simbad is not None and len(simbad) > 0:
                break

        while simbad is None or len(simbad) == 0:
            print(f"Provide alternative name for {name} in Simbad.")
            print('In some cases try replacing "HD" by "HD " or vice versa.')
            if ra is not None and dec is not None:
                print(f'Type "sky" to query {radius} around input ra/dec.')
            check = input('Alternative (or hit enter to skip): ').strip()
            if not check:
                print(f"Skipping source: {name}\n")
                break
            if check == 'sky' and ra is not None and dec is not None:
                simbad = Simbad.query_region(SkyCoord(ra, dec, unit='deg'), radius=radius)
            else:
                try:
                    simbad = Simbad.query_object(check)
                except Exception:
                    simbad = None

        if otypes and simbad is not None and len(simbad) > 0:
            # list all rows with different main_id
            main_ids = list(set(simbad['main_id']))
            for main_id in main_ids:
                simbad_i = simbad[simbad['main_id'] == main_id]
                # merge all the otypes in one string
                all_otypes = [simbad_i['otypes.otype'][i] for i in range(len(simbad_i))]
                # remove all rows with the same main_id and add one with all the otypes
                simbad = simbad[simbad['main_id'] != main_id]
                new_row = simbad_i[0]
                new_row['otypes.otype'] = '|'.join(all_otypes)
                simbad.add_row(new_row)

        if simbad is not None and len(simbad) > 1:
            msg.warn('More than one Simbad result, choosing the brightest source...')
            simbad.sort('V')
            simbad = Table(simbad[0])

        if simbad is not None and 'SCRIPT_NUMBER_ID' in simbad.colnames:
            simbad.remove_column('SCRIPT_NUMBER_ID')

        return simbad

    if ra is not None and dec is not None:
        simbad = Simbad.query_region(SkyCoord(ra, dec, unit='deg'), radius=radius)
        if simbad is None or len(simbad) == 0:
            msg.warn('No objects found.')
            return None
        if len(simbad) > 1:
            msg.warn('More than one Simbad result, choosing the brightest source...')
            simbad.sort('V')
            simbad = Table(simbad[0])
        return simbad

    msg.error('Name or RA+DEC is required for the query...')
    return None


def query_Gaia(name=None, gaia='dr3', ra=None, dec=None, radec=None, radius=2,
               limit_one_source=True, get_zp=False):
    '''
    Function to query an object in Gaia database.

    Parameters
    ----------
    name : str, optional
        Enter the name (ID) of the source to query.

    gaia : str, optional
        Enter the Gaia data you want to query, either 'dr2' or 'dr3' (default)

    ra : float, optional
        Enter the right ascension of the source, in degrees.

    dec : float, optional
        Enter the declination of the source, in degrees.

    radec : str, optional
        Enter the string with the coordinates as hh(:)mm(:)ss.sss(:) +-dd(:)mm(:)ss(:)

    radius : int, optional
        Enter an integer with the radius for the sky search in arcseconds. Default is 2.

    limit_one_source : boolean, optional
        If True, only one source will be returned. Default is True.

    get_zp : boolean, optional
        If True, the zero point will be also calculated. Default is False.

    Returns
    -------
    Queried object in Table format.
    '''

    # Gaia Query
    from astroquery.gaia import Gaia

    # Gaia Zero point offset
    from edr3_zp import zpt

    Gaia.ROW_LIMIT = -1
    zpt.load_tables()

    gaia_lower = gaia.lower()
    if gaia_lower not in ('dr2', 'dr3'):
        msg.error('Input Gaia catalog is not dr2 / dr3. Exiting...')
        return None

    Gaia.MAIN_GAIA_TABLE = f"gaia{gaia_lower}.gaia_source"

    if name is not None:
        RADEC = SkyCoord.from_name(name)
        while RADEC is None:
            print(f"Provide alternative name for {name} in Simbad.")
            print('In some cases try replacing "HD" by "HD " or vice versa.')
            check = input('Alternative (or hit enter to skip): ').strip()
            if not check:
                print(f"Skipping source: {name}\n")
                break
            try:
                RADEC = SkyCoord.from_name(check)
            except Exception:
                RADEC = None
    elif ra is not None and dec is not None:
        RADEC = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
    elif radec is not None:
        RADEC = SkyCoord(radec, unit=(u.hour, u.deg), frame='icrs')
    else:
        msg.error('Name or RA+DEC is required for the query...')
        return None

    query = Gaia.cone_search_async(coordinate=RADEC, radius=u.Quantity(radius, u.arcsec))
    query = query.results

    if len(query) == 0:
        msg.error(f"Gaia query failed for target {name or (RADEC.ra, RADEC.dec)}")
        return None

    if limit_one_source and len(query) > 1:
        msg.warn('More than one Gaia result')
        if name in query['designation']:
            msg.info('Choosing the same Gaia source ID as the input...')
            query = query[query['designation'] == name]

        else:
            msg.warn('Input name not in Gaia source ID, choosing the brigtest source...')
            query.sort('phot_g_mean_mag')
            query = Table(query[0])

    if get_zp and gaia_lower == 'dr3':
        if query['astrometric_params_solved'] <= 3:
            msg.info(f'Astrometric parameters solved is <3. Skipping source_id = {query["source_id"]} for zero point calculation.')
            query['zp_offset'] = -9999
        else:
            query['zp_offset'] = zpt.get_zpt(
                query['phot_g_mean_mag'], query['nu_eff_used_in_astrometry'],
                query['pseudocolour'], query['ecl_lat'], query['astrometric_params_solved']
            )

    query['ID'] = name if name is not None else '--'
    for extra_col in ('designation', 'phot_variable_flag', 'libname_gspphot'):
        if extra_col in query.colnames:
            query.remove_column(extra_col)

    return query


def show_header(fitsfile):
    '''
    Function to show the header of a fits file separated in rows.

    Parameters
    ----------
    fitsfile : str
        Enter the name of the fits file.

    Returns
    -------
    Nothing, but the header is printed.
    '''

    matches = findstar(fitsfile)
    if not matches:
        msg.error('File not found.')
        return None
    if len(matches) > 1:
        msg.error('Only one file can be shown at a time.')
        return None

    with fits.open(Path(matches[0])) as hdul:
        hdul.info()
        print(hdul[0].header.items)
        return hdul[0].header


def check_fits(list_name, max_dist=90):
    '''
    Function to detect errors in the filenames and headers of the fits files.

    Parameters
    ----------

    list_name : str
        Enter the input list, either name(s)/FITS of the source(s) separated by coma,
        or a .txt/.lst file containing the source names or files.

    max_dist : int/float
        Enter the threshold distance in arcseconds to consider the spectra to point
        to a different source, and not caused by telescope pointing offset.
        Default is 90 arcseconds.

    Returns
    -------
    Nothing, but a file with the errors found is generated.
    '''

    dir_spectra = findstar(list_name)
    if not dir_spectra:
        return None

    bar = pb.ProgressBar(maxval=len(dir_spectra),
                         widgets=[pb.Bar('=', '[', ']'), ' ', pb.Percentage(), ' ', pb.ETA()])
    bar.start()

    type_errors = {'Simbad': [], 'filename': [], 'duplicate': [], 'baddate': [], 'radec0': [], 'distance': []}
    IDs = {'F': [], 'M': [], 'N': []}
    date_0_long = []

    for i, spectrum in enumerate(dir_spectra):
        bar.update(i)
        p = Path(spectrum)
        filename = p.name
        id_star = filename.split('_')[0]

        # Problems quering in Simbad
        simbad = query_Simbad(id_star) if id_star not in type_errors['Simbad'] else []
        if not simbad and id_star not in type_errors['Simbad']:
            type_errors['Simbad'].append(id_star)

        with fits.open(p) as hdu:
            hdu.verify('fix')
            header = hdu[0].header

        if '_F_' in filename:
            IDs['F'].append(header['ARCFILE'])
        elif '_M_' in filename:
            IDs['M'].append(str(header['UNSEQ']))
        elif '_N_' in filename:
            IDs['N'].append(header['FILENAME'])

        # Catch files with wrong object name comparing filename and header
        if id_star != header['OBJECT'].replace(' ', '').upper():
            type_errors['filename'].append(filename)

        # Get the date and time from the header
        if '_N_' in filename or '_M_' in filename:
            date_0 = header['DATE-AVG']
            date_0_short = str(date_0[0:4]) + str(date_0[5:7]) + str(date_0[8:10])
        elif '_F_' in filename:
            date_0 = header['ARCFILE']
            date_0_short = str(date_0[6:10]) + str(date_0[11:13]) + str(date_0[14:16])

        # Get the date from the filename
        date_filename = filename.split('_')[1]
        # Catch two+ files with same date in header but different filenames
        # e.g. HD111111 and HDE111111
        if date_0 not in date_0_long:
            date_0_long.append(date_0)
        else:
            type_errors['duplicate'].append(filename)

        # Catch wrong dates in the filenames when comparing with the date in the header
        if date_filename != date_0_short:
            type_errors['baddate'].append(filename)

        # Get the coordinates from the header
        RA_0, DEC_0 = header['RA'], header['DEC']
        # Catch missing RA/DEC:
        if '0.0000' in str(RA_0) or '0.0000' in str(DEC_0):
            type_errors['radec0'].append(filename)

        # Catch spectra of different object than the object queried in Simbad
        if simbad and ('0.0000' not in str(RA_0) and '0.0000' not in str(DEC_0)):
            c1 = SkyCoord(RA_0, DEC_0, unit=u.deg)
            c2 = SkyCoord(simbad['ra'][0], simbad['dec'][0], unit=u.deg)
            difcoord = round(c1.separation(c2).arcsec, 3)
            if difcoord > max_dist:
                type_errors['distance'].append(f"{filename} -> {difcoord}")

    bar.finish()

    # Wrtiting the errors to the file:
    error_descriptions = [
        '# Problem querying these objects in Simbad. Correct the name and re-run to check object coordinates.\n',
        '# File names which has a wrong name for the object when compared to the header.\n',
        '# Duplicated files with same full date but different name in the file name.\n',
        '# Files with different date in the header and the file name.\n',
        '# Spectra with no RA or DEC values in the header. Correct the header and re-run to check object coordinates.\n',
        '# Distance from Simbad query is above the threshold for these spectra.\n'
    ]

    errors_path = maindir / 'lists' / 'Errors_DB.txt'
    errors_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Possible errors (if any) are written in {errors_path}")
    with errors_path.open('w') as errorsDB:
        for description, error in zip(error_descriptions, type_errors):
            if not type_errors[error]:
                continue
            errorsDB.write(description)
            for item in type_errors[error]:
                errorsDB.write(f"{item}\n")

    ids_path = maindir / 'lists' / 'IDs_spectra.txt'
    ids_path.parent.mkdir(parents=True, exist_ok=True)
    # Writing the file unique identificators to a file
    print(f"Unique IDs of the spectra are written in {ids_path}")
    with ids_path.open('w') as IDs_spectra:
        for description, tel in zip(['# FEROS files\n', '# HERMES files\n', '# FIES files\n'], IDs):
            if not IDs[tel]:
                continue
            IDs_spectra.write(description)
            for item in IDs[tel]:
                IDs_spectra.write(f"{item}\n")


def fix_fits(list_name, coords_Simbad=True, radius=60, savepath=''):
    '''
    Function to detect and correct errors in the filenames/headers allowing
    the user to save a corrected FITS file.
    Run this function for the stars with errors from the checknames function.

    Parameters
    ----------
    list_name : str
        Enter the input list, either name(s)/FITS of the source(s) separated by coma,
        or a .txt/.lst file containing the source names or files.

    coords_Simbad : boolean, optional
        If True, the coordinates of the header are queried in Simbad to find sources
        within the radius parameter. Default is True.

    radius : int
        Search radius where to find the source by the coordinates, in arcsec.
        Default is 60 arcseg.

    savepath : str, optional
        If set, updated fits files will be saved in that location (full path must be set).

    Returns
    -------
    Nothing.
    '''

    dir_spectra = findstar(list_name)
    if not dir_spectra:
        return None

    for spectrum in dir_spectra:
        spectrum_path = Path(spectrum)
        filename = spectrum_path.name
        id_star = filename.split('_')[0]

        with fits.open(spectrum_path) as hdu:
            hdu.verify('fix')
            header = hdu[0].header

            print(f"\nName (filename) / Name (object)\n{id_star} / {header['OBJECT']}\n")

            RADEC = SkyCoord(str(header['RA']), str(header['DEC']), unit=u.deg)
            if 'OBJ_RA' in header:
                RADEC_OBJ = SkyCoord(str(header['OBJ_RA']), str(header['OBJ_DEC']), unit=u.deg)
            elif 'OBJRA' in header:
                RADEC_OBJ = SkyCoord(str(header['OBJRA']), str(header['OBJDEC']), unit=(u.hour, u.deg))
            else:
                RADEC_OBJ = RADEC

            print('From header:\n--------------------------')
            print(f"RADEC / RADEC_OBJ\n{RADEC.ra} / {RADEC_OBJ.ra}")
            print(f"{RADEC.dec} / {RADEC_OBJ.dec}\n")

            simbad = query_Simbad(id_star)
            if simbad is not None and len(simbad) > 0:
                print('From Simbad (querying ID):\n--------------------------')
                RADEC_SB = SkyCoord(simbad['ra'], simbad['dec'], unit=(u.hour, u.deg))[0]
                try:
                    RADEC_SC = SkyCoord.from_name(id_star)
                except Exception:
                    RADEC_SC = SkyCoord('0 0', unit=u.deg)

                print(f"RADEC (Simbad) / RADEC (SkyCoord)\n{RADEC_SB.ra} / {RADEC_SC.ra}")
                print(f"{RADEC_SB.dec} / {RADEC_SC.dec}\n")
                print(f"SpC (header) / SpC (Simbad)\n{header.get('I-SPC')} ({header.get('I-SPCREF')}) / {simbad['sp_type'][0]}\n")

            if coords_Simbad:
                region = query_Simbad(ra=RADEC.ra, dec=RADEC.dec, radius=radius * u.arcsec)
                if region is not None:
                    print('From Simbad (querying header RA,DEC):\n-------------------------------------')
                    print('Sources found:')
                    for src in region:
                        print(f"{src['MAIN_ID']} ({src['sp_type']})")
                    print('\n')

            user = '-'
            change = False
            while user:
                print("To change something from the header type: <keyword>:<value>")
                print("Type +<name> to obtain a list of alternative names for <name> in Simbad.")
                user = input("Or hit return to continue. Input: ").strip()

                if user.startswith('+') and ':' not in user:
                    query_Simbad(name=user.replace('+', ''), other_names=True)

                if ':' in user:
                    k, v = user.split(':', 1)
                    header[k.strip()] = v.strip()
                    change = True

            if change:
                save = input('FITS file has been modified, do you want to save it? y/n/<filename.fits>: ').strip()
                out_path = spectrum_path.with_name(save) if save.endswith('.fits') else spectrum_path

                if save.lower() != 'n':
                    if savepath:
                        out_path = Path(savepath) / out_path.name

                    out_path.parent.mkdir(parents=True, exist_ok=True)
                    hdu.writeto(out_path, output_verify='ignore', overwrite=True)
