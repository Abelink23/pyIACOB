from db import *

import shutil
import tarfile

from astroquery.eso import Eso
# https://astroquery.readthedocs.io/en/latest/eso/eso.html

eso = Eso()
eso_user = input('Provide the username of your ESO account: ')
eso.login(username=eso_user, store_password=True)
eso.ROW_LIMIT = -1

def get_feros(fname, radius=2, program=None, limit=9999):
    '''
    Parameters
    ----------
    fname : str
        Enter the filename containing the source names [.lst/.txt/.fits(table)].

    radius : int/float, optional
        Max distance in arcmin to the queried result sources. Default is 2 arcmin.

    program : str, optional
        Program ID to search within the ESO archive. Default is None.

    limit : int, optional
        Maximum number of datasets to download per source. Default is 9999.

    Returns
    -------
    None (but it downloads the fits files).
    '''

    if fname.endswith('.lst') or fname.endswith('.txt'):
        sources = findlist(fname)
    elif fname.endswith('.fits'):
        colname = input('Enter the column name with target names: ')
        if colname not in findtable(fname).colnames:
            return print('Column name not found in the table. Check the column names and try again.')
        sources = findtable(fname)[colname]

    # create the FEROS/raw directory if it does not exist
    if not Path(datadir).joinpath('FEROS/raw').exists():
        Path(datadir).joinpath('FEROS/raw').mkdir(parents=True, exist_ok=True)

    db = open(datadir.joinpath('FEROS/raw/downloaded_sources.txt'), 'a')

    msg.info(f'Files will be downloaded to {datadir}/FEROS/raw/')

    bar = pb.ProgressBar(maxval=len(sources), term_width=80, redirect_stdout=True,
            widgets=[pb.Bar('=','[',']'),' ',pb.Percentage()])
    bar.start()

    i = 0
    for source in sources:

        bar.update(i)
        i += 1

        if type(source) is not str:
            name = str(source)

        name = source.strip()

        print('\nSearching data for %s...\n' % name)

        query = eso.query_surveys(surveys='FEROS', cache=False, target=name)
        if query is None:
            db.write(name+', 0, 0, Query returned no results.\n')
            continue

        if program is not None and 'Run/Program ID' in query.colnames:
            query = query[query['Run/Program ID'] == program]

        simbad = query_Simbad(name)
        if simbad is None:
            db.write(name+', 0, 0, Simbad return no results.\n')
            continue

        coords_S = SkyCoord(simbad['ra'], simbad['dec'], unit=(u.hourangle,u.deg), frame='icrs')

        datasets = []
        for entry in query:
            try:
                if re.sub('[ /_/]','',str(entry['Object'])).replace('HD-','HD') == name:
                    datasets.append(entry['ARCFILE'])
                else:
                    coords_F = SkyCoord(entry['RA'], entry['DEC'], unit=(u.deg,u.deg), frame='icrs')
                    if coords_S.separation(coords_F).arcmin[0] < radius:
                        datasets.append(entry['ARCFILE'])
            except:
                db.write('Source %s had problems during the query.\n' % name); continue

        if len(datasets) == 0:
            db.write(name+', 0, 0, No matches in FEROS database.\n'); continue

        if len(datasets) > limit:
            db.write(name+', 0, 0, >%i datasets for the source.\n'); continue

        spectra = 0
        for dataset in datasets:
            miliseconds = round(float(dataset[-5:])+0.001, 3)

            if miliseconds == 1.0:
                db.write('Please check: data from %s is probably missing.\n' % name)

            miliseconds = '%.3f' % miliseconds
            dataset_f = dataset[:22] + miliseconds

            print(dataset_f, dataset)
            # check if the dataset has already been downloaded
            if datadir.joinpath('FEROS/raw', dataset_f+'.fits').exists():
                db.write(f'Dataset {dataset_f} already exists. Skipping download.\n')
                spectra += 1
                continue

            try:
                data_files = eso.retrieve_data(dataset_f, destination=datadir.joinpath('FEROS/raw/'))
            except:
                db.write(f'ERROR: Dataset {dataset_f} could not be retrieved.\n')
                continue

            try:
                datadir.joinpath('FEROS/raw', dataset_f+'.xml').unlink()
            except:
                None

            folder_name = datadir.joinpath('FEROS/raw', dataset_f)

            # Extract the tar file and move the fits files to the FEROS/raw directory
            tar = tarfile.open(str(folder_name)+'.tgz')
            tar.extractall(folder_name)
            tar.close()
            Path(str(folder_name)+'.tgz').unlink()

            fits_files = [file.name for file in Path(folder_name).iterdir() if file.is_file() and file.suffix == '.fits']
            for fit_file in fits_files:
                if not fit_file.endswith('1081.fits'):
                    Path(folder_name).joinpath(fit_file).unlink()
                else:
                    spectra = spectra + 1
                    print(f'Found a valid spectrum: {fit_file}')
                    shutil.move(folder_name.joinpath(fit_file), datadir.joinpath('FEROS/raw/'+fit_file))

            shutil.rmtree(folder_name)

        db.write(name+', '+str(len(datasets))+', '+str(spectra)+',--\n')

    bar.finish()

    db.close()

    return msg.g('\n Data retrieval completed. Check the folder.')
