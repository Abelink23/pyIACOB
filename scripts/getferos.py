import sys; sys.path.append('../')

from db import *

import shutil
import tarfile

from astroquery.eso import Eso
# https://astroquery.readthedocs.io/en/latest/eso/eso.html
eso = Eso()
eso.login("Abelink", store_password = True)
eso.ROW_LIMIT = -1


def getferos(fname,colname=None,radius=2,limit=9999):
    '''
    Parameters
    ----------
    fname : str
        Enter the filename containing the source names [.lst/.txt/.fits].

    colname : str, optional
        Column name for the sources if .fits file used as input for table.

    radius : int/float, optional
        Max distance in arcmin to the queried result sources. Default is 2 arcmin.

    limit : int, optional
        Maximum number of datasets to download per source. Default is 9999.
    Returns: None (but it downloads the fits files).
    '''

    path = maindir

    file = search(fname,path)

    if fname.endswith('.lst') or fname.endswith('.txt'):
        sources = findlist(fname)
    elif fname.endswith('.fits'):
        if colname == None:
            while colname == '':
                colname = input('Enter the column name with target names: ')
        sources = findtable(fname)[colname]

    bar = pb.ProgressBar(maxval=len(sources),
                         widgets=[pb.Bar('=','[',']'),' ',pb.Percentage()])
    bar.start()

    db = open(datadir+'FEROS/raw/downloaded_sources.txt','a')

    i = 0
    for source in sources:

        bar.update(i); i = i + 1

        if type(source) != str: name = str(source)
        name = source.strip(); print('\nSearching data for %s...\n' % name)

        query = eso.query_surveys('FEROS',cache=False,target=name)
        if query is None:
            db.write(name+', 0, 0, Query returned no results.\n'); continue

        simbad = SB(name)
        if simbad is None:
            db.write(name+', 0, 0, Simbad return no results.\n'); continue

        coords_S = SkyCoord(simbad['RA'],simbad['DEC'],unit=(u.hourangle,u.deg),frame='icrs')

        datasets = []
        for entry in query:

            try:
                if re.sub('[ /_/]','',str(entry['Object'])).replace('HD-','HD') == name:
                    datasets.append(entry['ARCFILE'])

                else:
                    coords_F = SkyCoord(entry['RA'],entry['DEC'],unit=(u.deg,u.deg),frame='icrs')

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
            miliseconds = round(float(dataset[-5:])+0.001,3)

            if miliseconds == 1.0:
                db.write('Please check: data from %s is probably missing.\n' % name)

            miliseconds = '%.3f' % miliseconds
            dataset_f = dataset[:22] + miliseconds

            #if os.path.exists(datadir+dataset_f) == True: continue

            try:
                data_files = eso.retrieve_data(dataset_f,destination=datadir+'FEROS/raw')
            except:
                db.write('ERROR: Dataset %s could not be retrieved.\n' %dataset_f); continue

            try:
                os.remove(datadir+'FEROS/raw/'+dataset_f+'.xml')
            except:
                None#;db.write('Could not remove dataset %s.xml\n' % dataset_f)

            folder_name = datadir+'FEROS/raw/'+dataset_f

            try:
                tar = tarfile.open(folder_name+'.tgz')
                tar.extractall(folder_name)
                tar.close()
                os.remove(folder_name+'.tgz')
            except:
                db.write('Tar file for %s had problems during download. Denied access?\n' % name)
                continue


            fits_files = [file for file in os.listdir(folder_name) if file.endswith('.fits')]
            for fit_file in fits_files:
                if not fit_file.endswith('1081.fits'):
                    os.remove(folder_name+'/'+fit_file)
                else:
                    spectra = spectra + 1
                    shutil.move(folder_name+'/'+fit_file,datadir+'FEROS/raw/'+fit_file)

            shutil.rmtree(folder_name)

        db.write(name+', '+str(len(datasets))+', '+str(spectra)+',--\n')

    bar.finish()

    db.close()

    return print('\n Data retrieval completed. Check the folder.')
