from spec import *
from scipy.io.idl import readsav


def ib_input(table='IACOB_O9BAs_SNR20.fits', output_name='input_IB.txt'):

    '''
    Function to generate the input table for IACOB broad given an input table with the
    target stars.

    Parameters
    ----------
    table : str, optional
        Enter the input table contaning a column 'ID' with the identifier of the stars.

    output_name : str, optional
        Enter the name for the output table. Default is 'input_IB'.

    Returns
    -------
    Nothing but the IACOB broad input file is generated.
    '''

    table = findtable(table)
    input_IB = findtable(output_name)
    input_IB['i'] = input_IB['i'].astype(str)

    i = 0
    with open(maindir+'tables/input_IB_new.txt', 'w') as file1,\
         open(maindir+'tables/input_IB_updated.txt', 'w') as file2:

        file1.write('#i ID path spectrum resol line QIB\n')
        file2.write('#i ID path spectrum resol line QIB\n')

        bar = pb.ProgressBar(maxval=len(table),
                            widgets=[pb.Bar('=','[',']'),' ',pb.Percentage()])
        bar.start()

        for row,j in zip(table,range(len(table))):

            # Use next line to skip stars not processed in RVEWFW()
            #if row['QSiIII'].tolist() ==  None: continue

            id = row['ID'].strip()

            match_IB = input_IB[input_IB['ID'] == id]

            star = spec(id,SNR='best')

            line = 'SiIII4567' # SiIII4567 / SiII6347

            # Exists in the old table?
            if len(match_IB) == 1:
                    # Exist with the same filename (i.e. is the best one)?
                    if star.filename == match_IB['spectrum']:
                        # If it has been processed, add "#"
                        if not search(id+'_'+match_IB['line'][0]+'.ps', ibdir) == None:
                            idx = '#'+str(j)
                        else:
                            idx = str(j)
                            print('\nWARNING: %s not processed or is missing the .ps\n' %id)

                        file2.write('%s %s %s %s %i %s %i\n' % \
                        (idx,match_IB['ID'][0],match_IB['path'][0],match_IB['spectrum'][0],
                             match_IB['resol'][0],match_IB['line'][0],match_IB['QIB'][0]))

                    # Pick the line used in the last best spectrum
                    else: line = match_IB['line'][0]

            # Then there is a new one or a better one:
            if len(match_IB) == 0 or (len(match_IB) == 1 and not star.filename == match_IB['spectrum']):
                file1.write('%i %s %s %s %i %s %i\n' % \
                (i,star.id_star,star.spectrum.split(star.filename)[0],\
                star.filename,star.resolution,line,-1))
                i = i + 1

                # ...and add it to the new main list.
                file2.write('%s %s %s %s %i %s %i\n' % \
                (str(j),star.id_star,star.spectrum.split(star.filename)[0],\
                star.filename,star.resolution,line,-1))

            bar.update(j)

        bar.finish()

    return 'DONE'


def ib_results(input_table='input_IB.txt', check_best=True, format='fits'):

    '''
    Function to generate a table with the results from IACOB-broad given an input
    table containing the name of the stars.

    Parameters
    ----------
    input_table : str, optional
        Name of the input table contaning the list of stars to search.

    check_best : boolean, optional
        True if each spectra from the .xdr file is checked against the best
        spectrum in the database. Default is True.

    format : str, optional
        Enter the output format for the table: 'fits' (default), 'ascii' or 'csv'.

    Returns
    -------
    Nothing but the output table with the IACOB broad results is generated.
    '''

    table = findtable(input_table)

    bar = pb.ProgressBar(maxval=len(table),
                         widgets=[pb.Bar('=','[',']'),' ',pb.Percentage()])
    bar.start()

    data_rows = Table()
    for row,i in zip(table,range(len(table))):
        id_star = row['ID']
        filename = row['spectrum']
        line = row['line']
        QIB = row['QIB']

        match = []
        for file in os.listdir(ibdir):
            if file.endswith('_resul.xdr') and file.split('_')[0] == id_star:
                match.append(ibdir+file)

        if len(match) == 0:
            print('\nWARNING: No .xdr file found for %s. Continuing...' % id_star)
            continue
        elif len(match) > 1:
            print('\nWARNING: More than one .xdr (line) file found for %s:' % id_star)
            print([xdr.split('%s_' %id_star)[1].split('_resul.xdr')[0] for xdr in match])
            print('Using the one provided in the input table...')

            match = [xdr for xdr in match if line in xdr]

        idldata = readsav(match[0])

        if check_best == True and filename != spec(id_star,SNR='best').filename:
            print('\nWARNING: %s does not match with best spectrum available.'
            % filename)

        data_row = Table(idldata.d)
        data_row['filename'] = filename
        data_row['QIB'] = QIB

        data_rows = vstack([data_rows,data_row])

        parameters = [j for j in idldata.d.dtype.names]

        bar.update(i)

    bar.finish()

    # Renaming the columns and splitting some of them
    data_rows.rename_column('STAR','ID')
    data_rows.rename_column('LINE','line_IB')
    data_rows.rename_column('EW','EW_IB')
    data_rows.rename_column('SNR','SNR_IB')
    data_rows.rename_column('VFT','vsini_FT')
    data_rows.rename_column('VMFT','vmac_FT')
    data_rows['vmac_FT_eDW'] = [i[0] for i in data_rows['EVMFT']]
    data_rows['vmac_FT_eUP'] = [i[1] for i in data_rows['EVMFT']]
    data_rows.rename_column('VSGOF','vsini_GF')
    data_rows['vsini_GF_eDW'] = [i[0] for i in data_rows['EVSGOF']]
    data_rows['vsini_GF_eUP'] = [i[1] for i in data_rows['EVSGOF']]
    data_rows.rename_column('VMGOF','vmac_GF')
    data_rows['vmac_GF_eDW'] = [i[0] for i in data_rows['EVMGOF']]
    data_rows['vmac_GF_eUP'] = [i[1] for i in data_rows['EVMGOF']]

    # Post calculation modifications:
    data_rows['vmac_GF_eUP'][data_rows['vmac_GF'] > 120] = np.nan
    data_rows['vmac_GF_eDW'][data_rows['vmac_GF'] > 120] = np.nan
    data_rows['vmac_GF'][data_rows['vmac_GF'] > 120] = 120
    data_rows['SNR_IB'] = [int(round(row)) for row in data_rows['SNR_IB']]

    # Saving the results:
    names = ['ID','filename','vsini_FT','vmac_FT','vmac_FT_eDW','vmac_FT_eUP','vsini_GF',
             'vsini_GF_eDW','vsini_GF_eUP','vmac_GF','vmac_GF_eDW','vmac_GF_eUP',
             'line_IB','EW_IB','SNR_IB','QIB']

    output = Table(rows=data_rows[names], names=(names))

    full_path = maindir + 'tables/IB_results_new.' + format
    if format == 'ascii':
        format += '.fixed_width_two_line'

    output.write(full_path, format=format, overwrite=True)

    return 'DONE'
