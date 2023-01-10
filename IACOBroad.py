from spec import *
from scipy.io.idl import readsav


def ib_input(table, folder, old_input=None, new_master=False, txt=False):

    '''
    Function to generate the input table for IACOB broad given an input table with the
    target stars or specific spectrum.

    Parameters
    ----------
    table : str
        Enter the input table contaning a column 'ID' or 'Ref_file' with the identifier
        of the stars. If the ID is provided, then the best SNR spectrum is used for each
        star (prioritizing HERMES+FEROS).

    folder : str
        Name of the sub-folder under the IACOB-Broad directory (ibdir) where to
        search for the solution files.

    old_input : str, optional
        Name of a previous IACOB-Broad input table. If used, it will be used to create
        an updated verstion of it. Default is None.

    new_master : boolean, optional
        If True, it will create a new master table, adding the new rows to the old input
        table, or updating the old one if that one has to be redone. Default is False.

    txt : boolean, optional
        If True, it will parse it to spec() in case that the input spectrum is in ascii.

    Returns
    -------
    Nothing but the IACOB broad input file called 'input_IB_todo.txt' is generated.
    If old_input is provided, it will also generate a table called 'input_IB_todo-updated.txt'
    with updated rows in case that the spectrum has been already processed.
    '''

    table = findtable(table)
    if 'Ref_file' in table.colnames:
        table['ID'] = table['Ref_file']

    file1 = open(maindir+'tables/input_IB_todo.txt', 'w')
    file1.write('#i ID path Ref_file resol line QIB\n')

    if old_input is not None:
        old_input = findtable(old_input)
        old_input['i'] = old_input['i'].astype(str)
        file2 = open(maindir+'tables/input_IB_todo-updated.txt', 'w')
        file2.write('#i ID path Ref_file resol line QIB\n')

    i = 0

    bar = pb.ProgressBar(maxval=len(table),
                        widgets=[pb.Bar('=','[',']'),' ',pb.Percentage()])
    bar.start()

    for row,j in zip(table,range(len(table))):
        
        id = row['ID'].strip()
        
        star = spec(id, SNR='bestHF', txt=txt)
        
        star.get_spc()
        
        if star.SpC != '':
            spt = spc_code(star.SpC)[0]
            if not np.isnan(spt):
                if spt < 2.4:
                    line = 'SiIII4567'
                elif spt >= 2.4:
                    line = 'SiII6371' # 6347 Has MgII contamination
            else:
                line = 'SiIII4567'
        else:
            line = 'SiIII4567'
        
        # Exists in the old table?
        if old_input is not None:
            
            if 'Ref_file' in row.colnames:
                match_old = old_input[old_input['Ref_file'] == star.filename]
            else:
                match_old = old_input[old_input['ID'] == star.id_star]
            
            if len(match_old) >= 1:
            
                # If Ref_file is used as input, but there are several matches in old table:
                if 'Ref_file' in row.colnames:
                    match_old = match_old[match_old['Ref_file'] == star.filename]
                # If ID is used as input, and the best SNR spectrum is in the old table:
                elif star.filename in match_old['Ref_file']:
                    match_old = match_old[match_old['Ref_file'] == star.filename]
                # If ID is used as input, and the best SNR spectrum is not in the old table:
                # Despite the Ref_file is different, picks the old line used there
                else:
                    line = match_old[0]['line'][0]
            
            # For the same Ref_file as in the old table...
            if len(match_old) == 1:
                # If it has been already processed, add "#"
                if not search(star.id_star+'_'+match_old['line'][0]+'.ps', ibdir+folder) == None:
                    idx = '#'+str(j)
                else:
                    idx = str(j)
                    print('\nWARNING: %s not processed or is missing the .ps\n' % id)
                file2.write('%s %s %s %s %i %s %i\n' % \
                (idx,match_old['ID'][0],match_old['path'][0],match_old['Ref_file'][0],
                match_old['resol'][0],match_old['line'][0],match_old['QIB'][0]))
            
            # Then there is a new one or a better one:
            elif len(match_old) == 0:
                file2.write('%s %s %s %s %i %s %i\n' % \
                (str(j),star.id_star,star.fullpath.split(star.filename)[0],\
                star.filename,star.resolution,line,-1))
        
        file1.write('%i %s %s %s %i %s %i\n' % \
        (i,star.id_star,star.fullpath.split(star.filename)[0],\
        star.filename,star.resolution,line,-1))
        i = i + 1
        
        bar.update(j)
        
        bar.finish()
    
    if old_input is not None:
        file2.close()

    file1.close()
    
    # If new_master is True, it will create a new master table with the new input
    if old_input is None and new_master == True:
        print('Cannot create a new master table without an old one.\n')
    elif old_input is not None and new_master == True:
        table = findtable('input_IB_updated.txt')
        table['i'] = table['i'].astype(str)
        for row in table:
            if row['Ref_file'] in old_input['Ref_file']:
                old_input_i = old_input[old_input['Ref_file'] == row['Ref_file']]
                # To catch potential issues (not really necessary)
                if row['QIB'] != old_input_i['QIB'] or row['line'] != old_input_i['line']:
                    print('Something went wrong for:', row['Ref_file'], row['QIB'], row['line'], old_input_i['QIB'], old_input_i['line'])
                    continue
            else:
                old_input = vstack([old_input, row])
        
        old_input['i'] = np.arange(len(old_input)).astype(str)
        
        file3 = open(maindir+'tables/input_IB_new.txt', 'w')
        file3.write('#i ID path Ref_file resol line QIB\n')
        for row in old_input:
            file3.write('%s %s %s %s %s %s %i\n' % \
            (row['i'],row['ID'],row['path'],row['Ref_file'],row['resol'],row['line'],row['QIB']))
        file3.close()
    
    return 'DONE'


def ib_results(table, folder, check_best=False, format='fits'):

    '''
    Function to generate a table with the results from IACOB-broad given an input
    table containing the name of the stars.

    Parameters
    ----------
    table : str
        Name of the input table contaning the list of stars to search.

    folder : str
        Name of the sub-folder under the IACOB-Broad directory (ibdir) where to
        search for the solution files.

    check_best : boolean, optional
        True if each spectra from the .xdr file is checked against the best
        spectrum in the database. Default is False.

    format : str, optional
        Enter the output format for the table: 'fits' (default), 'ascii' or 'csv'.

    Note: the option for using this over ascii input spectra is not yet implemented.

    Returns
    -------
    Nothing but the output table with the IACOB broad results is generated.
    '''

    table = findtable(table)

    bar = pb.ProgressBar(maxval=len(table),
                         widgets=[pb.Bar('=','[',']'),' ',pb.Percentage()])
    bar.start()

    data_rows = Table()
    for row,i in zip(table,range(len(table))):
        id_star = row['ID']
        filename = row['Ref_file']
        line = row['line']
        QIB = row['QIB']

        match = []
        if not folder.endswith('/'):
            folder += '/'
        for file in os.listdir(ibdir+folder):
            if file.endswith('_resul.xdr') and file.split('_')[0] == id_star:
                match.append(ibdir+folder+file)

        if len(match) == 0:
            print('\nWARNING: No .xdr file found for %s. Continuing...' % id_star)
            continue
        elif len(match) > 1:
            print('\nWARNING: More than one .xdr (line) file found for %s:' % id_star)
            print([xdr.split('%s_' %id_star)[1].split('_resul.xdr')[0] for xdr in match])
            print('Using the one provided in the input table...')

            match = [xdr for xdr in match if line in xdr]

        idldata = readsav(match[0])

        if check_best == True and filename != spec(id_star,SNR='bestHF').filename:
            print('\nWARNING: %s does not match with best spectrum available.'
            % filename)

        data_row = Table(idldata.d)
        data_row['Ref_file'] = filename
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
        # No values above 130 km/s were found in Simon-Diaz (2017)
        # we prevent the parameter to take values too large.
    data_rows['vmac_GF'][data_rows['vmac_GF'] > 130] = 130 
    #data_rows['vmac_GF_eUP'][data_rows['vmac_GF'] > 180] = np.nan
    #data_rows['vmac_GF_eDW'][data_rows['vmac_GF'] > 180] = np.nan
    data_rows['SNR_IB'] = [int(round(row)) for row in data_rows['SNR_IB']]

    # Saving the results:
    names = ['ID','Ref_file','vsini_FT','vmac_FT','vmac_FT_eDW','vmac_FT_eUP','vsini_GF',
             'vsini_GF_eDW','vsini_GF_eUP','vmac_GF','vmac_GF_eDW','vmac_GF_eUP',
             'line_IB','EW_IB','SNR_IB','QIB']

    output = Table(rows=data_rows[names], names=(names))

    full_path = maindir + 'tables/IB_results_new.' + format
    if format == 'ascii':
        format += '.fixed_width_two_line'

    output.write(full_path, format=format, overwrite=True)

    return 'DONE'
