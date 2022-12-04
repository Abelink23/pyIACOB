from spec import *
from scipy.io.idl import readsav


grids_dic = {
'all': ('Grids coverage','dodgerblue',0,
[[4.543,4.290,4.290,4.146,4.146,4.543,4.543],[2.391,2.391,3.092,3.092,4.391,4.391,2.391]]),
'nlte_10.1.6_SOLAR_expoclump_2019-10-24': ('BSgs_CNOSiMg','b',1,
[[4.190,4.477,4.477,4.190,4.190],[3.785,3.785,4.391,4.391,3.785]]),
'nlte_10.1.6_bdwarfs_SOLAR_2020-01-29': ('BDws_CNOSIMg','orange',2,
[[4.290,4.543,4.543,4.290,4.290],[2.391,2.391,3.889,3.889,2.391]]),
'nlte_10.4.7_OB.Sg_SOLAR_2021-01-23': ('OBSgs_hot_NOSi','g',3,
[[4.399,4.544,4.544,4.399,4.399],[3.488,3.488,4.386,4.386,3.488]]),
'nlte_10.4.7_late.bsgs_SOLAR_expoclump_NOSi.djl_2021-02-06': ('BSgs_cool_NOSi','r',4,
[[4.146,4.322,4.322,4.146,4.146],[3.092,3.092,4.391,4.391,3.092]]),
'astar2013_SOLAR_2_LMC_4_grid_2019-10-24_2019-10-24': ('ASgs_CNOMgSTiFe_Kurucz','purple',5,
[[3.900,4.114,4.114,3.900,3.900],[3.142,3.142,4.292,4.292,3.142]]),
'nlte_10.4.7_bsgs_SOLAR_expoclump_n12345o123c234mg2si234djl_v1_2021-05-05': ('BSg_CNOSiMg','DeepPink',6,
[[4.146,4.477,4.477,4.146,4.146],[3.392,3.392,4.386,4.386,3.392]]),
'nlte_10.4.7_bsgs_SOLAR_expoclump_n12345o123c234mg2si234djl_v1ehot_2022-01-19': ('O9BSg_CNOSiMg','turquoise',7,
[[4.146,4.543,4.543,4.146,4.146],[3.54,3.54,4.394,4.394,3.54]]),
'nlte_10.4.7_obgiants_SOLAR_noclump_n12345o123c234mg2si234djl_v1ehot_2022-02-21': ('O9BGs_CNOSiMg','lime',8,
[[4.204,4.543,4.543,4.204,4.204],[2.937,2.937,3.791,3.791,2.937]]),
}
''' IMPORTANT NOTE:
-   O9BSg_CNOSiMg has the upper limit of the logL (lgf) constraint <= 3.54 (1.85)
    also the upper limit of logQs <= -12.5
-   O9BGs_CNOSiMg has Teff > 16000K, logQs <= -12.5 and lgf <= 2.4

    These are cuts a posteriori in user_defined_model_set_fastwind_....pro in PRO_BASE_USER
    and are introduced manually in gen_gridlimits()
'''

def gen_gridlimits(models_dir=mauidir+'MODELS/'):

    '''
    Function to generate a fit table containing for each MAUI grid, the boundaries of
    each of the parameterm. The table is needed for the rest of the programs to work.

    Parameters
    ----------
    models_dir : str, optional
        Enter the directory where the model files are.

    Returns
    -------
    Nothing but the file is generated.
    '''

    param_dic = {
    'Teff' : ('Teff_UP','Teff_DW'),
    'logg' : ('logg_UP','logg_DW'),
    'lgf' : ('lgf_UP','lgf_DW'),
    'beta' : ('beta_UP','beta_DW'),
    'xf' : ('Micro_UP','Micro_DW'),
    'clf' : ('fcl_UP','fcl_DW'),
    'vclf' : ('vcl_UP','vcl_DW'),
    'logqs' : ('logQs_UP','logQs_DW'),
    'He' : ('He_UP','He_DW'),
    'C' : ('C_UP','C_DW'),
    'N' : ('N_UP','N_DW'),
    'O' : ('O_UP','O_DW'),
    'Mg' : ('Mg_UP','Mg_DW'),
    'Si' : ('Si_UP','Si_DW'),
    'S' : ('S_UP','S_DW'),
    'Ti' : ('Ti_UP','Ti_DW'),
    'Fe' : ('Fe_UP','Fe_DW')
    }

    param_other = {
    'micro' : 'xf',
    'HE' : 'He',
    'CARB' : 'C',
    'NIT' : 'N',
    'OXY' : 'O',
    'MAG' : 'Mg',
    'SUL' : 'S',
    'TIT' : 'Ti',
    'IRON' : 'Fe'
    }

    data_rows = []
    for file in os.listdir(models_dir):
        if not file.startswith('._') and file.endswith('.idl'):
            idldata = readsav(models_dir+file)

            data_row = []; data_row.extend([file.split('.idl')[0]])
            parameters = [i.decode() for i in idldata.param_labl]

            # To fix for the B8-As grid with different param labels:
            for i,par_name in zip(range(len(parameters)),parameters):
                if par_name in param_other: parameters[i] = param_other[par_name]

            for par_name in param_dic:

                if not par_name in parameters:
                    if not par_name == 'lgf':
                        data_row.extend([np.nan,np.nan])
                    continue

                idx = parameters.index(par_name)
                if not parameters[0] == 'Teff':
                    print('WARNING: lgf will be wrong!')

                if par_name == 'Teff':
                    data_row.extend(
                    [1e-4*idldata.param[idx].max(),
                    1e-4*idldata.param[idx].min()])

                elif par_name == 'logg':
                    # NOTE: with this we add loggf from logg to compare with the solution
                    # files from MAUI.
                    data_row.extend([idldata.param[idx].max(),idldata.param[idx].min()])

                    data_row.extend([
                    (idldata.param[1] -4*np.log10(idldata.param[0]) + 16).max(),
                    (idldata.param[1] -4*np.log10(idldata.param[0]) + 16).min(),
                    ])

                else:
                    data_row.extend([idldata.param[idx].max(),idldata.param[idx].min()])

            data_rows.append(tuple(data_row))

    names = ['Grid_name']
    for i in param_dic: names = names + [j for j in param_dic[i]]

    output = Table(rows=data_rows, names=(names))

    # Ad-hoc for one particular grid:
    if 'nlte_10.4.7_bsgs_SOLAR_expoclump_n12345o123c234mg2si234djl_v1ehot_2022-01-19' in output['Grid_name']:
        output['Teff_DW'][output['Grid_name']=='nlte_10.4.7_bsgs_SOLAR_expoclump_n12345o123c234mg2si234djl_v1ehot_2022-01-19'] = 1.4
        output['lgf_UP'][output['Grid_name']=='nlte_10.4.7_bsgs_SOLAR_expoclump_n12345o123c234mg2si234djl_v1ehot_2022-01-19'] = 1.85
        output['logQs_UP'][output['Grid_name']=='nlte_10.4.7_bsgs_SOLAR_expoclump_n12345o123c234mg2si234djl_v1ehot_2022-01-19'] = -12.5

    if 'nlte_10.4.7_obgiants_SOLAR_noclump_n12345o123c234mg2si234djl_v1ehot_2022-02-21' in output['Grid_name']:
        output['Teff_DW'][output['Grid_name']=='nlte_10.4.7_obgiants_SOLAR_noclump_n12345o123c234mg2si234djl_v1ehot_2022-02-21'] = 1.6
        output['lgf_UP'][output['Grid_name']=='nlte_10.4.7_obgiants_SOLAR_noclump_n12345o123c234mg2si234djl_v1ehot_2022-02-21'] = 2.40
        output['logQs_UP'][output['Grid_name']=='nlte_10.4.7_obgiants_SOLAR_noclump_n12345o123c234mg2si234djl_v1ehot_2022-02-21'] = -12.5

    output.write(maindir + 'tables/MAUI_grid_limits.fits', format='fits', overwrite=True)

    return 'DONE'


def maui_input(table, table_IB='IB_results.fits', table_REF=None,
    output_name='MAUI_input', RV0tol=200, ascii=False, txt=False):

    '''
    Function to generate the input table for MAUI given an input table with the
    target stars and quality flags for the line fittings.
    Optionally, the function allows to generate the input ascii spectra for MAUI
    subtracting the individual radial velocity.

    NOTE - Stars in the input table MUST also be in the other tables with the same name.

    Parameters
    ----------
    table : str
        Enter the input table contaning a column 'ID' with the name of the stars.

    table_IB : str
        Input table containing IACOB-broad parameters for the stars in 'table'.

    table_REF : str, optional
        Input table containing RVs, EWs, FWs and SNRs for the stars in 'table'.

    output_name : str, optional
        Enter the name for the output table. Default is 'MAUI_input'.

    RV0tol : int, optional
        Enter the input radial velocity tolerance for the radial velocity correction.

    ascii : boolean, optional
        If True, ascii files will be created for each of the input sources.
        Default is False.

    txt : boolean, optional
        If True, it will parse it to spec() in case that the input spectrum is in ascii.

    Returns
    -------
    Nothing but the MAUI input file is generated.
    '''

    if type(table) is type(Table()): pass # In case the input table is already a table
    else: table = findtable(table) # file where star names and quality flags are

    inpt_id = ''
    if len(table.colnames) > 1:
        print('If you want to use a specific column as the ID, type it here.')
        inpt_id = input('It must contain the full filename.\nOtherwise, press enter: \n')
        if inpt_id != '' and inpt_id in table.colnames:
            table['ID'] = table[inpt_id]
        else:
            print('Input column not found. Exiting...\n')
            return None

    table_IB = findtable(table_IB)

    if table_REF is not None:
        table_REF = findtable(table_REF)

    maui_txt = open(maindir + 'lists/%s.txt' % output_name,'a')
    maui_txt.write(
        "{:<40}".format('fullname')+"{:<48}".format('filename')+"{:<6}".format('vrad')+\
        "{:<6}".format('vsini')+"{:<7}".format('evsini')+"{:<5}".format('vmac')+\
        "{:<7}".format('evmac')+"{:<7}".format('R')+"{:<4}".format('SNR')+\
        ' ;#\n')

    quit = ''
    for row in table:

        do_ascii = ascii

        if quit == 'quit': break

        source = row['ID'].strip()

        if inpt_id != '':
            SNR = None
        else: 
            SNR = 'bestHF'

        star = spec(source, SNR=SNR, txt=txt)

        # Filer based on properties from the main table:
        if 'SB' in table.columns and (row['SB'] == 'SB2' or row['SB'] == 'SB3'):
            print('Skipping SB2 source %s' % source)
            continue

        # Checks if the star is in the IACOB-Broad table:
        match_IB = table_IB[table_IB['ID']==star.id_star]
        if len(match_IB) == 0:
            print('Info: Missing information in IB table for %s, skipping...\n' % source)
            continue

        # Checks if the star is in the IACOB-Broad table:
        if table_REF is not None:
            match_REF = table_REF[[i.strip()==star.id_star for i in table_REF['ID']]]

            if len(match_REF) == 0:
                print('Info: Missing information in RVEWFW table for %s, skipping...\n' % source)
                continue

            # Filter based on Si lines properties:
            #if ('QSiII' in table.columns and row['QSiII']<3) \
            #  or match_REF['EWSiII']<50 or np.isnan(match_REF['EWSiII']) \
            #  or match_REF['depSiII']<3/match_REF['SNR_B']:
            #    SiIIFG = 0
            #else: SiIIFG = 1

            #if ('QSiIII' in table.columns and row['QSiIII']<3) \
            #  or match_REF['EWSiIII1']<50 or np.isnan(match_REF['EWSiIII1']) \
            #  or match_REF['depSiIII1']<3/match_REF['SNR_B']:
            #    SiIIIFG = 0
            #else: SiIIIFG = 1

            #if SiIIIFG == 0 and SiIIFG == 0:
            #    print('Info: No SiIII/SiII found for %s, skipping...\n' % id)
            #    continue

        else:
            match_REF = Table([[star.snrcalc('B')]], names=('SNR_B',))

        if str(match_IB['Ref_file'][0].split('.fits')[0]) not in star.filename:
            print('WARNING: Different files from best SNR and from IB results for %s' % id)
            print(star.filename,' vs ',match_IB['Ref_file'][0])
            if ascii == True:
                do_file = input('Which ascii do you want to create 1 or 2: ')
                if int(do_file) == 1: star.filename = star.filename
                elif int(do_file) == 2: star.filename = match_IB['Ref_file'][0]

        # Checks if the ascii file for the star already exists to not repeat it
        if do_ascii == True and not search(star.filename.split('.')[0]+'_RV.ascii',\
        os.path.expanduser('~')+'/Documents/MAUI/SPECTRA/') is None:
            print('Info: ascii file for %s already exists' % star.filename.split('.')[0])
            do_ascii = False

        if do_ascii == True:

            gen_ascii(star.filename, db_table=table, spt='auto', rv_corr=True, RV0tol=200,
                cosmetic=True, cosmic=True, degrade=None, show_plot=True)

        # MAUI input last modifications:
        match_IB['evsini'] = abs(match_IB['vsini_GF_eDW'] + match_IB['vsini_GF_eUP'])/2

        try:
            if match_IB['vsini_GF'] < 5:
                match_IB['vsini_GF'] = 0
                match_IB['evsini'] = 0
        except:
            print(star.filename, match_IB)

        match_IB['evmac'] = abs(match_IB['vmac_GF_eDW'] + match_IB['vmac_GF_eUP'])/2

        # Based on comments from Simon-Diaz 
        if (match_IB['vmac_GF'][0] < 15) or (match_IB['vsini_GF'] > 130):
            match_IB['vmac_GF'] = 0
            match_IB['evmac'] = 0

        if match_REF['SNR_B'][0] > 200:
            match_REF['SNR_B'][0] = 200

        if star.filename.endswith('_RV.ascii'):
            apen = '.ascii'
        else:
            apen = '_RV.ascii'

        maui_txt.write(
            "{:<40}".format(star.filename.split('.')[0])\
            +"{:<48}".format(star.filename.split('.')[0] + apen) + "{:<6}".format('0.0d0')\
            +"{:<6}".format(str(int(round(match_IB['vsini_GF'][0],0))))\
            +"{:<7}".format(str(int(round(match_IB['evsini'][0]))))\
            +"{:<5}".format(str(int(round(match_IB['vmac_GF'][0]))))\
            +"{:<7}".format(str(int(round(match_IB['evmac'][0]))))\
            +"{:<7}".format(str(star.resolution)+'.')\
            +"{:<5}".format(str(int(round(match_REF['SNR_B'][0],0))))\
            +' ;#\n')

    maui_txt.close()

    if ascii == True:
        print('Remember to move the new ascii into "ASCII_ARCHIVE" folder')

    return 'DONE'


class solution_idl():
    def __init__(self, idlfile, solution='max', grids_table='MAUI_grid_limits.fits'):

        '''
        Parameters
        ----------
        idlfile : str
            Enter the input spectrum full path to the .idl file.

        solution : str, optional
            Choose between 'max'/'smooth' for taking the solution values from the maximum
            value of the distribution or from the maximum value after gaussian smoothing.
            Default is 'max'.

        grids_table : str, optional
            Provide the name of the table containing the information of the MAUI grids.
            Default is

        '''

        # Load the .idl file
        idldata = readsav(idlfile)
        #for i in idldata.keys():
        #    try: print(i,len(idldata[i]))
        #    except: print(i)

        # Identifiers
        self.filename = idldata.aa.label[0].decode() + '.fits'
        self.id_star = self.filename.split('_')[0]

        # Resolution
        self.resolution = int(idldata.obsdat.spectrum[0].reso_for_obs)

        # vsini and vmac used in the input for MAUI
        self.vsini = idldata.obsdat.spectrum[0].VSINI[0]
        self.vmac = idldata.obsdat.spectrum[0].MACRO[0]

        # Grid used
        self.gridname = idldata.modelgridname.decode()

        # Observed and synthetic spectrum
        self.obswave = idldata.xobs_out
        self.obsflux = idldata.yobs_out

        self.synwave = idldata.xx_mod
        try:
            self.synflux = idldata.spec_prim
        except:
            print(idlfile)

        self.synconv = idldata.sol_conv # This is flux. Combine with self.obswave

        # Delta lambda of spectrum
        self.dx = (idldata.xx_mod[-1]-idldata.xx_mod[0])/len(idldata.xx_mod)

        # Photometric information
        if idldata.phot_prim.dtype == 'int16':
            self.BC = self.B_V0 = np.nan
        else:
            self.BC = round(idldata.phot_prim[0], 3)
            self.B_V0 = round(idldata.phot_prim[2], 3)

        # Load table with the boundaries of each MAUI grid:
        try:
            grids = findtable(grids_table)
        except:
            print('Cannot find table with the MAUI grid boudaries (rquiered)')

        grid  = grids[grids['Grid_name'] == idldata.modelgridname.decode()]
        if len(grid) == 0:
            print('Cannot find the used grid in the grid database.')

        # QSA Parameters
        param_lst = ['Teff','logg','lgf','He','Micro','logQs','beta',
            'C','N','O','Mg','Si','S','Fe','Ti','fcl','vcl']

        # Fill the class with the parameter values:
        self.parameters = [j.decode() for j in idldata.solution.var_label[0]] # 11/13 param.
        if 'SOL_LOGG' in idldata.solution.dtype.names:
            self.parameters.append('logg') # assumes lgf is always

        for par_name in param_lst:

            if not par_name in self.parameters:
                setattr(self, 'l_'+par_name, '')
                [setattr(self, par_name+suffix, np.nan) for suffix in ['','_eUP','_eDW']]
                continue

            idx = self.parameters.index(par_name)

            if par_name == 'logg':
                if solution == 'max': # maximum of distribution without smoothing
                    sol_max = idldata.solution[0].sol_logg[0].map
                elif solution == 'smooth': # maximum of distribution without smoothing
                    sol_max = idldata.solution[0].sol_logg[0].map_smooth
                hpd_dw = idldata.solution[0].sol_logg[0].hpdint[0]
                hpd_up = idldata.solution[0].sol_logg[0].hpdint[1]

            else:
                if solution == 'max': # maximum of distribution without smoothing
                    sol_max = idldata.solution[0].sol_max[0][idx]
                elif solution == 'smooth': # maximum of distribution with smoothing
                    sol_max = idldata.solution[0].sol_max[1][idx]
                hpd_dw = idldata.solution[0].hpd_interval[0][idx]
                hpd_up = idldata.solution[0].hpd_interval[1][idx]

            # logQs is given as logQs-10:
            if par_name == 'logQs':
                sol_max -= 10
                hpd_up -= 10
                hpd_dw -= 10

            # Use 60% of the range to be considered as a degenerated case
            if abs(hpd_up-hpd_dw) > abs(grid[par_name+'_UP']-grid[par_name+'_DW'])[0]*0.6:
                label, err_dw, err_up = 'd', hpd_dw, hpd_up
                # Alternatively given by lower/upper limit of the grid for the param
                # label, err_dw, err_up = 'd', grid[par_name+'_DW'][0], grid[par_name+'_UP'][0]

            elif round(hpd_dw * 0.99, 3) < grid[par_name+'_DW'][0]:
                label, err_dw, err_up = '<', grid[par_name+'_DW'][0], hpd_up

            elif round(hpd_up * 1.01, 3) > grid[par_name+'_UP'][0]:
                label, err_dw, err_up = '>', hpd_dw, grid[par_name+'_UP'][0]

            else:
                label, err_dw, err_up = '=', abs(sol_max - hpd_dw), abs(sol_max - hpd_up)

            setattr(self, 'l_'+par_name, label)
            [setattr(self, par_name+suffix, value) for suffix,value in
                zip(['','_eUP','_eDW'],[round(sol_max,5),round(err_up,5),round(err_dw,5)])]

        # In most cases l_logg should be 'd' but is not due to 60% is not enough
        if 'Teff' in self.parameters and 'lgf' in self.parameters:
            if getattr(self, 'l_Teff') == 'd' or getattr(self, 'l_lgf') == 'd':
                self.l_logg = 'd'

        # Lines used in the analysis
        line_weights = []; line_names = []; line_windows = []
        for i in ['balmer','helium1','helium2','silicon']:
            lines = idldata.obsdat.spectrum[0][i][0]
            line_weights += [j for j in lines.weight[0]]
            line_names += [j.decode() for j in lines.line[0]]
            line_windows += [[j,k] for j,k in zip(lines.lmin[0],lines.lmax[0])]

        self.line_weights = line_weights
        self.line_names = line_names
        self.line_windows = line_windows


def maui_results(input_list, output_dir, check_best=False, last_only=False, FR=False,
    pdfplots=False, pdflines='diag', grids_table='MAUI_grid_limits.fits', grid_only=[],
    output_table=True, format_table='fits'):

    '''
    Function to generate a table with the results from MAUI given an input table
    containing the ID of the stars and filename to search in the MAUI-SOLUTION directory.

    Parameters
    ----------
    input_list : str
        Input star, list of stars, or table contaning the 'ID' or 'filename' to search.
        E.g. 'HD2905', 'HD2905,HD7902', 'table.txt/fits', '*' (to select all .idl in output_dir)

    output_dir : str
        Enter the path to the OUTPUT folder where the SOLUTION sub-folder containing the
        *solution*.idl files are. # Note: this could be your "mauidir" variable itself.

    check_best : boolean, optional
        True if each spectra from the input_table is checked against the best
        spectrum in the database. Default is False.

    last_only : boolean, optional
        True if only the last analysed .idl results want to be kept.

    FR : boolean, optional
        True if the analyses correspond to Fast Rotator stars for which whe MAUI windows
        are different. Default is False.

    pdfplots : boolean, optional
        If True, a pdf comparing the synthetic diagnostic lines with the original is made.

    pdflines : str, optional
        Choose between 'diag'/'all'/'def' to select the lines to be used in the pdf plots.
        Default is 'diag'.

    grids_table : str, obtional
        Name of the table containing the limits of the grids in MAUI.

    grid_only : list, optional
        List of grid names to limit the output to those results analysed with such grid.

    output_table : boolean, optional
        If True, the table with the output data will be created. Default is True.

    format_table : str, optional
        Enter the output format for the table: 'fits' (default), 'ascii' or 'csv'.

    Returns
    -------
    Nothing but the output table (+PDFs) with the MAUI results are generated.
    '''

    timenow = time.strftime('%Y%m%d_%H%M%S')

    # Create the input list from a table.
    # NOTE: A column named 'ID' or 'filename' is required. If 'filename' if provided, then
    # the specific file among other possible solution files for the same source will be used.
    if input_list.endswith(('.txt','.fits')):
        stars = findtable(input_list)
        try:
            stars = stars['filename']
        except:
            stars = stars['ID']

    # Create the input list from a list which has to be X.lst within list folder
    elif input_list.endswith('.lst'):
        stars = findlist(input_list)

    # Get the input list by the .idl files in the chosen SOLUTION folder
    elif input_list == '*':
        stars = [file.split('_sqexp_mat1_')[1].split('.idl')[0][0:-11] for file in \
            os.listdir(output_dir + 'SOLUTION/') if file.endswith('.idl')]
        stars = list(set(stars)) # To avoid duplicates if same file had several analyses

    # Create the input list from a string with one ID or IDs separated by coma
    else:
        stars = input_list.split(',')

    # Load the table with the information of the different grids
    grids = findtable(grids_table)

    # Set the progress bar
    bar = pb.ProgressBar(maxval=len(stars),
                         widgets=[pb.Bar('=','[',']'),' ',pb.Percentage()])
    bar.start()

    # Create pdf file to save the plots of the results
    if pdfplots == True:

        from matplotlib.backends.backend_pdf import PdfPages

        pdf_solution = PdfPages(maindir + 'plots/MAUI/MAUI_results_lines_%s.pdf' % timenow)
        pdf_makchain = PdfPages(maindir + 'plots/MAUI/MAUI_results_chain_%s.pdf' % timenow)

        plt.rcParams.update({
            'xtick.labelsize' : 6,
            'ytick.labelsize' : 6,
            'axes.titlesize' : 6,
            'axes.titlepad' : 3
            })

    # Parameters with errors
    param_err = ['Teff','logg','lgf','He','Micro','logQs','beta',
        'C','N','O','Mg','Si','S','Fe','Ti','fcl','vcl']

    # Iteration through the list of stars
    data_rows = []
    for i in range(len(stars)):
        name = stars[i]

        matches = []
        for file in os.listdir(output_dir + 'SOLUTION/'):
            if file.endswith('.idl'):

                # If input list contain the filename then this is used
                if '.ascii' in name and file.split('_sqexp_mat1_')[1][:-14] + 'RV.ascii' == name:
                    matches.append(output_dir + 'SOLUTION/' + file)

                # while if the ID of the star is used, it is still contained as _ID_ in the file
                elif '_' + name + '_' in file.split('_sqexp_mat1')[1][:-14]:
                    matches.append(output_dir + 'SOLUTION/' + file)

        if len(matches) == 0:
            print('\nWARNING: No .idl file found for %s. Continuing...' % name)
            continue

        if len(matches) > 1 and last_only == True:
            matches = [sorted(matches, key=lambda x: int(x[-14:-4].replace('-','')), reverse=True)[0]]

        for match in matches:

            # Load the idl class for the file
            try:
                star = solution_idl(match)
            except:
                print('ERROR: Problem loading file: %s, skipping...' % match)
                continue

            # Skip the used grid if not selected from the grid_only keyword
            if grid_only != [] and not star.gridname in grid_only: continue

            # Find the best SNR spectra in the DB
            if check_best == True or pdfplots == True:
                best_SNR = spec(star.id_star, SNR='bestHF')

            # Check if the input file matches with the best SNR spectra available
            if check_best == True and star.filename != best_SNR.filename:
                print('\nWARNING: %s does not match with best spectrum available.' % star.filename)

            # Generate the table with the results of each of the parameters (basic and with errors)
            data_row = [star.id_star,star.filename,star.gridname,star.BC,star.B_V0]

            [data_row.extend([
                getattr(star, 'l_'+par_name),
                getattr(star, par_name),
                getattr(star, par_name+'_eUP'),
                getattr(star, par_name+'_eDW'),
                ]) for par_name in param_err]

            # Append the raw to the table
            data_rows.append(tuple(data_row))

            if pdfplots == True:

                # PLOT OF SPECTRAL LINES OUT OF THE IDL SOLUTION FILES
                if pdflines == 'diag':
                    mask = [i == 1. for i in star.line_weights]

                    line_names = np.asarray(star.line_names)[mask].tolist()
                    lines_lwl,lines_rwl = np.asarray(star.line_windows).T
                    lines_lwl = np.asarray(lines_lwl)[mask].tolist()
                    lines_rwl = np.asarray(lines_rwl)[mask].tolist()
                    line_colors = ['g']*len(line_names)

                elif pdflines == 'all':
                    line_names = star.line_names
                    lines_lwl,lines_rwl = np.asarray(star.line_windows).T
                    #lines_lwl = np.asarray(lines_lwl).tolist()
                    #lines_rwl = np.asarray(lines_rwl).tolist()
                    line_colors = ['g' if i == 1 else 'r' for i in star.line_weights]

                else:
                    line_names = [
                        r'H$_{\delta}$',r'H$_{\gamma}$',r'H$_{\beta}$',r'H$_{\alpha}$',
                        'HeI 5016','HeI 5876','HeII 4542','HeII 5412',
                        'SiIV 4116','SiIII 4552','SiIII 4568/75','SiII 6347',
                        'NII 3995','CII 4267','OII 4662','MgII 4881',
                        ]

                    lines_lamb = [
                        4101.735,4340.463,4861.325,6562.8,
                        5015.678,5875.62,4541.591,5411.52,
                        4116.103,4552.622,4571.2985,4130.89,
                        3994.997,4267.183,4661.632,4481.126]

                    lines_lwl = [i-10 for i in lines_lamb]
                    lines_rwl = [i+10 for i in lines_lamb]

                    line_colors = ['g']*len(line_names)

                nrows, ncols = even_plot(len(line_names))

                fig, ax = plt.subplots(nrows, ncols, figsize=(13,8))

                fig_title = ''
                for par_name in param_err:
                    val = getattr(star, par_name)
                    if np.isnan(val) == True: continue

                    label = getattr(star, 'l_'+par_name)
                    if label == 'd': label = '=d.'
                    if par_name == 'Teff':
                        val = val*1e4
                        Teff = val

                    fig_title += par_name+label+str(round(val,2))+'  '

                fig.suptitle(star.id_star + ' -- ' + match.split('emulated_solution_mcmc_sqexp_mat1_')[-1]
                    + ' -- ' + star.gridname + '\n' + fig_title, fontsize=8)

                axs = ax.flatten()
                for ax_i,line_lwl,line_rwl,line_name,c in zip(axs,lines_lwl,lines_rwl,line_names,line_colors):
                    #mask = (star.synwave > line_lwl) & (star.synwave < line_rwl)
                    #ax_i.plot(star.synwave[mask], star.synflux[mask], color='gray', lw=.3)

                    mask = (star.obswave > line_lwl) & (star.obswave < line_rwl)
                    ax_i.plot(star.obswave[mask], star.obsflux[mask], color='k', lw=.7)
                    if ax_i.get_ylim()[0] > 0.9:
                        ax_i.set_ylim(bottom=0.895)

                    if FR == False:
                        blocked = [None if (
                            (i >= 6521. and i <= 6532.) or # Halpha, exclude CII lines in the red wing
                            (i > 6575.75 and i < 6585.29) or # Halpha
                            (i > 4344.42 and i < 4355.83) or # Hgamma, exclude OII lines in the red wing
                            (i >= 4083.20 and i <= 4085.50) or # Hdelta, exclude metal lines
                            (i >  4086.29 and i <  4090.53) or
                            (i >  4092.27 and i <  4093.99) or
                            (i >= 4096.5  and i <= 4098.) or
                            (i >= 3953.09 and i <= 3955.05) or # Hepsil, exclude metal lines
                            (i >  3960.70 and i <  3962.24) or
                            (i >  3963.38 and i <  3965.76) or
                            (i >= 3967.00 and i <= 3968.50) or
                            (i >= 3972.53 and i <= 3974.10) or
                            (i >= 4923.5 and i <= 4926.0) or # HeI 4922
                            (i >= 5017.5 and i <= 5019.5) or # HeI 5015
                            (i >= 4544.0 and i <= 4546.0) or # HeII 4541, exclude AlIII lines
                            (i >  4478.87 and i <  4480.20) or # MgII, exclude the AlIII blend
                            (i >= 4128.71 and i <= 4130.10) or # SiIII 4130
                            (i >= 4131.40 and i <= 4134.08) or
                            (i >= 4547.60 and i <= 4550.99) or # SiIII 4552
                            (i >= 4554.20 and i <= 4558.96) or
                            (i >= 4563.08 and i <= 4565.97) or # SiIII 4567
                            (i >= 4569.925 and i <= 4571.35) or
                            (i >= 4571.30 and i <= 4573.01) or # SiIII 4575
                            (i >= 4575.78 and i <= 4579.32) or
                            (i >= 4114.71 and i <= 4114.81) or # SiIV 4116
                            (i >= 4117.7 and i <= 4125.10)
                            ) else np.asarray(ax_i.get_ylim()).mean() for i in star.obswave[mask]]

                    elif FR == True:
                        blocked = [None if (
                            (i >= 6521. and i <= 6532.) or # Halpha, exclude CII lines in the red wing
                            (i > 6575.75 and i < 6585.29) or # Halpha
                            (i > 4344.42 and i < 4355.83) or # Hgamma, exclude OII lines in the red wing
                            (i >= 4083.20 and i <= 4090.53) or # Hdelta, exclude metal lines
                            (i >  4092.27 and i <  4093.99) or
                            (i >= 4096.5  and i <= 4098.) or
                            (i >= 3953.09 and i <= 3955.05) or # Hepsil, exclude metal lines
                            (i >  3960.70 and i <  3962.24) or
                            (i >  3963.38 and i <  3965.76) or
                            (i >= 3967.00 and i <= 3968.50) or
                            (i >= 3972.53 and i <= 3974.10) or
                            (i >= 4923.5 and i <= 4926.0) or # HeI 4922
                            (i >= 5017.5 and i <= 5019.5) or # HeI 5015
                            (i >= 4544.0 and i <= 4546.0) or # HeII 4541, exclude AlIII lines
                            (i >= 4128.71 and i <= 4130.10) or # SiIII 4130
                            (i >= 4131.40 and i <= 4134.08) or
                            (i >= 4110.90 and i <= 4111.00) or # SiIV 4116
                            (i >= 4117.7 and i <= 4126.10)
                            ) else np.asarray(ax_i.get_ylim()).mean() for i in star.obswave[mask]]

                    # This is a visual trick to put the synthetic spectra where the normalization
                    # feature is putting it as it is not stored but is the way to see its effect
                    weight = [0 if i is None else 1 for i in blocked]
                    scale = np.sum(star.obsflux[mask]*star.synconv[mask]*weight)\
                                /np.sum(star.synconv[mask]*star.synconv[mask]*weight)
                    star.synconv_scaled = scale * star.synconv[mask]

                    ax_i.plot(star.obswave[mask], star.synconv_scaled, color=c, ls='--', lw=1)

                    ax_i.plot(star.obswave[mask], blocked, c='cyan', lw=.5, alpha=0.5)

                    ax_i.set_title(line_name)
                    ax_i.tick_params(direction='in', top='on', right='on')

                [fig.delaxes(axs[i]) for i in np.arange(len(line_names), len(axs), 1)]

                fig.tight_layout()
                fig.subplots_adjust(wspace=0.14, hspace=0.3)

                pdf_solution.savefig(fig); plt.close(fig)

                # PLOT OF PROBABILLITY DISTRIBUTIONS OUT OF THE IDL MARKOV-CHAIN FILES
                mcmc_file = match.split('emulated_solution_')[1]
                try:
                    mcmc_data = readsav(output_dir + 'MARKOV_CHAIN/' + mcmc_file)
                except:
                    print('%s associated mcmc file not found in MARKOV_CHAIN/ folder' % star.id_star)
                    continue

                parameters = [var.decode() for var in mcmc_data.varname]

                nrows, ncols = even_plot(len(parameters))

                fig, ax = plt.subplots(nrows, ncols, figsize=(13,8))
                fig.subplots_adjust(wspace=.5, hspace=.5)

                fig.suptitle(star.id_star + ' -- ' + match.split('emulated_solution_mcmc_sqexp_mat1_')[-1]
                    + ' -- ' + star.gridname + '\n' + fig_title, fontsize=8)

                axs = ax.flatten()

                for j in range(len(parameters)):

                    xin = mcmc_data.chain_final.T[j]*(mcmc_data.xmax[j] - mcmc_data.xmin[j]) + mcmc_data.xmin[j]

                    # Replace lgf by logg
                    #if 'lgf' in parameters and parameters.index('lgf') == i:
                    #    idx = parameters.index('Teff')
                    #    teff = mcmc_data.chain_final.T[idx]*(mcmc_data.xmax[idx] - mcmc_data.xmin[idx]) + mcmc_data.xmin[idx]
                    #    xin = xin + 4*np.log10(teff)
                    #    parameters[i] = 'logg'

                    iqr = np.quantile(xin, q=[.25, .75])
                    fd_bin = 2*np.diff(iqr)/(len(xin)**(0.3))

                    weights = np.ones_like(xin)/float(len(xin))
                    axs[j].hist(xin, bins=np.arange(min(xin), max(xin) + fd_bin[0], fd_bin[0]),
                        weights=weights, histtype='stepfilled', fc='gray', ec='g', lw=2, alpha=0.6)
                    ylim = axs[j].get_ylim()[1]
                    axs[j].plot([np.median(xin),np.median(xin)], [0,ylim], '--', c='orange')
                    axs[j].plot([iqr[0],iqr[0]], [0,ylim], '-.', c='orange')
                    axs[j].plot([iqr[1],iqr[1]], [0,ylim], '-.', c='orange')
                    axs[j].set_title(parameters[j])
                    axs[j].tick_params(direction='in', top='on')

                [fig.delaxes(axs[i]) for i in np.arange(len(parameters), len(axs), 1)]

                fig.tight_layout()

                pdf_makchain.savefig(fig); plt.close(fig)

        bar.update(i)

    if pdfplots == True:
        pdf_solution.close()
        pdf_makchain.close()

    bar.finish()

    # Saving the results:
    # - Add names for the basic columns
    names = ['ID','filename','Grid_name','BC','BV_0']

    # - Add names for the parameter column names
    for i in range(len(param_err)):
        names += ['l_'+param_err[i],param_err[i],param_err[i]+'_eUP',param_err[i]+'_eDW']

    output = Table(rows=data_rows, names=(names))

    full_path = (maindir + 'tables/MAUI_results_%s.' + format_table) % timenow

    if format_table == 'ascii':
        format_table += '.fixed_width_two_line'

    if output_table == True:
        output.write(full_path, format=format_table, overwrite=True)

    return 'DONE'


def gen_stars_in_grids(input_table, table_results):

    '''
    Function to generate MAUI-input txt lists containing the stars that lie within the
    boundaries of each MAUI grid.

    Parameters
    ----------
    input_table : str
        Table containing the input stars under an 'ID' column name.

    table_results : str
        Table containing the MAUI results for the stars in the previous table.

    Returns
    -------
    Nothing, but the output lists are created.
    '''

    import matplotlib.path as mpath

    table = findtable(input_table)
    results = findtable(table_results)

    log_Teff = np.asarray(4+np.log10(results['Teff']))
    log_LLsol = np.asarray(5.39-results['lgf'])

    points = np.column_stack([log_Teff,log_LLsol])

    # Change ...grids_dic][1:] for individual grids
    for name,_,_,box in [grids_dic[i] for i in grids_dic][1:]:
        verts = np.array([box[0],box[1]]).T
        path = mpath.Path(verts)
        inout = path.contains_points(points)
        log_Teff_in,log_LLsol_in = points[inout].T

        results_in = results[path.contains_points(points)]['ID']
        table_red = table[[i['ID'] in results_in for i in table]]

        maui_input(table=table_red, output_name=name, ascii_0=False)


def phot_table(input_table):

    '''
    Function to generate an ascii table which contains the input for an IDL program
    that computes the photometric parameters for stars analysed with MAUI

    Parameters
    ----------
    input_table : str
        Enter the name of the table (located in table/ folder) that contains the results
        of the stars from MAUI.

    Returns
    -------
    Nothing but the ascii table is generated.
    '''

    table = findtable(input_table)[0:10]

    cols = ['Teff','lgf','He','Micro','logQs','beta','C','N','O','Mg','Si','fcl','vcl','Grid_name']
    if any([i not in table.colnames for i in cols]):
        print('WARNING: missing column names. Exiting...\n')
        return None

    table['Teff'] = [('%.3f') % i +'d0' for i in table['Teff']]
    table['logf'] = [('%.3f') % i +'d0' for i in table['lgf']]
    table['He'] = [('%.2f') % i for i in table['He']]
    table['mic'] = [('%.0f') % i+'.' for i in table['Micro']]
    table['logQ+10'] = [('%.2f') % (i+10) +'d0' for i in table['logQs']]
    table['beta'] = [('%.2f') % i for i in table['beta']]
    for i in ['C','N','O','Mg','Si']:
        table[i] = [('%.2f') % i +'d0' for i in table[i]]
    table['fcl'] = [('%.2f') % i if i != 'nan' else 1.00 for i in table['fcl']]
    table['vcl'] = [('%.0f') % i+'.' if i != 'nan' else 100 for i in table['vcl']]

    table = table['ID','Teff','logf','He','mic','logQ+10','beta','C','N','O','Mg','Si','fcl','vcl','Grid_name']
    table.write(maindir+'tables/%s_phot.txt' % input_table.split('.')[0], format='ascii.fixed_width_two_line')


def gen_synthetic(output_dir, save_dir='server', lwl=3900, rwl=5080):

    '''
    Function to generate the .dat files of the synthetic spectra generated by
    MAUI in all the .idl files found in the SOLUTION folder.

    Parameters
    ----------
    output_dir : str
        Enter the path to the OUTPUT folder where the SOLUTION sub-folder containing the
        *solution*.idl files are. # Note: this could be your "mauidir" variable itself.

    save_dir : str, optional
        Enter the directory where to save all the output spectra in ascii format.

    lwl : float, optional
        Sets the start wavelenght of the output spectrum.

    rwl : float, optional
        Sets the end wavelenght of the output spectrum.

    Returns
    -------
    Nothing but the ascii .dat files are generated.
    '''

    if save_dir == 'local':
        save_dir = datadir + 'ASCII/Synthetic_MAUI/'

    elif save_dir == 'server':
        save_dir = '/net/nas/proyectos/hots/masblue/obs_iac/spec_opt/IACOB_DB/ASCII/SYNTHETIC/'

    for file in os.listdir(output_dir + 'SOLUTION/'):
        if not file.startswith('._') and file.endswith('.idl'):
            star = solution_idl(output_dir + 'SOLUTION/' + file)

            star_db = spec(star.id_star, SNR='bestHF')

            if star.filename != star_db.filename:
                print('\nWARNING: %s does not match with best spectrum available.'
                % star.filename)

            else:
                #star.filename = star.filename.replace(str(star.resolution),'85000')
                new_star = '%s_red%i.dat' % (star.filename,grids_dic[star.gridname][2])
                np.savetxt(save_dir + new_star, np.c_[star.synwave,star.synflux],
                    fmt=('%.4f','%.6f'))

                star_idl = spec(new_star, txt=True)
                star_idl.txtwaveflux(lwl, rwl)
                #plt.plot(star_idl.wave, star_idl.flux, 'r', lw=.3) # plot to check
                star_idl.degrade(profile='rotmac', vsini=star.vsini, vmac=star.vmac)
                #plt.plot(star_idl.wave, star_idl.flux, 'g', lw=.3) # plot to check

                np.savetxt(save_dir + new_star, np.c_[star_idl.wave,star_idl.flux],
                    fmt=('%.4f','%.6f'))

def even_plot(n):
    '''
    Function to determine the number of rows, columns needed to have a square plot.

    Parameters
    ----------
    n : int
        Number of elements to be splited.

    Returns
    -------
        Number of rows and columns
    '''

    nrows, ncols = int(np.ceil(np.sqrt(n))), round(n/np.ceil(np.sqrt(n))+0.4)

    return nrows,ncols
