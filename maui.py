from turtle import color

from spec_posproc import *
from scipy.io import readsav

from hardcoded_maui import *

def gen_gridlimits(models_dir=mauidir+'MODELS/'):

    '''
    Function to generate a fit table containing for each MAUI grid, the boundaries of
    each of the parameters. The table is needed for the rest of the programs to work.

    Parameters
    ----------
    models_dir : str, optional
        Enter the directory where the MAUI model files are (e.g. MAUI/MODELS/).

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
        'fcl' : ('fcl_UP','fcl_DW'),
        'vcl' : ('vcl_UP','vcl_DW'),
        'logQs' : ('logQs_UP','logQs_DW'),
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
        'vclf': 'vcl',
        'clf' : 'fcl',
        'logqs' : 'logQs',
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
        if not file.startswith('.') and file.endswith('.idl'):
            soldata = readsav(models_dir+file)

            data_row = []; data_row.extend([file.split('.idl')[0]])
            parameters = [i.decode() for i in soldata.param_labl]

            # Obtain indices of the Teff and logg parameters
            idx_Teff = parameters.index('Teff')
            idx_logg = parameters.index('logg')

            # To fix for the B8-As grid with different param labels:
            for i,par_name in zip(range(len(parameters)),parameters):
                if par_name in param_other:
                    parameters[i] = param_other[par_name]

            for par_name in param_dic:

                if not par_name in parameters:
                    if not par_name == 'lgf':
                        data_row.extend([np.nan,np.nan])
                    continue

                idx = parameters.index(par_name)
                if not parameters[0] == 'Teff':
                    msg.warn('lgf will be wrong!')

                if par_name == 'Teff':
                    data_row.extend(
                    [1e-4*soldata.param[idx].max(),
                    1e-4*soldata.param[idx].min()])

                elif par_name == 'logg':
                    # NOTE: with this we add loggf from logg to compare with the solution
                    # files from MAUI.
                    data_row.extend([soldata.param[idx].max(),soldata.param[idx].min()])

                    data_row.extend([
                    (soldata.param[idx_logg] -4*np.log10(soldata.param[idx_Teff]) + 16).max(),
                    (soldata.param[idx_logg] -4*np.log10(soldata.param[idx_Teff]) + 16).min(),
                    ])

                else:
                    data_row.extend([soldata.param[idx].max(),soldata.param[idx].min()])

            # if the row has floats, round it to 3 decimals:
            data_row = [round(i,3) if type(i) in [float,np.float32,np.float64] else i for i in data_row]

            data_rows.append(tuple(data_row))

    names = ['Grid_name']
    for i in param_dic: names = names + [j for j in param_dic[i]]

    output = Table(rows=data_rows, names=(names))

    # Ad-hoc changes in grids:
    for grid_key, limits in ad_hoc_limits.items():
        if grid_key in output['Grid_name']:
            for param, value in limits.items():
                output[param][output['Grid_name']==grid_key] = value

    output.write(maindir + 'tables/MAUI_grid_limits.fits', format='fits', overwrite=True)

    return 'DONE'


def maui_input(table, table_IB='IB_results.fits', output_name='MAUI_input',
    RV0tol=200, ascii=False, spectra_path=mauidir+'/SPECTRA/', orig='IACOB'):

    '''
    Function to generate the input table for MAUI given an input table with the
    target stars and quality flags for the line fittings.
    Optionally, the function allows to generate the input ascii spectra for MAUI
    subtracting the individual radial velocity.

    NOTE - Stars in the input table MUST also be in the other tables with the same name.

    Parameters
    ----------
    table : str
        Enter the input table containing a column 'ID' with the name of the stars.

    table_IB : str
        Input table containing IACOB-broad parameters for the stars in 'table'.

    output_name : str, optional
        Enter the name for the output table. Default is 'MAUI_input'.

    RV0tol : int, optional
        Enter the input radial velocity tolerance for the radial velocity correction.

    ascii : boolean, optional
        If True, ascii files will be created for each of the input sources.
        Default is False.

    spectra_path : str, optional
        Enter the path to the spectra inside the maui folder defined in paths.txt.
        Note: it only works if the keyword 'ascii' is True.

    orig : str, optional
        See spec() function for more information. Default is 'IACOB'.

    Returns
    -------
    Nothing but the MAUI input file is generated.
    '''

    if type(table) is type(Table()):
        print('Input table is already a table. Disabling generation of ascii files (not implemented yet)...\n')
        ascii = False
        pass # In case the input table is already a table
    else:
        db_table = table
        table = findtable(table) # file where star names and quality flags are

    inpt_id = ''
    if len(table.colnames) > 1:
        print('If you want to use a specific column for the ID / reference file, type it here.')
        inpt_id = input('Otherwise, press enter: \n')
        if inpt_id != '' and inpt_id in table.colnames:
            table['ID'] = table[inpt_id]
        else:
            print('Input column not found. Exiting...\n')
            return None

    table_IB = findtable(table_IB)

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
            snr = None
        else:
            snr = 'bestHF'

        star = spec(source, snr=snr, orig=orig)

        # Filer based on properties from the main table:
        if 'SB' in table.columns and (row['SB'] == 'SB2' or row['SB'] == 'SB3'):
            msg.info('Skipping SB2 source %s' % source)
            continue

        # Checks if the star is in the IACOB-Broad table:
        match_IB = table_IB[table_IB['ID']==star.id_star]
        if len(match_IB) == 0:
            msg.warn('Missing information in IB table for %s, skipping...\n' % source)
            continue

        # Obtains the SNR in the 4000-5000 range
        SNR_B = star.snrcalc('B')

        if str(match_IB['Ref_file'][0].split('.fits')[0]) not in star.filename:
            msg.warn('Different files from best SNR and from IB results for %s' % id)
            msg.warn(star.filename,' vs ',match_IB['Ref_file'][0])
            if ascii == True:
                do_file = input('Which ascii do you want to create 1 or 2: ')
                if int(do_file) == 1: star.filename = star.filename
                elif int(do_file) == 2: star.filename = match_IB['Ref_file'][0]

        # Checks if the ascii file for the star already exists to not repeat it
        if do_ascii == True and not search(star.filename.split('.')[0]+'_RV.ascii',\
            spectra_path) is None:
            msg.info('ascii file for %s already exists' % star.filename.split('.')[0])
            do_ascii = False

        if do_ascii == True:

            gen_ascii(star.filename, db_table=db_table, spt='auto', rv_corr=True, RV0tol=RV0tol,
                export_rv=True, cosmetic=True, cosmic=True, degrade=None, show_plot=True)
            plt.close('all')

        # MAUI input final modifications:
        if 'vsini_GF_eDW' in match_IB.columns and 'vsini_GF_eUP' in match_IB.columns:
            match_IB['evsini'] = abs(match_IB['vsini_GF_eDW'] + match_IB['vsini_GF_eUP'])/2
        else:
            match_IB['evsini'] = 0.10*match_IB['vsini_GF']

        if match_IB['vsini_GF'] < 5:
            msg.info('vsini_GF <5 km/s for %s. Setting value and error to 0.' % star.filename)
            match_IB['vsini_GF'] = 0
            match_IB['evsini'] = 0

        if match_IB['vsini_GF'] < 40:
            msg.info('vsini<40km/s for %s -> Following Simon-Diaz+14 this value is likely an upper limit.' % star.filename)

        if 'vmac_GF_eDW' in match_IB.columns and 'vmac_GF_eUP' in match_IB.columns:
            match_IB['evmac'] = abs(match_IB['vmac_GF_eDW'] + match_IB['vmac_GF_eUP'])/2
        else:
            match_IB['evmac'] = 0.10*match_IB['vmac_GF']

        # From Simon-Diaz+14 (vmac <15 for Gs and vmac <40 for SGs are not reliable -> we measure the micro)
        if (match_IB['vmac_GF'][0] < 15) or (match_IB['vsini_GF'] > 130):
            msg.info('vmac <15 or vsini >130 km/s for %s. Setting vmac and error to 0.' % star.filename)
            match_IB['vmac_GF'] = 0
            match_IB['evmac'] = 0

        if SNR_B > 200:
            SNR_B = 200

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
            +"{:<5}".format(str(int(round(SNR_B))))\
            +' ;#\n')

    maui_txt.close()

    if ascii == True:
        print('Remember to move the new ascii into "ASCII_ARCHIVE" folder')

    return 'DONE'


def split_MAUI_input(path_to_file, number_per_file=10):

    '''
    Function to split the MAUI input file into several files with a maximum number
    of entries per file.

    Parameters
    ----------
    path_to_file : str
        Enter the full path to the MAUI input file.

    number_per_file : int, optional
        Enter the maximum number of entries per file. Default is 10.

    Returns
    -------
    Nothing but the MAUI input files are generated.
    '''

    filename = path_to_file.split('/')[-1][:-4]
    path = path_to_file.split(filename)[0]

    with open(path_to_file, 'r') as file:
        file = file.read().splitlines()
        file = [i for i in file if not i.startswith('fullname')]

    n = 0; m = 1; new_file = []
    for line in file:

        new_file.append(line)
        n = n + 1

        # if the last file has less than number_per_file entries, the last lines will be added to the last new file
        if len(file) % number_per_file != 0 and m == len(file) // number_per_file:
            number_per_file = number_per_file + len(file) % number_per_file

        if n == number_per_file:
            with open(path+filename+'_%i.txt' % m, 'w') as f:
                f.write('%s\n' % str(n))
                f.write('fullname                                filename                                        vrad  vsini evsini vmac evmac  R      SNR ;# \n')
                for item in new_file:
                    f.write("%s\n" % item)
            new_file = []
            n = 0; m += 1


class solution_maui():
    def __init__(self, solfile, mcmcfile=None, solution='smooth'):

        '''
        Parameters
        ----------
        solfile : str
            Enter the input spectrum full path to the solution*.idl file.

        mcmcfile : str, optional
            Enter the input spectrum full path to the mcmc*.idl file.

        solution : str, optional
            Choose between 'max'/'smooth' for taking the solution values from the maximum
            value of the distribution or from the maximum value after gaussian smoothing.
            Default is 'smooth'.

        '''

        # Load the solution*.idl file
        soldata = readsav(solfile)
        #for i in soldata.keys():
        #    try: print(i,len(soldata[i]))
        #    except: print(i)

        # Identifiers
        self.filename = soldata.aa.label[0].decode() + '.fits'
        self.id_star = self.filename.split('_')[0]

        # Resolution
        self.resolution = soldata.obsdat.spectrum[0].reso_for_obs
        self.resolution = self.resolution[0] if type(self.resolution) in [np.ndarray, list] else self.resolution

        # vsini and vmac used in the input for MAUI
        self.vsini = soldata.obsdat.spectrum[0].VSINI[0]
        self.vmac = soldata.obsdat.spectrum[0].MACRO[0]

        # Grid used
        self.gridname = soldata.modelgridname.decode()

        # Observed and synthetic spectrum
        self.obswave = soldata.xobs_out
        self.obsflux = soldata.yobs_out

        # Synthetic wavelengths for spec_prim
        self.synwave = soldata.xx_mod
        # Synthetic not convolved flux of the solution
        try:
            self.synflux = soldata.spec_prim
        except:
            print(solfile)

        # The convolved flux of the solution. Combine with self.obswave
        self.synconv = soldata.sol_conv

        # Delta lambda of spectrum
        self.dx = (soldata.xx_mod[-1]-soldata.xx_mod[0])/len(soldata.xx_mod)

        # Photometric information
        if soldata.phot_prim.dtype in ['int16','int32','float32','float64']:
            self.BC = self.B_V0 = np.nan
        else:
            self.BC = round(soldata.phot_prim[0], 3)
            self.B_V0 = round(soldata.phot_prim[2], 3)

        # Load the mcmc*.idl file
        if mcmcfile is not None:
            mcmcdata = readsav(mcmcfile)
        else:
            mcmcdata = None

        # Fill the class with the parameter values:
        self.parameters = [j.decode() for j in soldata.solution.var_label[0]] # 11/13 param.
        if 'SOL_LOGG' in soldata.solution.dtype.names:
            self.parameters.append('logg') # assumes lgf is always

        # Take hardcoded QSA Parameters
        for par_name in dic_maui_param:

            if not par_name in self.parameters:
                setattr(self, 'l_'+par_name, '')
                [setattr(self, par_name+suffix, np.nan) for suffix in ['','_eUP','_eDW']]
                continue

            idx = self.parameters.index(par_name)

            if par_name == 'logg':
                sol_max = soldata.solution[0].sol_logg[0].map
                sol_smooth = soldata.solution[0].sol_logg[0].map_smooth
                if solution == 'max': # maximum of distribution without smoothing
                    sol_final = sol_max
                elif solution == 'smooth': # maximum of distribution without smoothing
                    sol_final = sol_smooth
                hpd_dw = soldata.solution[0].sol_logg[0].hpdint[0]
                hpd_up = soldata.solution[0].sol_logg[0].hpdint[1]

                idx_tef = self.parameters.index('Teff')
                idx_lgf = self.parameters.index('lgf')

                if mcmcdata is not None:
                    chain_teff = mcmcdata.chain_final.T[idx_tef]*(mcmcdata.xmax[idx_tef] - mcmcdata.xmin[idx_tef]) + mcmcdata.xmin[idx_tef]
                    chain_lgf = mcmcdata.chain_final.T[idx_lgf]*(mcmcdata.xmax[idx_lgf] - mcmcdata.xmin[idx_lgf]) + mcmcdata.xmin[idx_lgf]

                    chain = [(chain_lgf +4*np.log10(chain_teff*1e4) - 16).min(), (chain_lgf +4*np.log10(chain_teff*1e4) - 16).max()]

            else:
                sol_max = soldata.solution[0].sol_max[0][idx]
                sol_smooth = soldata.solution[0].sol_max[1][idx]
                if solution == 'max': # maximum of distribution without smoothing
                    sol_final = sol_max
                elif solution == 'smooth': # maximum of distribution with smoothing
                    sol_final = sol_smooth
                hpd_dw = soldata.solution[0].hpd_interval[0][idx]
                hpd_up = soldata.solution[0].hpd_interval[1][idx]

                if mcmcdata is not None:
                    chain = [mcmcdata.xmin[idx], mcmcdata.xmax[idx]]

            if abs(sol_max - sol_smooth) > 0.10*abs(sol_max):
                msg.warn('max vs smooth values differ by more than 10%% or parameter %s in %s.' % (par_name, self.filename))

            # logQs is given as logQs-10:
            if par_name == 'logQs':
                sol_final -= 10
                hpd_up -= 10
                hpd_dw -= 10
                chain = [i-10 for i in chain] if mcmcdata is not None else None

            delta = 0.01
            if np.sign(sol_max) == -1:
                delta = -delta

            if mcmcdata is not None:
                # Use 70% of the range to be considered as a degenerated case
                if abs(hpd_up-hpd_dw) > abs(max(chain)-min(chain))*0.70:
                    label, err_dw, err_up = 'd', hpd_dw, hpd_up
                # If the solution is close to the limits of the chain (i.e., the grid),
                # we consider the difference with those limits as the error bars and we
                # label it the parameter as upper or lower limit depending on the case,
                elif round(hpd_dw*(1-delta), 3) < min(chain):
                    label, err_dw, err_up = '<', abs(sol_final - min(chain)), abs(sol_final - hpd_up)

                elif round(hpd_up*(1+delta), 3) > max(chain):
                    label, err_dw, err_up = '>', abs(sol_final - hpd_dw), abs(sol_final - max(chain))

                else:
                    label, err_dw, err_up = '=', abs(sol_final - hpd_dw), abs(sol_final - hpd_up)

            else:
                # if no mcmc data is available, we cannot check for degeneracy or limits,
                # so we just provide the solution and the hpd interval as error bars
                label, err_dw, err_up = '?', abs(sol_final - hpd_dw), abs(sol_final - hpd_up)

            setattr(self, 'l_'+par_name, label)
            [setattr(self, par_name+suffix, value) for suffix,value in
                zip(['','_eUP','_eDW'],[round(sol_final,5),round(err_up,5),round(err_dw,5)])]

        # In most cases l_logg should be 'd' but is not due to 70% is not enough
        if 'Teff' in self.parameters and 'lgf' in self.parameters:
            if getattr(self, 'l_Teff') == 'd' or getattr(self, 'l_lgf') == 'd':
                self.l_logg = 'd'

        # Lines used in the analysis
        line_weights = []; line_names = []; line_windows = []
        for i in ['balmer','helium1','helium2','silicon']:
            lines = soldata.obsdat.spectrum[0][i][0]
            line_weights += [j for j in lines.weight[0]]
            line_names += [j.decode() for j in lines.line[0]]
            line_windows += [[j,k] for j,k in zip(lines.lmin[0],lines.lmax[0])]

        self.line_weights = line_weights
        self.line_names = line_names
        self.line_windows = line_windows


def maui_results(input_list, output_dir, check_best=False, last_only=False, solution='max', FR=False,
    do_logg=False, do_pdf=False, pdflines='diag', grid_only=[], output_table=True, format_table='fits',
    black_theme=False):

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

    solution : str, optional
        Choose between 'max'/'smooth' to select either the maximum of the probability distribution,
        or the maximum of a smoothed distribution, for the solution provided in the output table.

    FR : boolean, optional
        True if the analyses correspond to Fast Rotator stars for which whe MAUI windows
        are different. Default is False.

    do_logg : boolean, optional
        If True, the logg is calculated from the Teff and lgf values replacing 'lgf' in the
        output plots. Default is False.

    do_pdf : boolean, optional
        If True, a pdf comparing the synthetic diagnostic lines with the original is made.

    pdflines : str/float, optional
        Choose between 'diag'/'all'/'def'/'<line>,<line>',float to select the lines to
        be used in the pdf plots. If <line> option is used, the line wavelengths must
        be a string with lines separated by commas or a single float.
        Default is 'diag'.

    grid_only : list, optional
        List of grid names to limit the output to those results analysed with such grid.

    output_table : boolean, optional
        If True, the table with the output data will be created. Default is True.

    format_table : str, optional
        Enter the output format for the table: 'fits' (default), 'ascii' or 'csv'.

    black_theme : boolean, optional
        If True, the plots will be generated with a black background. Default is False.

    Returns
    -------
    Nothing but the output table (+PDFs) with the MAUI results are generated.
    '''

    timenow = time.strftime('%Y%m%d_%H%M%S')

    if output_dir[-1] != '/':
        output_dir += '/'

    # Create the input list from a table.
    # NOTE: A column named 'ID' or 'filename' is required. If 'filename' if provided, then
    # the specific file among other possible solution files for the same source will be used.
    if input_list.endswith(('.txt','.fits')):
        stars = findtable(input_list)
        try:
            stars = [i.replace('.fits','') for i in stars['filename']]
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

    # Sort the list of stars by name
    stars = sorted(stars)

    # Set the progress bar
    bar = pb.ProgressBar(maxval=len(stars),
                         widgets=[pb.Bar('=','[',']'),' ',pb.Percentage(),' ',pb.ETA()])
    bar.start()

    # Create pdf file to save the plots of the results
    if do_pdf == True:

        from matplotlib.backends.backend_pdf import PdfPages
        if black_theme == True:
            plt.style.use('dark_background')

        if not os.path.exists(output_dir + 'plots/'):
            os.makedirs(output_dir + 'plots/')

        pdf_solution = PdfPages(output_dir + 'plots/MAUI_results_lines_%s.pdf' % timenow)
        pdf_makchain = PdfPages(output_dir + 'plots/MAUI_results_chain_%s.pdf' % timenow)

        plt.rcParams.update({
            'xtick.labelsize' : 6,
            'ytick.labelsize' : 6,
            'axes.titlesize' : 6,
            'axes.titlepad' : 3
            })

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
            msg.warn('No solution*.idl file found for %s. Continuing...' % name)
            continue

        if len(matches) > 1 and last_only == True:
            matches = [sorted(matches, key=lambda x: int(x[-14:-4].replace('-','')), reverse=True)[0]]

        for match in matches:
            # Search for the corresponding markov chain file
            mcmcfile = match.split('emulated_solution_')[1]
            # if the MARKOV_CHAIN folder does not exist it will continue but without the mcmc information
            if not 'MARKOV_CHAIN' in os.listdir(output_dir):
                msg.warn('MARKOV_CHAIN/ folder not found in output directory. Continuing without mcmc information...')
                mcmcfile = None
            # if the MARKOV_CHAIN folder exists but the specific mcmc file does not exist, the file is skipped
            if mcmcfile is not None and mcmcfile not in os.listdir(output_dir + 'MARKOV_CHAIN/'):
                msg.warn('Associated mcmc file for %s not found in MARKOV_CHAIN/ folder...' % match.split('SOLUTION/')[1])
                continue
            elif mcmcfile is not None:
                mcmcfile = output_dir + 'MARKOV_CHAIN/' + mcmcfile

            # Load the idl class for the file
            star = solution_maui(match, mcmcfile=mcmcfile, solution=solution)

            # Skip the used grid if not selected from the grid_only keyword
            if grid_only != [] and not star.gridname in grid_only:
                continue

            # Find the best SNR spectra in the DB
            if check_best == True:
                best_SNR = spec(star.id_star, snr='bestHF')

            # Check if the input file matches with the best SNR spectra available
            if check_best == True and star.filename != best_SNR.filename:
                msg.warn('%s does not match with best spectrum available.' % star.filename)

            # Generate the table with the results of each of the parameters (basic and with errors)
            data_row = [star.id_star, star.filename, star.gridname, star.BC, star.B_V0, star.vsini, star.vmac]

            [data_row.extend([
                getattr(star, 'l_'+par_name),
                getattr(star, par_name),
                getattr(star, par_name+'_eUP'),
                getattr(star, par_name+'_eDW'),
                ]) for par_name in dic_maui_param]

            # Append the raw to the table
            data_rows.append(tuple(data_row))

            if do_pdf == True:

                # PLOT OF SPECTRAL LINES OUT OF THE IDL SOLUTION FILES
                if pdflines == 'diag':
                    mask = [i > 0 for i in star.line_weights]
                    line_names = np.asarray(star.line_names)[mask].tolist()
                    lines_lwl,lines_rwl = np.asarray(star.line_windows).T
                    lines_lwl = np.asarray(lines_lwl)[mask].tolist()
                    lines_rwl = np.asarray(lines_rwl)[mask].tolist()
                    line_colors = ['g']*len(line_names)

                elif pdflines == 'all' or pdflines == '*':
                    line_names = star.line_names
                    lines_lwl,lines_rwl = np.asarray(star.line_windows).T
                    #lines_lwl = np.asarray(lines_lwl).tolist()
                    #lines_rwl = np.asarray(lines_rwl).tolist()
                    line_colors = ['g' if i > 0 else 'r' for i in star.line_weights]

                elif pdflines == 'def' or pdflines == 'default':
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

                else:
                    if type(pdflines) == float:
                        line_names = [str(pdflines)]
                    elif ',' in pdflines:
                        line_names = pdflines.split(',')
                    else:
                        line_names = [pdflines]

                    lines_lamb = [float(i) for i in line_names]
                    lines_lwl = [i-10 for i in lines_lamb]
                    lines_rwl = [i+10 for i in lines_lamb]
                    line_colors = ['g']*len(line_names)

                nrows, ncols = even_plot(len(line_names))

                fig, ax = plt.subplots(nrows, ncols, figsize=(13,8))

                fig_title = ''
                for par_name in dic_maui_param:
                    val = getattr(star, par_name)
                    if np.isnan(val) == True:
                        continue
                    label = getattr(star, 'l_'+par_name)
                    if label == 'd': label = '=d.'
                    if par_name == 'Teff':
                        val = val*1e4

                    fig_title += par_name + label + str(round(val,2)) + '  '

                fig.suptitle(star.id_star + ' -- ' + match.split('emulated_solution_mcmc_sqexp_mat1_')[-1]
                    + ' -- ' + star.gridname + '\n' + fig_title, fontsize=8)

                if nrows == ncols == 1:
                    axs = [ax]
                else:
                    axs = ax.flatten()

                for j,line_lwl,line_rwl,line_name,c in zip(range(len(line_names)),lines_lwl,lines_rwl,line_names,line_colors):
                    #mask = (star.synwave > line_lwl) & (star.synwave < line_rwl)
                    #ax_i.plot(star.synwave[mask], star.synflux[mask], color='gray', lw=.3)
                    mask = (star.obswave > line_lwl) & (star.obswave < line_rwl)
                    window_wave = star.obswave[mask]
                    window_flux = star.obsflux[mask]
                    axs[j].plot(window_wave, window_flux, color='k', lw=.7)

                    weight_mask = np.zeros(window_wave.shape, dtype=bool)
                    if FR == False:
                        for lmin, lmax in mask_maui_SR:
                            weight_mask |= (window_wave >= lmin) & (window_wave <= lmax)
                    elif FR == True:
                        for lmin, lmax in mask_maui_FR:
                            weight_mask |= (window_wave >= lmin) & (window_wave <= lmax)

                    # Where the condition is True, weight is 0. Otherwise, it's 1.
                    weight = np.where(weight_mask, 0, 1).tolist()

                    # This is a visual trick to place the synthetic spectra where it really is ifthe
                    # normalization option is used, as this is not stored in the solution*.idl file
                    scale = np.sum(window_flux*star.synconv[mask]*weight)\
                                /np.sum(star.synconv[mask]*star.synconv[mask]*weight)
                    star.synconv_scaled = scale * star.synconv[mask]

                    axs[j].plot(window_wave, star.synconv_scaled, color=c, ls='--', lw=1)

                    if axs[j].get_ylim()[0] > 0.875:
                        axs[j].set_ylim(bottom=0.875)
                    if axs[j].get_ylim()[1] < 1.025:
                        axs[j].set_ylim(top=1.025)

                    # Plot the region with weight = 1
                    ymean = np.asarray(axs[j].get_ylim()).mean()
                    weight = [None if i==0 else ymean for i in weight]
                    axs[j].plot(window_wave, weight, c='dodgerblue', lw=1, alpha=0.5)

                    axs[j].set_title(line_name)
                    axs[j].tick_params(direction='in', top='on', right='on')
                    axs[j].minorticks_on()

                [fig.delaxes(axs[i]) for i in np.arange(len(line_names), len(axs), 1)]

                fig.tight_layout()
                fig.subplots_adjust(wspace=0.15, hspace=0.3)

                pdf_solution.savefig(fig); plt.close(fig)

                # Plot of probabillity distributions out of the idl markov-chain files
                mcmcdata = readsav(mcmcfile)

                parameters = [var.decode() for var in mcmcdata.varname]

                nrows, ncols = even_plot(len(parameters))

                fig, ax = plt.subplots(nrows, ncols, figsize=(13,8))
                fig.subplots_adjust(wspace=.5, hspace=.5)

                fig.suptitle(star.id_star + ' -- ' + match.split('emulated_solution_mcmc_sqexp_mat1_')[-1]
                    + ' -- ' + star.gridname + '\n' + fig_title, fontsize=8)

                axs = ax.flatten()

                for j in range(len(parameters)):

                    chain = mcmcdata.chain_final.T[j]*(mcmcdata.xmax[j] - mcmcdata.xmin[j]) + mcmcdata.xmin[j]

                    if parameters[j] == 'logQs':
                        chain -= 10

                    # Replace lgf by logg
                    if do_logg == True and 'lgf' in parameters and parameters.index('lgf') == j:
                        idx = parameters.index('Teff')
                        teff = mcmcdata.chain_final.T[idx]*(mcmcdata.xmax[idx] - mcmcdata.xmin[idx]) + mcmcdata.xmin[idx]
                        chain = chain + 4*np.log10(teff)
                        parameters[j] = 'logg'

                    iqr = np.quantile(chain, q=[.25, .75])
                    fd_bin = 2*np.diff(iqr)/(len(chain)**(0.3))

                    # Change the original fd_bin value to the nearest one to fit within the range of the data evenly
                    fd_bin = (max(chain) - min(chain))/np.round((max(chain) - min(chain))/fd_bin[0])

                    weights = np.ones_like(chain)/float(len(chain))
                    axs[j].hist(chain, bins=np.arange(min(chain), max(chain) + fd_bin, fd_bin),
                        weights=weights, histtype='stepfilled', fc='gray', ec='g', lw=1, alpha=0.6)

                    par_val = getattr(star, parameters[j])

                    # plot the values of sol_max and the hdp intervals
                    axs[j].axvline(par_val, ls='--', c='r', label='sol_%s + HDP' % solution)

                    if not 'd' in getattr(star, 'l_'+parameters[j]):
                        axs[j].axvline(par_val-getattr(star, parameters[j]+'_eDW'), ls=':', lw=2,c='r')                        
                        axs[j].axvline(par_val+getattr(star, parameters[j]+'_eUP'), ls=':', lw=2, c='r')

                    # plot the median and the IQR intervals
                    axs[j].axvline(np.median(chain), ls='--', c='orange', label='median + IQR')
                    axs[j].axvline(iqr[0], ls=':', c='orange')
                    axs[j].axvline(iqr[1], ls=':', c='orange')

                    # add the legend
                    if j == 0:
                        fig.legend(fontsize=8, loc='upper left', handlelength=3)

                    # add parameter in bold  and the sol_ and median to the title
                    axs[j].set_title(r"$\bf{%s}$ -- sol_%s = %.5f, median = %.5f" % (parameters[j],\
                        solution, par_val, np.median(chain)), fontsize=8)
                    axs[j].tick_params(direction='in', top='on', right='on')
                    axs[j].minorticks_on()

                [fig.delaxes(axs[i]) for i in np.arange(len(parameters), len(axs), 1)]

                fig.tight_layout()

                pdf_makchain.savefig(fig); plt.close(fig)

        bar.update(i)

    if do_pdf == True:
        pdf_solution.close()
        pdf_makchain.close()

    if black_theme == True:
        plt.style.use('default')

    bar.finish()

    # Saving the results:
    # - Add names for the basic columns
    names = ['ID','filename','Grid_name','BC','BV_0','vsini','vmac']

    # - Add names for the parameter column names
    for i in range(len(dic_maui_param)):
        names += ['l_'+ list(dic_maui_param.keys())[i],list(dic_maui_param.values())[i],\
                        list(dic_maui_param.values())[i]+'_eUP',list(dic_maui_param.values())[i]+'_eDW']

    output = Table(rows=data_rows, names=(names))

    # remove all the columns with only 'nan' values
    for col in output.colnames:
        # if the column is of a numeric type and all values are nan, or is a string column and all values are empty, remove it
        if ((output[col].dtype.kind in 'if' and np.isnan(output[col]).all()) or\
            (output[col].dtype.kind in 'U' and (output[col] == '').all())):
            output.remove_column(col)

    full_path = (output_dir + 'MAUI_results_%s.' + format_table) % timenow

    if format_table == 'ascii':
        format_table += '.fixed_width_two_line'
        full_path = full_path.replace('.ascii', '.txt')

    if output_table == True:
        output.write(full_path, format=format_table, overwrite=True)

    # print the results in the terminal if only one star is in the output table
    if len (output) == 1:
        for row in output:
            msg.info('\nFile & grid: %s -- %s' % (row['filename'], row['Grid_name']))
            print('Param  :  l  sol      err_d    err_u')
            for par_name in dic_maui_param:
                if par_name in output.colnames:
                    label = row['l_'+par_name]
                    if label == 'd': label = '=d.'
                    val = row[par_name]
                    err_dw = row[par_name+'_eDW']
                    err_up = row[par_name+'_eUP']
                    print('%-6s :  %s  %.5f  %.5f  %.5f' % (par_name, label, val, err_dw, err_up))

    msg.bold('g','Finished!')

    return None


def compare_results(table_1, table_2, path_t1=None, path_t2=None, par_name='*', exclude=None,
                    save_plot=False):
    '''
    Function to compare the results of two tables containing the results of MAUI analyses.

    Parameters
    ----------
    table_1 : str
        First table to compare. E.g. 'table_1.txt/fits'.

    table_2 : str
        Second table to compare. E.g. 'table_2.txt/fits'.

    path_t1 : str, optional
        Path to the first table. Default is None.

    path_t2 : str, optional
        Path to the second table. Default is None.

    par_name : str, optional
        Name of the parameter to compare.
        Default is '*' to compare all parameters in common between the two tables.

    exclude : list, optional
        List of stars to exclude from the comparison. Default is None.

    save_plot : bool, optional
        Whether to save the plot. Default is False.

    Returns
    -------
    Nothing but a plot comparing the results of the two tables is generated.
    '''

    t1 = findtable(table_1, path=path_t1)
    t2 = findtable(table_2, path=path_t2)

    if t1 is None or t2 is None:
        print('ERROR: Problem loading the tables. Exiting...')
        return None

    t = join(t1, t2, keys='ID', table_names=['t1','t2'], join_type='inner')
    if len(t) == 0:
        print('ERROR: No common IDs between the two tables. Exiting...')
        return None

    if exclude is not None and type(exclude) == str:
        exclude = exclude.split(',')
        t = t[~np.isin(t['ID'], exclude)]

    if par_name == '*':
        par_name = [i for i in t1.colnames if i in t2.colnames and 
                    i not in ['ID','filename','Grid_name'] and 
                    (not np.ma.is_masked(t1[i]) and not np.ma.is_masked(t2[i])) and 
                    not i.endswith(('_eUP','_eDW')) and not i.startswith('l_')]
        par_name = [i for i in par_name if i not in ['vsini','vmac']]
    else:
        par_name = par_name.split(',') if ',' in par_name else [par_name]
        par_name = [i for i in par_name if i in t1.colnames and i in t2.colnames and 
                    (not np.ma.is_masked(t1[i]) and not np.ma.is_masked(t2[i])) and 
                    not i.endswith(('_eUP','_eDW')) and not i.startswith('l_')]

    n_pars = [i for i in par_name if i in t1.colnames and i in t2.colnames]
    n_rows, n_cols = even_plot(len(n_pars))

    plt.style.use('dark_background')

    fig1, ax1 = plt.subplots(n_rows, n_cols, figsize=(12,7.5))
    fig1.subplots_adjust(wspace=0.3, hspace=0.3)

    for i, par in enumerate(n_pars):
        # print the name of the star if the parameter values are outside the MAUI uncertainties
        if par in dic_maui_uncertainties:
            outliers = t[abs(t[par+'_t2']-t[par+'_t1']) > dic_maui_uncertainties[par]]['ID'].tolist()
            if len(outliers) > 0:
                print(f'Stars with difference in {par} outside uncertainties: {outliers}')

        ax_i = ax1.flatten()[i]
        x = t[par+'_t1']
        y = t[par+'_t2']
        ax_i.scatter(x, y-x, color='w', s=20, alpha=0.7)
        ax_i.plot([min(x), max(x)], [0,0], '--', color='b')
        if par in dic_maui_uncertainties:
            ax_i.axhline(dic_maui_uncertainties[par], ls=':', color='r')
            ax_i.axhline(-dic_maui_uncertainties[par], ls=':', color='r')
        
        ax_i.set_xlabel(par + ' (t1)')
        ax_i.set_ylabel(par + ' (t2) - ' + par + ' (t1)')
        ax_i.set_title(par)
        ax_i.tick_params(direction='in', top='on', right='on')
        ax_i.minorticks_on()
    
    fig1.suptitle('Comparison between %s (t1) and %s (t2)' % (table_1, table_2), fontsize=10)
    [fig1.delaxes(ax1.flatten()[i]) for i in np.arange(len(n_pars), len(ax1.flatten()), 1)]
    fig1.tight_layout()
    plt.show(block=False)

    fig2, ax2 = plt.subplots(1, 3, figsize=(12,3.5))
    fig2.subplots_adjust(wspace=0.3, hspace=0.3)
    if 'Teff' in t1.colnames and 'logg' in t1.colnames and 'Teff' in t2.colnames and 'logg' in t2.colnames:
        ax2[0].scatter(t['Teff_t1']-t['Teff_t2'], t['logg_t1']-t['logg_t2'], color='w', s=20)
        ax2[0].axhline(0, ls='--', color='b'); ax2[0].axvline(0, ls='--', color='b')
        ax2[0].plot([-dic_maui_uncertainties['Teff'],-dic_maui_uncertainties['Teff'],+dic_maui_uncertainties['Teff'],
                    +dic_maui_uncertainties['Teff'],-dic_maui_uncertainties['Teff']],\
                    [-dic_maui_uncertainties['logg'],+dic_maui_uncertainties['logg'],+dic_maui_uncertainties['logg'],
                    -dic_maui_uncertainties['logg'],-dic_maui_uncertainties['logg']],\
                    ls=':', color='r')
        ax2[0].set_xlabel('Teff (t1) - Teff (t2)')
        ax2[0].set_ylabel('logg (t1) - logg (t2)')
        ax2[0].tick_params(direction='in', top='on', right='on')
    else:
        ax2[0].set_visible(False)
    
    if 'He' in t1.colnames and 'Micro' in t1.colnames and 'He' in t2.colnames and 'Micro' in t2.colnames:
        ax2[1].scatter(t['He_t1']-t['He_t2'], t['Micro_t1']-t['Micro_t2'], color='w', s=20)
        ax2[1].axhline(0, ls='--', color='b'); ax2[1].axvline(0, ls='--', color='b')
        ax2[1].plot([-dic_maui_uncertainties['He'],-dic_maui_uncertainties['He'],+dic_maui_uncertainties['He'],
                    +dic_maui_uncertainties['He'],-dic_maui_uncertainties['He']],\
                    [-dic_maui_uncertainties['Micro'],+dic_maui_uncertainties['Micro'],+dic_maui_uncertainties['Micro'],
                    -dic_maui_uncertainties['Micro'],-dic_maui_uncertainties['Micro']],\
                    ls=':', color='r')
        ax2[1].set_xlabel('He (t1) - He (t2)')
        ax2[1].set_ylabel('Micro (t1) - Micro (t2)')
        ax2[1].tick_params(direction='in', top='on', right='on')
    else:
        ax2[1].set_visible(False)
    
    if 'Si' in t1.colnames and 'Micro' in t1.colnames and 'Si' in t2.colnames and 'Micro' in t2.colnames:
        ax2[2].scatter(t['Si_t1']-t['Si_t2'], t['Micro_t1']-t['Micro_t2'], color='w', s=20)
        ax2[2].axhline(0, ls='--', color='b'); ax2[2].axvline(0, ls='--', color='b')
        ax2[2].plot([-dic_maui_uncertainties['Si'],-dic_maui_uncertainties['Si'],+dic_maui_uncertainties['Si'],
                    +dic_maui_uncertainties['Si'],-dic_maui_uncertainties['Si']],\
                    [-dic_maui_uncertainties['Micro'],+dic_maui_uncertainties['Micro'],+dic_maui_uncertainties['Micro'],
                    -dic_maui_uncertainties['Micro'],-dic_maui_uncertainties['Micro']],\
                    ls=':', color='r')
        ax2[2].set_xlabel('Si (t1) - Si (t2)')
        ax2[2].set_ylabel('Micro (t1) - Micro (t2)')
        ax2[2].tick_params(direction='in', top='on', right='on')
    else:
        ax2[2].set_visible(False)

    fig2.tight_layout()
    plt.show(block=False)

    if save_plot == True:
        tail_name = (table_1.split('_')[-1].replace('.fits', ''), table_2.split('_')[-1].replace('.fits', ''))
        fig1.savefig(maindir + 'plots/MAUI/comparison_%s_vs_%s.pdf' % tail_name, format='pdf')
        fig2.savefig(maindir + 'plots/MAUI/comparison_deltas_%s_vs_%s.pdf' % tail_name, format='pdf')

    plt.style.use('default')

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

    # Change ...dic_maui_grids][1:] for individual grids
    for name,_,_,box in [dic_maui_grids[i] for i in dic_maui_grids][1:]:
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
        print(color.error+'missing column names. Exiting...\n' + color.end)
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


def gen_synthetic(output_dir, convolve=True, lwl=3900, rwl=5080):

    '''
    Function to generate the .dat files of the synthetic spectra generated by
    MAUI in all the .idl files found in the SOLUTION folder.

    Parameters
    ----------
    output_dir : str
        Enter the path to the OUTPUT folder where the SOLUTION sub-folder containing the
        *solution*.idl files are. # Note: this could be your "mauidir" variable itself.

    colvolve : boolean, optional
        If True, the synthetic spectra are convolved with the rotational and macroturbulent
        broadening profiles, plus the instrumental profile. Default is True.

    lwl : float, optional
        Sets the start wavelenght of the output spectrum.

    rwl : float, optional
        Sets the end wavelenght of the output spectrum.

    Returns
    -------
    Nothing but the ascii .dat files are generated.
    '''

    save_dir = datadir + 'ASCII/Synthetic_MAUI/'

    if not output_dir.endswith('/'):
        output_dir += '/'
    if output_dir.endswith('SOLUTION/'):
        output_dir = output_dir[:-9]

    for file in os.listdir(output_dir + 'SOLUTION/'):
        if not file.startswith('._') and file.endswith('.idl'):
            star = solution_maui(output_dir + 'SOLUTION/' + file,
                                 output_dir + 'MARKOV_CHAIN/' + file.replace('emulated_solution_',''))

            #star.filename = star.filename.replace(str(star.resolution),'85000')
            new_star = '%s_red%i.dat' % (star.filename[:-5],dic_maui_grids[star.gridname][2])

            # Check if the file already exists and ask if it should be overwritten
            if os.path.exists(save_dir + new_star):
                print('File %s already exists. Overwrite? (y/n)' % new_star)
                answer = input()
                if answer == 'n':
                    continue

            # Save the non-convolved synthetic spectrum in the ASCII folder
            np.savetxt(save_dir + new_star, np.c_[star.synwave,star.synflux],
                       fmt=('%.4f','%.6f'), header='lambda    flux', comments='')

            star_idl = spec(new_star, orig='txt', delimiter=' ')
            star_idl.waveflux(lwl, rwl, delimiter=' ')
            #plt.plot(star_idl.wave, star_idl.flux, 'r', lw=.3) # plot to check
            if convolve == True:
                star_idl.degrade(resol=star_idl.resolution, profile='rotmac',
                                 vsini=star.vsini, vmac=star.vmac)
                #plt.plot(star_idl.wave, star_idl.flux, 'g', lw=.3) # plot to check

            # Save the convolved synthetic spectrum in the ASCII folder
            np.savetxt(save_dir + new_star, np.c_[star_idl.wave,star_idl.flux],
                       fmt=('%.4f','%.6f'), header='lambda    flux', comments='')


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
