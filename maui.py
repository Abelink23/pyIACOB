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
'nlte_10.4.7_bsgs_SOLAR_expoclump_n12345o123c234mg2si234djl_v1_2021-05-05.idl': ('BSg_CNOSiMg','DeepPink',6,
[[4.148,4.477,4.477,4.148,4.148],[3.392,3.392,4.386,4.386,3.392]])
}

def mauipath(path=None):

    '''
    Function to set the main directory of MAUI.

    Parameters
    ----------
    path : str, optional
        If 'default' or 'def', it will choose the default path (see below).

    Returns
    -------
    Selected MAUI main directory path.

    Notes
    -----
    As default, it is configured for one user specific paths.
    This must be modified according to your particular need or case.
    '''

    if path in ['def','default']:

        if platform.system() == 'Darwin':
            defmainpath = '/Users/abelink/Documents/MAUI/'

        elif 'iac.es' in platform.uname().node:
            defmainpath = '/net/nas/proyectos/hots/masblue/maui2021/'

        mainpath = defmainpath

    elif path == None:

        mainpath = input('Working directory path (default is %s) : ' % defmainpath)

        if mainpath == '':
            mainpath = defmainpath

    else: mainpath = path

    return mainpath

mauidir = mauipath('def')


def maui_input(table='IACOB_O9BAs_SNR20.fits', output_name='MAUI_input',
    RV0tol=200,ascii_0=False):

    '''
    Function to generate the input table for MAUI given an input table with the
    target stars and quality flags for the line fittings.
    Optionally, the function allows to generate the input ascii spectra for MAUI
    subtracting the individual radial velocity.

    IMPORTANT NOTE 1 - Make sure you update some other input tables contaning:
        1) RV, EWs, FWs --> table_REF('')
        2) IACOB-broad parameters --> table_IB('')
        3) Results from MAUI --> results('')
        (They are listed at the beginning of the funcion.)

    IMPORTANT NOTE 2 - The stars in the input table MUST also be in the other
    tables with the same name except table 3) Results from MAUI

    Parameters
    ----------
    table : str, optional
        Enter the input table contaning a column 'ID' with the name of the stars.

    output_name : str, optional
        Enter the name for the output table. Default is 'MAUI_verX'.

    RV0tol : int, optional
        Enter the input radial velocity tolerance for the radial velocity correction.

    ascii_0 : boolean, optional
        If True, ascii files will be created for each of the input sources.
        Default is False.

    Returns
    -------
    Nothing but the MAUI input file is generated.
    '''

    if type(table) is type(Table()): pass # In case the input table is already a table
    else: table = findtable(table) # file where star names and quality flags are
    table_REF = findtable('RVEWFWs_O9BAs.fits') # file where RVs, EWs and FWs are
    table_IB = findtable('IB_results_ver5.txt') # file where vsini and vmac are
    results = findtable('MAUI_results.fits') # file with output from MAUI

    maui_txt = open(maindir + 'lists/%s.txt' % output_name,'a')
    maui_txt.write(
        "{:<40}".format('fullname')+"{:<48}".format('filename')+"{:<6}".format('vrad')+\
        "{:<6}".format('vsini')+"{:<7}".format('evsini')+"{:<5}".format('vmac')+\
        "{:<7}".format('evmac')+"{:<7}".format('R')+"{:<4}".format('SNR')+\
        ' ;# SpC       FW34-14 SiIII SiII l lTef l lgf  Grid\n')

    quit = ''
    for row in table:

        ascii = ascii_0

        if quit == 'quit': break

        id = row['ID'].strip()

        # Filer based on properties from the main table:
        if 'SB' in table.columns and 'SB2' in row['SB']: continue
        if 'CHb' in table.columns:
            if 'Em' in row['CHb']:
                if not 'Em(p)' in row['CHb']: continue
            if 'PCyg' in row['CHb']: continue
        if 'QIB' in table.columns and row['QIB'] < 2: continue

        match_REF = table_REF[[i.strip()==id for i in table_REF['ID']]]
        match_IB = table_IB[table_IB['ID']==id]

        if len(match_REF) == 0 or len(match_IB) == 0:
            print('Missing information in RVEWFW or IB tables for %s\n' % id); continue

        # Filter based on Si lines properties:
        if row['QSiII']<3 or match_REF['EWSiII']<50 or np.isnan(match_REF['EWSiII']) or\
            match_REF['depSiII']<3/match_REF['SNR_B']: SiIIFG = 0
        else: SiIIFG = 1

        if row['QSiIII']<3 or match_REF['EWSiIII1']<50 or np.isnan(match_REF['EWSiIII1']) or\
            match_REF['depSiIII1']<3/match_REF['SNR_B']: SiIIIFG = 0
        else: SiIIIFG = 1

        if SiIIIFG == 0 and SiIIFG == 0:
            print('No SiIII/SiII found for %s\n' % id); continue

        star = spec(id, SNR='best')

        if match_IB['filename'][0] != star.filename:
            print('Warning: Different files from best SNR and from IB results for %s' % name)
            print(star.filename,' vs ',match_IB['filename'][0])
            if ascii_0 == True:
                do_file = input('Which ascii do you want to create 1 or 2: ')
                if int(do_file) == 1: star.filename = star.filename
                elif int(do_file) == 2: star.filename = match_IB['filename'][0]

        if ascii_0 == True and not search(star.filename[:-5]+'_RV.ascii',\
        os.path.expanduser('~')+'/Documents/MAUI/ASCII/') is None: ascii = False

        # ----------------------------------------------------------------------
        # Extra information appended to the end of each row:
        match_results = results[results['ID']==id]
        if len(match_results) == 0:
            l_Tef = logTf = l_lgf = loggf = grid = 0
        else:
            l_Tef = match_results['l_Teff'][0]; logTf = 4+np.log10(match_results['Teff'][0])
            l_lgf = match_results['l_lgf'][0]; loggf = match_results['lgf'][0];
            grid = grids_dic[match_results['Model_name'][0]][0]
        # ----------------------------------------------------------------------

        if ascii == True:

            from RV import RV0

            # If RV SiIII is good enough, it uses it for the rv0 correction:
            if abs(match_REF['RVSiIII1']-match_REF['RVHb']) < 10:
                star.rv0 = match_REF['RVSiIII1']
                star.waveflux() # Applies the rv0 correction to the wavelenght vector.
                print('RV good enough.\n')

            # Otherwise it calculates the rv0 correction with the RV package:
            else:
                next = 'n'
                while next == 'n':

                    skip = input('%s - Hit return to continue, type "s" to skip: ' % id)
                    if skip == 's': break

                    if table['SpT_code'] <= 2.5: star.plotspec(4530,4590)
                    else: star.plotspec(6337.11,6357.11)

                    SpT = '-'
                    while SpT not in ['O','B','A']:
                        SpT = input('Choose SpT for the RV0 list of lines (default is B): ')
                        if SpT == '' or SpT == 'B': SpT = 'B'; spt_list = 'rv_Bs.lst'
                        elif SpT == 'A': SpT = 'A'; spt_list = 'rv_As.lst'
                        elif SpT == 'O': SpT = 'O'; spt_list = 'rv_Os.lst'
                    fun = '-'
                    while fun not in ['g','l','v','r','vr']:
                        fun = input('Choose function to fit between g/l/v/r/vr (default is g): ')
                        if fun == '': fun = 'g'
                    wid = '-'
                    while type(wid) is not float:
                        wid = input('Choose the initial width in angstroms (default is 15): ')
                        if wid == '': wid = 15.
                        else: wid = float(wid)

                    plt.close()

                    star.rv0 = RV0(spt_list,star.spectrum,ewcut=30,width=wid,tol=RV0tol,func=fun)
                    star.waveflux() # Applies the rv0 correction
                    #star.cosmic(sigclip=0.005)

                    if table['SpT_code'] <= 2.5: star.plotspec(4530,4590,poslines='OB')
                    else: star.plotspec(6337.11,6357.11,poslines='OB')

                    input(); plt.close('all')

                    next = input('Type "n" to repeat, hit return to move to the next star. ')

            star.export(tail='_RV',extension='.ascii')

        # MAUI input last modifications:
        if star.resolution == 67000: star.resolution = 85000

        if match_IB['vsini'] < 5:
            match_IB['vsini'] = 0
            match_IB['evsini'] = 0

        if match_IB['vmac'][0] < 10:
            match_IB['vmac'] = 0
            match_IB['evmac'] = 0

        if match_REF['SNR_B'][0] > 200:
            match_REF['SNR_B'] = 200

        maui_txt.write(
            "{:<40}".format(star.filename[:-5])+\
            "{:<48}".format(star.filename[:-5]+'_RV.ascii')+"{:<6}".format('0.0d0')+\
            "{:<6}".format(str(int(round(match_IB['vsini'][0],0))))+\
            "{:<7}".format(str(int(round(match_IB['evsini'][0]))))+\
            "{:<5}".format(str(int(round(match_IB['vmac'][0]))))+\
            "{:<7}".format(str(int(round(match_IB['evmac'][0]))))+\
            "{:<7}".format(str(star.resolution)+'.')+\
            "{:<5}".format(str(int(round(match_REF['SNR_B'][0],0))))+';# '+\
            "{:<12}".format(row['SpC'].strip().replace(' ',''))+' '+\
            "{:<8}".format(str(round(match_REF['FW34Hb'][0]-match_REF['FW14Hb'][0],1)))+\
            "{:<5}".format(str(SiIIIFG))+"{:<2}".format(str(SiIIFG))+' '+\
            str(l_Tef)+' '+"{:<4}".format(str(round(logTf,2)))+' '+\
            str(l_lgf)+' '+"{:<5}".format(str(round(loggf,2)))+\
            "{:<15}".format(str(grid))+\
            '\n')

    maui_txt.close()

    if ascii == True:
        print('Remember to move the new ascii into "ASCII_ARCHIVE" folder')

    return 'DONE'


def maui_results(input_table='MAUI_ver11.txt', solution_dir='server', check_best=True,
    pdfplots=False, grids_table='MAUI_grid_limits.fits', format='fits'):

    '''
    Function to generate a table with the results from MAUI given an input table
    containing the ID of the stars and filename to search in the MAUI-SOLUTION
    directory.

    Parameters
    ----------
    input_table : str, optional
        Name of the input table contaning the list of stars to search.

    solution_dir : str, optional
        Enter the directory where to locate the input_table and output fits.
        Shortcuts are 'local' and 'server'. Default is 'server'.

    check_best : boolean, optional
        True if each spectra from the input_table is checked against the best
        spectrum in the database. Default is True.

    pdfplots : boolean, optional
        If True, a pdf comparing the synthetic diagnostic lines with the original is made.

    grids_table : str, obtional
        Name of the table containing the limits of the grids in MAUI.

    format : str, optional
        Enter the output format for the table: 'fits' (default), 'ascii' or 'csv'.

    Returns
    -------
    Nothing but the output table (+PDFs) with the MAUI results are generated.
    '''

    # Set the paths to the solutions directory
    if solution_dir == 'local':
        solution_dir = mauidir + 'SOLUTION/'
    elif solution_dir == 'server':
        solution_dir = mauidir + 'RESULTS_BSGS_202101/SOLUTION/'

    # Load the tables
    table = findtable(input_table)
    grids = findtable(grids_table)

    # Create pdf file to save the plots of the results
    if pdfplots == True:

        from RV import RV0

        from matplotlib.backends.backend_pdf import PdfPages
        pp = PdfPages(maindir + 'plots/MAUI/%s.pdf' % input_table)

    param_lst = ['Teff','lgf','He','Micro','logQs','beta','C','N','O','Mg','Si',
                 'S','Fe','Ti','fcl','vcl']

    bar = pb.ProgressBar(maxval=len(table),
                         widgets=[pb.Bar('=','[',']'),' ',pb.Percentage()])
    bar.start()

    data_rows = []
    for row,i in zip(table,range(len(table))):
        filename = row['filename']

        match = []
        for file in os.listdir(solution_dir):
            if file.endswith('.idl') and \
            file.split('_sqexp_mat1_')[1][:-14] + 'RV.ascii' == filename:
                match.append(solution_dir + file)

        if len(match) == 0:
            print('\nWARNING: No .idl file found for %s. Continuing...' % filename)
            continue

        elif len(match) > 1:
            match = sorted(match, key=lambda x: int(x[-14:-4].replace('-','')), reverse=True)

        # Load the .idl file
        idldata = readsav(match[0])

        # Find which grid has been used for the analysis
        filename = idldata.aa[0][3].decode()
        id = filename.split('_')[0]
        grid  = grids[grids['Model_name'] == idldata.modelgridname.decode()]

        # Find the best SNR spectra in the DB
        if check_best == True or pdfplots == True:
            best_SNR = spec(id,SNR='best')

        # Check if the input file matches with the best SNR spectra available
        if check_best == True and filename != best_SNR.filename[:-5]:
            print('\nWARNING: %s does not match with best spectrum available.'
            % filename[:-5])

        # Generate the table with the results of each of the parameters
        data_row = []; data_row.extend([id])
        parameters = [j.decode() for j in idldata.solution.var_label[0]] # 11/13 parameters
        for par_name in param_lst:

            if not par_name in parameters:
                data_row.extend(['',np.nan,np.nan,np.nan]); continue

            idx = parameters.index(par_name)

            #median = idldata.solution[0].sol[idx] # not used
            err_sim = idldata.solution[0].err_sol[0][idx] # symetric error of distribution at 67%
            if err_sim > abs(grid[par_name+'_UP']-grid[par_name+'_DW'])*0.2:
                sol_max = np.nan
                err_dw = grid[par_name+'_DW'][0] # Lower limit of the grid for the param
                err_up = grid[par_name+'_UP'][0] # Upper limit of the grid for the param
                label = 'd'

            else:
                sol_max = idldata.solution[0].sol_max[0][idx] # maximum of distribution
                err_dw = abs(sol_max-idldata.solution[0].hpd_interval[0][idx])
                err_up = abs(sol_max-idldata.solution[0].hpd_interval[1][idx])

                if par_name == 'logQs': sol_max -= 10

                if err_up < err_sim/2 and err_dw > err_sim/2: label = '>'
                elif err_up > err_sim/2 and err_dw < err_sim/2: label = '<'
                else: label = '='

            data_row.extend([label,round(sol_max,5),round(err_up,5),round(err_dw,5)])

        # Add the name of the model used in the MAUI emulator
        data_row.extend([grid['Model_name'][0]])

        # Add photometric data
        data_row.extend([idldata.phot_prim[0],idldata.phot_prim[2]])

        # Append the raw to the table
        data_rows.append(tuple(data_row))

        if pdfplots == True:

            lines_name = [
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

            rv0 = RV0(4552.622, findstar(filename+'.fits'), ewcut=50, width=20, tol=150, func='r')
            best_SNR.rv0 = rv0
            best_SNR.waveflux()
            best_SNR.spc()

            fig, ax = plt.subplots(4, 4, figsize=(10,10), tight_layout=True)
            fig_title = ''
            for par,val in zip(parameters,idldata.solution[0].sol_max[0]):
                if par == 'Teff':
                    val = val*1e4
                    Teff = val
                if par == 'lgf':
                    par = 'logg'
                    val = val + 4*np.log10(Teff) - 16
                fig_title += par+'='+str(round(val,2))+'  '

            fig.suptitle(id + ' -- ' + best_SNR.SpC + ' -- ' + match[0].split('/')[-1]
                + ' -- ' + grid['Model_name'][0] + '\n' + fig_title, fontsize=8)

            axs = ax.flatten()
            for ax_i,line_lamb,line_name in zip(axs,lines_lamb,lines_name):
                mask = [(best_SNR.wave > line_lamb-10) & (best_SNR.wave < line_lamb+10)]
                ax_i.plot(best_SNR.wave[mask], best_SNR.flux[mask], color='gray')
                mask = [(idldata.xobs_out > line_lamb-10) & (idldata.xobs_out < line_lamb+10)]
                ax_i.plot(idldata.xobs_out[mask], idldata.yobs_out[mask], color='k')
                ax_i.plot(idldata.xobs_out[mask],idldata.sol_conv[mask], color='r', ls='--')
                ax_i.set_title(line_name)
                ax_i.tick_params(direction='in',top='on')

            pp.savefig(fig); plt.close(fig)

        bar.update(i)

    pp.close()

    bar.finish()

    # Saving the results:
    names = ['ID']

    # Add parameter column names
    for i in range(len(param_lst)):
        names += ['l_'+param_lst[i],param_lst[i],param_lst[i]+'_eUP',param_lst[i]+'_eDW']

    # Add model column name
    names += ['Model_name']

    # Add the photometric column names
    names += ['BC','(B-V)0']

    output = Table(rows=data_rows, names=(names))

    full_path = maindir + 'tables/MAUI_results.' + format

    if format == 'ascii':
        format += '.fixed_width_two_line'

    output.write(full_path, format=format, overwrite=True)

    return 'DONE'


def gen_gridlim(models_dir=mauidir+'MODELS/'):

    '''
    Function to generate fits tables for each MAUI grid with the limits for each
    parameter.

    Parameters
    ----------
    models_dir : str, optional
        Enter the directory where the model files are.

    Returns
    -------
    Nothing but the fits are generated.
    '''

    param_dic = {
    'Teff' : ('Teff_UP','Teff_DW'),
    'logg' : ('lgf_UP','lgf_DW'),
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
                    data_row.extend([np.nan,np.nan]); continue

                idx = parameters.index(par_name)
                if not parameters[0] == 'Teff': print('WARNING: lgf will be wrong!')

                if par_name == 'Teff':
                    data_row.extend(
                    [1e-4*idldata.param[idx].max(),
                    1e-4*idldata.param[idx].min()])

                elif par_name == 'logg':
                    # NOTE: with this we change from logg to loggf to compare with
                    # the solution files from MAUI.
                    data_row.extend([
                    5.39-(4*np.log10(idldata.param[0])-idldata.param[idx]-10.61).min(),
                    5.39-(4*np.log10(idldata.param[0])-idldata.param[idx]-10.61).max()])

                else:
                    data_row.extend([idldata.param[idx].max(),idldata.param[idx].min()])

            data_rows.append(tuple(data_row))

    names = ['Model_name']
    for i in param_dic: names = names + [j for j in param_dic[i]]

    output = Table(rows=data_rows, names=(names))
    output.write(maindir + 'tables/MAUI_grid_limits.fits', format='fits', overwrite=True)

    return('DONE')


def gen_stars_in_grids(table='IACOB_O9BAs_SNR20.fits', table_results='MAUI_results.fits'):

    '''
    Function to generate MAUI-input txt lists containing the stars that lie within the
    boundaries of each MAUI grid.

    Parameters
    ----------
    table : str, optional
        Table containing the input stars

    table_results : str, optional
        Table containing the MAUI results for the stars in the previous table.

    Returns
    -------
    Nothing, but the output lists are created.
    '''

    import matplotlib.path as mpath

    table = findtable(table)
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

        # QUITAR, ES SOLO PARA SERGIO DE CARA A ORDENAR POR FWHM Hb
        RVEW = findtable('RVEWFWs_O9BAs.fits')
        table_red=join(table_red,RVEW,keys='ID')
        table_red['3414'] = table_red['FW34Hb']-table_red['FW14Hb']
        table_red.sort('3414')
        table_red.reverse()

        maui_input(table=table_red, output_name=name, ascii_0=False)


class idl():
    def __init__(self, idlfile):

        '''
        Parameters
        ----------
        idlfile : str
            Enter the input spectrum full path to the .idl file.
        '''

        idldata = readsav(idlfile)
        #for i in idldata.keys():
        #    try: print(i,len(idldata[i]))
        #    except: print(i)

        self.filename = idldata.aa[0][3].decode()
        self.id = self.filename.split('_')[0]
        self.resolution = int(self.filename.split('_V')[-1][0:5])
        self.gridname = idldata.modelgridname.decode()
        self.synwave = idldata.xx_mod
        self.synflux = idldata.spec_prim
        self.dx = (idldata.xx_mod[-1]-idldata.xx_mod[0])/len(idldata.xx_mod)

        self.vsini = idldata.obsdat.spectrum[0].VSINI[0]
        self.vmac = idldata.obsdat.spectrum[0].MACRO[0]


def gen_synthetic(save_dir='server', lwl=3900, rwl=5080):

    '''
    Function to generate the .dat files of the synthetic spectra generated by
    MAUI in all the .idl files found in the SOLUTION folder.

    Parameters
    ----------
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

    solution_dir = mauidir + 'RESULTS_BSGS_202101/SOLUTION/'

    if save_dir == 'local':
        save_dir = datadir + 'ASCII/Synthetic_MAUI/'
        
    elif save_dir == 'server':
        save_dir = '/net/nas/proyectos/hots/masblue/obs_iac/spec_opt/IACOB_DB/ASCII/SYNTHETIC/'

    for file in os.listdir(solution_dir):
        if not file.startswith('._') and file.endswith('.idl'):
            idlspec = idl(solution_dir + file)

            star_db = spec(idlspec.id, SNR='best')

            if idlspec.filename != star_db.filename[:-5]:
                print('\nWARNING: %s does not match with best spectrum available.'
                % idlspec.filename[:-5])

            else:
                #idlspec.filename = idlspec.filename.replace(str(idlspec.resolution),'85000')
                new_idlspec = '%s_red%i.dat' % (idlspec.filename,grids_dic[idlspec.gridname][2])
                np.savetxt(save_dir + new_idlspec, np.c_[idlspec.synwave,idlspec.synflux],
                    fmt=('%.4f','%.6f'))

                star_idl = spec(new_idlspec, txt=True)
                star_idl.txtwaveflux(lwl, rwl)
                #plt.plot(star_idl.wave, star_idl.flux, 'r', lw=.3) # plot to check
                star_idl.degrade(profile='rotmac', vsini=idlspec.vsini, vmac=idlspec.vmac)
                #plt.plot(star_idl.wave, star_idl.flux, 'g', lw=.3) # plot to check

                np.savetxt(save_dir + new_idlspec, np.c_[star_idl.wave,star_idl.flux],
                    fmt=('%.4f','%.6f'))
