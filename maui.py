from spec import *
from scipy.io.idl import readsav

grids_dic = {
'all': ('Grids coverage','dodgerblue',0,
[[4.543,4.290,4.290,4.146,4.146,4.543,4.543],[2.391,2.391,3.092,3.092,4.391,4.391,2.391]]),
'nlte_10.1.6_SOLAR_expoclump_2019-10-24': ('BSgs_CNOSiMg_old','b',1,
[[4.190,4.477,4.477,4.190,4.190],[3.785,3.785,4.391,4.391,3.785]]),
'nlte_10.1.6_bdwarfs_SOLAR_2020-01-29': ('BDws_CNOSIMg_old','orange',2,
[[4.290,4.543,4.543,4.290,4.290],[2.391,2.391,3.889,3.889,2.391]]),
'nlte_10.4.7_OB.Sg_SOLAR_2021-01-23': ('OBSgs_hot_NOSi_new','g',3,
[[4.399,4.544,4.544,4.399,4.399],[3.488,3.488,4.386,4.386,3.488]]),
'nlte_10.4.7_late.bsgs_SOLAR_expoclump_NOSi.djl_2021-02-06': ('BSgs_cool_NOSi_new','r',4,
[[4.146,4.322,4.322,4.146,4.146],[3.092,3.092,4.391,4.391,3.092]]),
'astar2013_SOLAR_2_LMC_4_grid_2019-10-24_2019-10-24' : ('ASgs_new','purple',5,
[[3.900,4.114,4.114,3.900,3.900],[3.142,3.142,4.292,4.292,3.142]])
}

def mauipath(path=None):
    '''
    Function to set the main directory of MAUI.

    Parameters
    ----------
    path : str, optional
        If 'default' or 'def', it will choose the default path (see below).

    Returns: Selected MAUI main directory path.
    '''

    if path in ['def','default']:
        if platform.system() == 'Darwin':
            defmainpath = '/Users/abelink/Documents/MAUI/'
        elif platform.uname().node == 'msi':
            defmainpath = '/home/abelink/Documents/MAUI/'
        elif 'iac.es' in platform.uname().node:
            defmainpath = '/net/nas/proyectos/hots/masblue/maui2021/RESULTS_BSGS_202101/'

        mainpath = defmainpath

    elif path == None:
        mainpath = input('Working directory path (default is %s) : ' %defmainpath)
        if mainpath == '': mainpath = defmainpath

    else: mainpath = path

    return mainpath

mauidir = mauipath('def')


def maui_input(table='IACOB_O9BAs_SNR20.fits',output_name='MAUI_verX',RV0tol=200,ascii_0=False):
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
        Enter the input table contaning a column 'Name' with the name of the stars.

    output_name : str, optional
        Enter the name for the output table. Default is 'MAUI_verX'.

    RV0tol : int, optional
        Enter the input radial velocity tolerance for the radial velocity correction.

    ascii_0 : boolean, optional
        If True, ascii files will be created for each of the input sources.
        Default is False.

    Returns: Nothing but the MAUI input file is generated.
    '''

    table = findtable(table) # file where star names and quality flags are
    table_REF = findtable('RVEWFWs_O9BAs.fits') # file where RVs, EWs and FWs are
    table_IB = findtable('IB_results_ver5.txt') # file where vsini and vmac are
    results = findtable('MAUI_results.fits') # file with output from MAUI

    maui_txt = open(maindir+'lists/%s.txt' % output_name,'a')
    maui_txt.write(
        "{:<40}".format('starname')+"{:<48}".format('filename')+"{:<6}".format('vrad')+\
        "{:<6}".format('vsini')+"{:<7}".format('evsini')+"{:<5}".format('vmac')+\
        "{:<7}".format('evmac')+"{:<7}".format('R')+"{:<4}".format('SNR')+\
        ' ;# SpC       FW34-14 SiIII SiII l lTef l lgf  Grid\n')

    quit = ''
    for row in table:

        ascii = ascii_0

        if quit == 'quit': break

        name = row['Name'].strip()

        # Filer based on properties from the main table:
        if 'SB' in table.columns and 'SB2' in row['SB']: continue
        if 'CHb' in table.columns:
            if 'Em' in row['CHb']:
                if not 'Em(p)' in row['CHb']: continue
            if 'PCyg' in row['CHb']: continue
        if 'QIB' in table.columns and row['QIB'] < 2: continue

        match_REF = table_REF[[i.strip()==name for i in table_REF['Name']]]
        match_IB = table_IB[table_IB['Name']==name]

        if len(match_REF) == 0 or len(match_IB) == 0: continue

        # Filter based on Si lines properties:
        if row['QSiII']<3 or match_REF['EWSiII']<50 or np.isnan(match_REF['EWSiII']) or\
            match_REF['depSiII']<3/match_REF['SNR_B']: SiIIFG = 0
        else: SiIIFG = 1

        if row['QSiIII']<3 or match_REF['EWSiIII1']<50 or np.isnan(match_REF['EWSiIII1']) or\
            match_REF['depSiIII1']<3/match_REF['SNR_B']: SiIIIFG = 0
        else: SiIIIFG = 1

        if SiIIIFG == 0 and SiIIFG == 0:
            print('No SiIII/SiII found for %s\n' % name); continue

        star = spec(name,SNR='best')

        if match_IB['filename'][0] != star.file_name:
            print('Warning: Different files from best SNR and from IB results for %s' % name)
            print(star.file_name,' vs ',match_IB['filename'][0])
            if ascii_0 == True:
                do_file = input('Which ascii do you want to create 1 or 2: ')
                if int(do_file) == 1: filename = star.file_name
                elif int(do_file) == 2: filename = match_IB['filename'][0]

        if ascii_0 == True and not search(filename[:-5]+'_RV.ascii',\
        os.path.expanduser('~')+'/Documents/MAUI/ASCII/') is None: ascii = False

        # ----------------------------------------------------------------------
        # Extra information appended to the end of each row:
        match_results = results[results['Name']==name]
        if len(match_results) == 0:
            l_Tef = logTf = l_lgf = loggf = grid = 0
        else:
            l_Tef = match_results['l_Teff'][0]; logTf = 4+np.log10(match_results['Teff'][0])
            l_lgf = match_results['l_lgf'][0]; loggf = match_results['lgf'][0];
            grid = grids_dic[match_results['Model_name'][0]][0]
        # ----------------------------------------------------------------------

        if ascii == True:

            from RV import RV0

            # If RV SiIII is good enough, it uses it for the offset:
            if abs(match_REF['RVSiIII1']-match_REF['RVHb']) < 10:
                star.offset = float(match_REF['RVSiIII1'])*float(4500)*1000/cte.c
                star.waveflux(); print('RV good enough.\n')

            # Otherwise it calculates the offset with the RV0 program:
            else:
                next = 'n'
                while next == 'n':

                    skip = input('%s - Hit return to continue, type "s" to skip: ' % name)
                    if skip == 's': break

                    if match_REF['SpT_code'] <= 2.5: star.plotspec(4530,4590)
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

                    star.offset = RV0(spt_list,star.spectrum,ewcut=30,width=wid,tol=RV0tol,func=fun)
                    star.waveflux() # Applies the offset
                    #star.cosmic(sigclip=0.005)

                    if match_REF['SpT_code'] <= 2.5: star.plotspec(4530,4590,poslines='OB')
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
            "{:<40}".format(star.file_name[:-5])+\
            "{:<48}".format(star.file_name[:-5]+'_RV.ascii')+"{:<6}".format('0.0d0')+\
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

    if ascii == True: print('Remember to move the new ascii into "ASCII_ARCHIVE" folder')

    return 'DONE'


def gen_stars_in_grids(table='IACOB_O9BAs_SNR20.fits'):

    import matplotlib.path as mpath

    log_Teff = np.asarray(4+np.log10(results['Teff']))
    log_LLsol = np.asarray(5.39-results['lgf'])

    points = np.column_stack([log_Teff,log_LLsol])

    # Change ...grids_dic][1:] for individual grids
    for name,_,_,box in [grids_dic[i] for i in grids_dic][1:]:
        verts = np.array([box[0],box[1]]).T
        path = mpath.Path(verts)
        inout = path.contains_points(points)
        log_Teff_in,log_LLsol_in = points[inout].T

        results_in = results[path.contains_points(points)]['Name']
        table_red = table[[i['Name'].strip() in results_in for i in table]]

        maui_input(table=table_red,output_name=name,ascii_0=False)


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

        self.file_name = idldata.aa[0][3].decode()
        self.name_star = self.file_name.split('_')[0]
        self.resolution = int(self.file_name.split('_V')[-1][0:5])
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

    Returns: Nothing but the ascii .dat files are generated.
    '''

    solution_dir = mauidir+'SOLUTION/'

    if save_dir == 'local':
        save_dir = datadir+'ASCII/Synthetic_MAUI/'
    elif save_dir == 'server':
        save_dir = '/net/nas/proyectos/hots/masblue/obs_iac/spec_opt/IACOB_DB/ASCII/SYNTHETIC/'

    for file in os.listdir(solution_dir):
        if not file.startswith('._') and file.endswith('.idl'):
            idlspec = idl(solution_dir+file)

            star_db = spec(idlspec.name_star,SNR='best')

            if idlspec.file_name != star_db.file_name[:-5]:
                print('\nWARNING: %s does not match with best spectrum available.'
                % idlspec.file_name[:-5])

            else:
                #idlspec.file_name = idlspec.file_name.replace(str(idlspec.resolution),'85000')
                new_idlspec = '%s_red%i.dat' % (idlspec.file_name,grids_dic[idlspec.gridname][2])
                np.savetxt(save_dir+new_idlspec,np.c_[idlspec.synwave,idlspec.synflux],
                           fmt=('%.4f','%.6f'))

                star_idl = spec(new_idlspec,txt=True)
                star_idl.txtwaveflux(lwl,rwl)
                #plt.plot(star_idl.wave,star_idl.flux,'r',lw=.3) # plot to check
                star_idl.degrade(profile='rotmac',vsini=idlspec.vsini,vmac=idlspec.vmac)
                #plt.plot(star_idl.wave,star_idl.flux,'g',lw=.3) # plot to check

                np.savetxt(save_dir+new_idlspec,np.c_[star_idl.wave,star_idl.flux],
                           fmt=('%.4f','%.6f'))


def gen_table(tables_dir='server', input_table='MAUI_ver10.txt', check_best=True,
              grids_table='MAUI_grid_limits.fits', format='fits'):
    '''
    Function to generate a table with the results from MAUI given an input table
    containing the name of the stars and filename to search in the MAUI-SOLUTION
    directory.

    Parameters
    ----------

    tables_dir : str, optional
        Enter the directory where to locate the input_table and output fits.

    input_table : str, optional
        Name of the input table contaning the list of stars to search.

    check_best : boolean, optional
        True if each spectra from the input_table is checked against the best
        spectrum in the database. Default is True.

    grids_table : str, obtional
        Name of the table containing the limits of the grids in MAUI.

    format : str, optional
        Enter the output format for the table: 'fits' (default), 'ascii' or 'csv'.

    Returns: Nothing but the output table with the MAUI results is generated.
    '''

    solution_dir = mauidir+'SOLUTION/'

    if tables_dir == 'local':
        tables_dir = maindir+'tables/'
    elif tables_dir == 'server':
        tables_dir = '/net/nas/proyectos/hots/adeburgos/tables/'

    table = findtable(input_table)
    grids = findtable(grids_table)

    param_lst = ['Teff','lgf','He','Micro','logQs','beta','C','N','O','Mg','Si','S','Fe','Ti','fcl','vcl']

    bar = pb.ProgressBar(maxval=len(table),
                         widgets=[pb.Bar('=','[',']'),' ',pb.Percentage()])
    bar.start()

    data_rows = []
    for row,i in zip(table,range(len(table))):
        file_name = row['filename']

        match = []
        for file in os.listdir(solution_dir):
            if file.endswith('.idl') and \
               file.split('_sqexp_mat1_')[1][:-14]+'RV.ascii' == file_name:
                   match.append(solution_dir+file)

        if len(match) == 0:
            print('\nWARNING: No .idl file found for %s. Continuing...' % file_name)
            continue

        elif len(match) > 1:
            match = sorted(match, key=lambda x: int(x[-14:-4].replace('-','')),reverse=True)

        idldata = readsav(match[0])

        file_name = idldata.aa[0][3].decode()
        name_star = file_name.split('_')[0]
        grid  = grids[grids['Model_name'] == idldata.modelgridname.decode()]

        if check_best == True and file_name != spec(name_star,SNR='best').file_name[:-5]:
            print('\nWARNING: %s does not match with best spectrum available.'
            % file_name[:-5])

        data_row = []; data_row.extend([name_star])
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

        data_row.extend([grid['Model_name'][0]])
        data_rows.append(tuple(data_row))

        bar.update(i)
        time.sleep(0.1)

    bar.finish()

    names = ['Name']
    for i in range(len(param_lst)):
        names += ['l_'+param_lst[i],param_lst[i],param_lst[i]+'_eUP',param_lst[i]+'_eDW']
    names += ['Model_name']

    output = Table(rows=data_rows,names=(names))

    full_path = tables_dir+'MAUI_results'+'.'+format
    if format == 'ascii': format += '.fixed_width_two_line'
    output.write(full_path,format=format,overwrite=True)

    return 'DONE'


def gen_gridlim(tables_dir='local'):
    '''
    Function to generate fits tables for each MAUI grid with the limits for each
    parameter.

    Parameters
    ----------

    tables_dir : str, optional
        Enter the directory where to save the output fits files.

    Returns: Nothing but the fits are generated.
    '''

    models_dir = mauidir+'MODELS/'

    if tables_dir == 'local':
        tables_dir = maindir+'tables/'
    elif tables_dir == 'server':
        tables_dir = '/net/nas/proyectos/hots/adeburgos/tables/'

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

            idldata

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

    output = Table(rows=data_rows,names=(names))
    output.write(tables_dir+'MAUI_grid_limits.fits',format='fits',overwrite=True)

    return('DONE')
