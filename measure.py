from rv import *

from binarity import *


def measure(lines, table, output_table, colnames='lambda', rv_lines='rv_Bs.lst', rv_func='g', rv_tol=150,
    ewcut=10, tol=100, orig='IACOB', redo='n', show_plot=False, do_pdf=True):

    '''
    Function to interactively calculate and store radial velocity, equivalent
    width and full width at half maximum of stars a given list of input lines.

    Note: Feature/keyword to find binary stars is not implemented yet. Only for the Hb.

    Parameters
    ----------
    lines : str, list
        Enter the wavelenght(s) of the line(s) to fit, either in a coma-separated
        string, or in a .txt/.lst file containing the lines.

    table : str
        Name of the input table contaning the list of stars to analyze.

    output_table : str
        Name of the output (new) table contaning the results.

    colnames : str, optional
        Use 'lambda' if the table has the wavelenghts from the list of lines.
        Use 'label' if the table has the line labels from the list of lines.
        
    rv_lines : str, list
        Enter the wavelenght(s) of the line(s) to fit, either in a coma-separated
        string, or in a .txt/.lst file containing the lines.

    rv_func : str, optional
        Choose the function to fit the lines used for the initial RV:
        'g' Gaussian (default); 'l' Lorentzian; 'v' Voigt; 'r' Rotational.

    rv_tol : int, optional
        Tolerance for a line to be considered for the RV0 calculation.

    ewcut : float, optional
        EW threshold value for a line to be considered as detected. Default is 10.

    tol : int, optional
        Sets the tolerance [km/s] to shifting the spectrum in order to fit the lines 
        after being corrected from radial velocity. Default is 100.

    orig : str, optional
        If 'IACOB', it assumes that the spectrum comes from the IACOB database.
        If 'txt', it assumes that the spectrum comes from a two-columns file with
        wavelenght and flux with no header in the file. 
        If 'syn', it assumes that the spectrum comes from a synthetic spectrum.
        Default is IACOB.

    redo : str, optional
        Coma separated string with the list of stars for which repeat the analysis.

    show_plot : bool, optional
        If True, the plot of the spectrum and the fit is shown for each line.
        Default is False.

    do_pdf : bool, optional
        If True, the plots of the spectrum and the fittings are saved in a pdf file.

    Returns
    -------
    Nothing, but the output table with the results is created.
    '''

    if colnames not in ['lambda','label']:
        print('Bad input for "colnames" parameter. Exitting...')
        return None
    elif colnames == 'lambda':
        line_names = [str(round(i,1)) for i in findlines(lines)[0]]
    elif colnames == 'label':
        line_names = findlines(lines)
        if len(line_names) == 1:
            print('Lines file must include labels for the lines. Exitting...')
            return None
        else:
            line_names = line_names[1]

    lines = findlines(lines)[0]

    table = findtable(table)

    #===========================================================================
    output = findtable(output_table)

    if output is not None:
        # Remove entries to be re-done
        if redo != 'n':
            try:
                redo = redo.split(',')
            except:
                print('Bad input for "redo" parameter. Exitting...')
                return None

            for id in redo:
                output = output[output['ID'] != id]

    else:
        print('Creating new table: %s\n' % output_table)

        output = Table(
            names = (
                ['ID','Ref_file','SNR_B','SNR_V','SNR_R','RV0','eRV0']
                + [j+i for i in line_names for j in ['RV_','EW_','FW_','dep_','snr_']]),
            dtype = (
                ['S16','S50','int64','int64','int64','float64','float64']
                + ['float64','int64','float64','float64','int64']*len(line_names))
            )

    # To create the pdf file
    if do_pdf == True:

        from matplotlib.backends.backend_pdf import PdfPages

        if not os.path.exists(maindir + 'plots/measure/'):
            os.makedirs(maindir + 'plots/measure/')

        pdf_fitting = PdfPages(maindir + 'plots/measure/measure_lines_%s.pdf' % time.strftime('%Y%m%d_%H%M%S'))

        plt.rcParams.update({
            'xtick.labelsize' : 6,
            'ytick.labelsize' : 6,
            'axes.titlesize' : 6,
            'axes.titlepad' : 3
            })

        # To define the number of rows and columns for the plots
        nrows, ncols = even_plot(len(lines))

        # Avoid showing the spectrum for each line when is already saved in pdf
        show_plot = False


    quit = ''
    for source in table:
        if quit == 'q': break

        if 'Ref_file' in source.colnames:
            id = source['Ref_file']
            if source['Ref_file'] in [i.strip() for i in output['Ref_file']]:
                continue
        else:
            id = source['ID']
            if source['ID'] in [i.strip() for i in output['ID']]:
                continue

        if 'SpC' in source.colnames:
            spt  = source['SpC']
        else:
            spt = ''

        skip = input("%s (%s) - Hit return to continue, type 's' to skip: " % (id,spt))

        next = 'n'
        while next == 'n':

            if skip == 's': break

            star = spec(id, SNR='bestHF', orig=orig)

            snr_b = star.snrcalc(zone='B')
            snr_v = star.snrcalc(zone='V')
            snr_r = star.snrcalc(zone='R')

            #===================================================================
            print('\nAnalyzing Si III triplet...\n')
            star.plotspec(4510,4600)

            fun = '-'
            while fun not in ['g','r','vr_Z','vrg_Z']:
                fun = input('Choose function to fit between g,r,vr_Z,vrg_Z (default is g): ')
                if fun == '': fun = 'g'
            wid = '-'
            while type(wid) is not float:
                wid = input('Choose the initial width in angstroms (default is 15): ')
                if wid == '': wid = 15.
                else:
                    try: wid = float(wid)
                    except: wid = '-'

            plt.close('all')

            star.rv0, eRV0 = RV0(rv_lines, star.filename, func=rv_func, ewcut=30, tol=rv_tol, orig=orig)

            star.waveflux(min(lines)-30, max(lines)+30) # PONER MIN MAX EN FUNCION DE LOS LIM DE LINES
            star.cosmic()

            star.plotspec(4510,4600, lines='35-10K')

            input('Hit return to continue...')
            plt.close('all')

            T_source = Table(
                [[star.id_star],[star.filename],[snr_b],[snr_v],[snr_r],[round(star.rv0, 2)],[round(eRV0, 2)]],
                names=('ID','Ref_file','SNR_B','SNR_V','SNR_R','RV0','eRV0'))

            for i,line in zip(range(len(lines)),lines):

                if do_pdf == True and (i % 36) == 0:
                    if len(lines) <= 36:
                        fig, ax = plt.subplots(nrows, ncols, figsize=(13,8))
                    else:
                        fig, ax = plt.subplots(6, 6, figsize=(13,8))
                    axs = ax.flatten()

                fit = star.fitline(line, width=wid, tol=tol, func=fun, plot=show_plot, outfit=do_pdf)

                RV = round(fit['RV_kms']+star.rv0, 2)
                EW = fit['EW']
                FW = fit['FWHM']
                dep = fit['depth']
                snr = fit['snr']

                if EW != None and EW < ewcut:
                    EW = FW = np.nan

                # Threshold in the SNR of a line window for the SNR to be replaced by 100. This is used
                # to prevent the line properties to be exported if 3/snr(line) > depth(line)
                snrcut = 100 
                if snr > snrcut:
                    snr_min = snrcut
                else:
                    snr_min = snr

                if 3/snr_min > dep:
                    RV = EW = FW = dep = snr = np.nan

                for par,val in zip(['RV_','EW_','FW_','dep_','snr_'],[RV,EW,FW,dep,snr]):
                    T_source[par+line_names[i]] = val

                if do_pdf == True:

                    c_fit = 'g'
                    c_spec = 'k'

                    if fit['sol'] == 0:
                        c_spec = 'maroon'

                    elif fit['sol'] == 1:
                        if abs(fit['RV_kms']) > 0.5*tol:
                            c_fit = 'orange'
                        elif fit['EW'] < ewcut:
                            c_fit = 'gold'
                        elif 3/snr_min > dep:
                            c_fit = 'maroon'
                        else:
                            c_fit = 'g'

                    axs[i % 36].set_title(line_names[i])
                    axs[i % 36].tick_params(direction='in', top='on', right='on')

                    if not fit['sol'] == 0 and 'wave' in fit.keys():
                        axs[i % 36].plot(fit['wave'], fit['flux_norm'], c=c_spec, lw=0.5)
                        axs[i % 36].plot(fit['wave'], fit['flux_fit'], c=c_fit, ls='--', lw=1.5)

                    if i % 36 == 35 or i == len(lines)-1:
                        fig.suptitle(star.id_star + ' -- Filename: ' + star.filename + ' -- Lines for rv corrextion: ' + rv_lines, fontsize=8)

                        # Delete the unused axes
                        [fig.delaxes(axs[i]) for i in np.arange(i % 36 + 1, len(axs))]

                        fig.tight_layout()
                        fig.subplots_adjust(wspace=0.2, hspace=0.3)

                        pdf_fitting.savefig(fig)
                        plt.close(fig)

            next = input("\nRepeat / continue to the next star / save and exit ['n'/''/'q']: ")
            plt.close('all')

            if next == 'n':
                continue
            else:
                output = vstack([output,T_source])

                if next == 'q':
                    quit = next

    if do_pdf == True:
        pdf_fitting.close()

    output.write(maindir+'tables/'+output_table, format='fits', overwrite=True)

    return 'DONE'



def measure_Hb(table, output_table, rv_lines='rv_Bs.lst', rv_func='vrg_H', rv_tol=150, binarity=False, orig='IACOB', redo='n'):

    '''
    Function to interactively calculate and store radial velocity, equivalent
    width and full width at half maximum 3/4 and 1/4 from continuum for Hb line.

    Parameters
    ----------
    table : str
        Name of the input table contaning the list of stars to analyze.

    output_table : str
        Name of the output (new) table contaning the results.

    rv_lines : str, list
        Enter the wavelenght(s) of the line(s) to fit, either in a coma-separated
        string, or in a .txt/.lst file containing the lines.

    rv_func : str, optional
        Choose the function to fit the lines used for the initial RV:
        'g' Gaussian; 'l' Lorentzian; 'v' Voigt; 'r' Rotational; 
        'vr_H' Voigt profile with rotation;
        'vrg_H' (default) Voigt profile with rotation and extra gaussian.

    rv_tol : int, optional
        Tolerance for a line to be considered for the RV0 calculation.

    binarity : bool, optional
        If True, the function will also make use of the binarity module to help 
        you to identify possible binary stars. See binarity.py for more info.
        Default is False.

    orig : str, optional
        If 'IACOB', it assumes that the spectrum comes from the IACOB database.
        If 'txt', it assumes that the spectrum comes from a two-columns file with
        wavelenght and flux with no header in the file. 
        If 'syn', it assumes that the spectrum comes from a synthetic spectrum.
        Default is IACOB.

    redo : str, optional
        Coma separated string with the list of stars for which repeat the analysis.

    Returns
    -------
    Nothing, but the output table with the results is created.
    '''

    table = findtable(table)

    #===========================================================================
    output = findtable(output_table)

    if output is not None:
        # Remove entries to be re-done
        if redo != 'n':
            try:
                redo = redo.split(',')
            except:
                print('Bad input for "redo" parameter. Exitting...')
                return None

            for id in redo:
                output = output[output['ID'] != id]

    else:
        print('Creating new table: %s\n' % output_table)

        output = Table(
            names = (
                ['ID','Ref_file','SNR_B','SNR_V','SNR_R','RV0','eRV0']
                + ['RV_Hb','EW_Hb','FW_Hb','FW14_Hb','FW34_Hb','dep_Hb','gamma_Hb']),
            dtype = (
                ['S16','S50','int64','int64','int64','float64','float64']
                + ['float64','int64','float64','float64','float64','float64','float64'])
            )

    quit = ''
    for source in table:
        if quit == 'q': break

        if 'Ref_file' in source.colnames:
            id = source['Ref_file']
            if source['Ref_file'] in [i.strip() for i in output['Ref_file']]:
                continue

        else:
            id = source['ID']
            if source['ID'] in [i.strip() for i in output['ID']]:
                continue

        if 'SpC' in source.colnames:
            spt  = source['SpC']
        else:
            spt = ''

        skip = input("%s (%s) - Hit return to continue, type 's' to skip: " % (id,spt))

        if binarity == True:
            findSB(id.split('/')[-1].split('_')[0],SNR=20,degrade=40000,vspace=0)
            plt.figure()

        star = spec(id, SNR='bestHF', orig=orig)

        snr_b = star.snrcalc(zone='B')
        snr_v = star.snrcalc(zone='V')
        snr_r = star.snrcalc(zone='R')

        star.rv0, eRV0 = RV0(rv_lines, star.filename, func=rv_func, ewcut=30, tol=rv_tol, orig=orig)
 
        next = 'n'
        while next == 'n':

            if skip == 's': break

            T_source = Table(
                [[star.id_star],[star.filename],[snr_b],[snr_v],[snr_r],[round(star.rv0, 2)],[eRV0]],
                names=('ID','Ref_file','SNR_B','SNR_V','SNR_R','RV0','eRV0'))

            #=======================================================================
            print('\nPrinting H alpha line to detect winds...\n')

            star.plotspec(6522.80,6602.80, lines='35-10K')

            mngr = plt.get_current_fig_manager()
            x,y,dx,dy = mngr.window.geometry().getRect()
            # Get the current window size and position
            mngr.window.setGeometry(x+1/2*dx, y, dx, dy) #900

            plt.figure()

            #=======================================================================
            print('\nAnalyzing H beta line...\n')

            star.waveflux(4801,4921)
            star.cosmic()

            star.plotspec(4821,4901, lines='35-10K')

            mngr = plt.get_current_fig_manager()
            x,y,dx,dy = mngr.window.geometry().getRect()
            mngr.window.setGeometry(x-1/2*dx, y, dx, dy) # -900

            fun = '-'; iter = 3
            while fun not in ['vr_H','vrg_H']:
                fun = input('Choose function to fit between vr_H/vrg_H (default is vrg_H): ')
                if fun == '': fun = 'vrg_H'; iter = 1

            wid = '-'
            while type(wid) is not float:
                wid = input('Choose the initial width in angstroms (default is 50): ')
                if wid == '': wid = 50.
                else:
                    try: wid = float(wid)
                    except: wid = '-'

            plt.close('all')
            fit = star.fitline(4861.325, width=wid, func=fun, iter=1, fw3414=True, 
                               info=True, outfit=True, plot=True)

            if fit['sol'] != 0:
                RV = round(fit['RV_kms']+star.rv0, 2)
                EW = fit['EW']
                FW = fit['FWHM']
                dep = fit['depth']
                if fit['sol'] != 0:
                    gam = round(fit['gamma'], 3)
                else: gam = np.nan

                for par,val in zip(['RV_Hb','EW_Hb','FW_Hb','dep_Hb','gamma_Hb'],[RV,EW,FW,dep,gam]):
                    T_source[par] = val

                try:
                    wave = fit['wave']; flux_fit = fit['flux_fit']
                    lowval = (max(flux_fit) + 3*min(flux_fit))/4
                    uppval = (3*max(flux_fit) + min(flux_fit))/4

                    for par,val in zip(['FW14_Hb','FW34_Hb'],[lowval,uppval]):
                        medpos = [np.where(flux_fit <= val)[0][value] for value in (0,-1)]
                        try: l_val = np.interp(val,[flux_fit[medpos[0]],flux_fit[medpos[0]-1]],
                                                       [wave[medpos[0]],wave[medpos[0]-1]])
                        except: l_val = wave[medpos[0]]
                        try: r_val = np.interp(val,[flux_fit[medpos[1]],flux_fit[medpos[1]+1]],
                                                      [wave[medpos[1]],wave[medpos[1]+1]])
                        except: r_val = wave[medpos[1]]

                        T_source[par] = round(r_val-l_val, 3)

                except:
                    print('Problem calculating the FW at 1/4 and 3/4 of the Hb line.')
                    T_source['FW14_Hb'] =  T_source['FW34_Hb'] = T_source['gamma_Hb'] = np.nan

            else:
                print('Hb line could not be fitted...')

            next = input("\nRepeat Hb / continue to the next star / save and exit ['n'/''/'q']: ")
            plt.close('all')

            if next == 'n':
                continue
            else:
                output = vstack([output,T_source])

                if next == 'q': quit = next

    output.write(maindir+'tables/'+output_table, format='fits', overwrite=True)

    return 'DONE'



def auto_measure(lines, table='IACOB_new_N+M_ToDo.fits', output_table='new_RVEWFW.fits',
    func='g', width=20, ewcut=10, tol=150, orig='IACOB', table_fit_param=False):

    '''
    Function to automatically calculate and store radial velocity, equivalent
    width and full width at half maximum of stars for input lines.

    Parameters
    ----------
    lines : str, list
        Enter the wavelenght(s) of the line(s) to fit, either in a coma-separated
        string, or in a .txt/.lst file containing the lines.

    table : str
        Name of the input table contaning the list of stars to analyze.

    output_table : str
        Name of the output (new) table contaning the results.

    ewcut : float, optional
        EW threshold value for a line to be considered as detected. Default is 10.

    txt : bolean, optional
        If True, it assumes spectrum from a two-columns file with wavelenght and flux.

    Other parameters : optional
        See help for see spec and spec.fitline

    Returns
    -------
    Nothing, but the output table with the results is created.
    '''

    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages(maindir+'tmp_plots/auto_RVEWFW.pdf')

    lines,elements,_ = findlines(lines)

    table = findtable(table)

    IDs = []; A1 = []; x0 = []; sig = []; sig1 = []; gamma = []; vsini = []; A2 = []; sig2 = []

    '''============================ Create table ============================'''

    columns = ['RV_','EW_','FW_','dep_','snr_','FW14_','FW34_']

    names = ['ID']+([j+i for i in [str(round(k)) for k in lines] for j in columns])
    dtypes =['S%i' % len(sorted(table['ID'],key=len)[-1]+' ')]
    dtypes += ['float64']*len(columns)*len(lines)

    output = Table(names=(names),dtype=(dtypes))

    if func in ['vrg_H','vrg_Z']: iter = 1
    else: iter = 3

    nrows = int(len(lines)/3)
    if len(lines) % 3 != 0.0: nrows += 1
    if len(lines) < 3: ncols = len(lines)
    else: ncols = 3

    bar = pb.ProgressBar(maxval=len(table),
                         widgets=[pb.Bar('=','[',']'),' ',pb.Percentage()])
    bar.start()

    for row,i in zip(table,range(len(table))):

        id = row['ID'].strip()

        #if int(re.findall('[0-9]+',id)[2]) < 100:  # TOREMOVE
        #    width = 10; func = 'g' # TOREMOVE
        #if 100<= int(re.findall('[0-9]+',id)[2]) < 175: # TOREMOVE
        #    width = 15; func = 'vr' # TOREMOVE
        #if int(re.findall('[0-9]+',id)[2]) >= 175 # TOREMOVE
        #    width = 15; func = 'r' # TOREMOVE
        #else: continue

        #if int(re.findall('[0-9]+',id)[2]) >= 175: continue

        star = spec(id, SNR='bestHF', orig=orig)

        fig = plt.figure(figsize=(12,10))
        fig.suptitle(id, fontsize=9)

        row_data = []; nplot = 1
        for line,element in zip(lines,elements):

            print('\nAnalyzing %s...\n' % (str(element)+' - '+str(line)))

            #star.cosmic()

            fit = star.fitline(line, width=width, func=func, iter=iter, outfit=True)

            if fit['sol'] == 0:
                row_data += [np.nan]*7

            else:

                plt.subplot(nrows, ncols, nplot)

                row_data += fit['RV_kms']+star.rv0, fit['EW'], fit['FWHM'], fit['depth'], fit['snr']

                wave,flux_norm,flux_fit = fit['wave'], fit['flux_norm'], fit['flux_fit']

                if table_fit_param == True:
                    IDs.append(id)
                    if func in ['vrg_H','vrg_Z']:
                        x0.append(fit['lam0'])
                        A1.append(fit['A1'])
                        sig1.append(fit['sigma1'])
                        gamma.append(fit['gamma'])
                        vsini.append(fit['vsini'])
                        A2.append(fit['A2'])
                        sig2.append(fit['sigma2'])
                    elif func in ['vr_H','vr_Z']:
                        A1.append(fit['A'])
                        x0.append(fit['lam0'])
                        sig.append(fit['sigma'])
                        gamma.append(fit['gamma'])
                        vsini.append(fit['vsini'])
                    elif func in ['r']:
                        A1.append(fit['A'])
                        x0.append(fit['lam0'])
                        sig.append(fit['sigma'])
                        vsini.append(fit['vsini'])
                    elif func in ['v']:
                        A1.append(fit['A'])
                        x0.append(fit['lam0'])
                        sig.append(fit['sigma'])
                        gamma.append(fit['gamma'])
                    elif func in ['g','l']:
                        A1.append(fit['A'])
                        x0.append(fit['lam0'])
                        sig.append(fit['sigma'])

                lowval = (max(flux_fit) + 3*min(flux_fit))/4
                uppval = (3*max(flux_fit) + min(flux_fit))/4

                for val in lowval,uppval:
                    medpos = [np.where(flux_fit <= val)[0][value] for value in (0,-1)]
                    try: l_val = np.interp(val,[flux_fit[medpos[0]],flux_fit[medpos[0]-1]],
                                                   [wave[medpos[0]],wave[medpos[0]-1]])
                    except: l_val = wave[medpos[0]]
                    try: r_val = np.interp(val,[flux_fit[medpos[1]],flux_fit[medpos[1]+1]],
                                                  [wave[medpos[1]],wave[medpos[1]+1]])
                    except: r_val = wave[medpos[1]]
                    row_data += [round(r_val-l_val,3)]

                plt.plot(wave, flux_norm, 'b', lw=.5)
                plt.plot(wave, flux_fit, 'g', lw=.5)

                plt.title('%s | RV: %d | EW: %d | FWHM: %.2f | SNR: %d' %
                (line,fit['RV_kms']+star.rv0,fit['EW'],fit['FWHM'],fit['snr']), fontsize=7, pad=4)

                plt.tick_params(direction='in', top='on')
                plt.ylim(ymax=1.03, ymin=0.4)

                nplot += 1

        fig.tight_layout()
        pp.savefig(fig)
        plt.close(fig)

        output.add_row(([id]+row_data))

        bar.update(i)

    pp.close()

    output.write(maindir+'tables/'+output_table, format='fits', overwrite=True)

    bar.finish()

    if table_fit_param == True:
        table = Table(); table['ID'] = IDs
        if func in ['vrg_H','vrg_Z']:
            table['A1'] = A1; table['x0'] = x0; table['sig1'] = sig1; table['gamma'] = gamma
            table['vsini'] = vsini; table['A2'] = A2; table['sig2'] = sig2
        else:
            table['A1'] = A1; table['x0'] = x0; table['sig'] = sig
            if func in ['vr_H','vr_Z']:
                table['gamma'] = gamma; table['vsini'] = vsini
            elif func == 'r':
                table['vsini'] = vsini
            elif func == 'v':
                table['gamma'] = gamma
        table.write(maindir+'tables/param_'+output_table, format='fits', overwrite=True)

    return 'DONE'


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


