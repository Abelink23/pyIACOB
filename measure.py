from RV import *


def measure(lines, table, output_table, RV0lines='rv_Bs.lst', RV0fun='g', RV0tol=150,
    ewcut=10, snrcut=100, tol=100, redo='n'):

    '''
    Function to interactively calculate and store radial velocity, equivalent
    width and full width at half maximum of stars for SiIII and Hb lines.
    Note: input table must contain columns "ID" and "SpC".

    Parameters
    ----------
    lines : str, list
        Enter the wavelenght(s) of the line(s) to fit, either in a coma-separated
        string, or in a .txt/.lst file containing the lines.

    table : str
        Name of the input table contaning the list of stars to analyze.

    output_table : str
        Name of the output (new) table contaning the results.

    RV0lines : str, list
        Enter the wavelenght(s) of the line(s) to fit, either in a coma-separated
        string, or in a .txt/.lst file containing the lines.

    RV0tol : int, optional
        Tolerance for a line to be considered for the RV0 calculation.

    ewcut : float, optional
        EW threshold value for a line to be considered as detected. Default is 10.

    snrcut : float/int, optional
        Cut in the SNR of a line window for the SNR to be replaced by 100. It is then used
        to prevent the line properties to be exported if 3/ snr(line) > depth(line)

    redo : str, optional
        Coma separated string with the list of stars for which repeat the analysis.

    Other parameters : optional
        See help for see spec and spec.fitline

    Returns
    -------
    Nothing, but the output table with the results is created.
    '''

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
                ['ID','Ref_file','SNR_B','SNR_V','SNR_R','RV0']
                + [j+i for i in [str(round(k)) for k in lines] for j in ['RV_','EW_','FW_','dep_','snr_']]
                + ['RV_Hb','EW_Hb','FW_Hb','FW14_Hb','FW34_Hb','dep_Hb','gamma_Hb']),
            dtype = (
                ['S16','S50','int64','int64','int64','float64']
                + ['float64','int64','float64','float64','int64']*len(lines)
                + ['float64','int64','float64','float64','float64','float64','float64'])
            )

    quit = ''
    for source in table:
        if quit == 'q': break

        if 'Ref_file' in source.colnames:
            id = source['Ref_file']
        else:
            id = source['ID']

        if 'SpC' in source.colnames:
            spt  = source['SpC']
        else:
            spt = ''

        if source['ID'] in [i.strip() for i in output['ID']]: continue

        skip = input("%s (%s) - Hit return to continue, type 's' to skip: " % (id,spt))

        next = 'n'
        while next == 'n':

            if skip == 's': break

            star = spec(id, SNR='bestMF')

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

            plt.close()

            star.rv0 = RV0(RV0lines, star.filename, func=RV0fun, ewcut=30, tol=RV0tol)
            star.waveflux(min(lines)-30, max(lines)+30) # PONER MIN MAX EN FUNCION DE LOS LIM DE LINES
            star.cosmic()
            star.plotspec(4510,4600, lines='35-10K')

            input(); plt.close()

            T_source = Table(
                [[star.id_star],[star.filename],[snr_b],[snr_v],[snr_r],[round(star.rv0, 2)]],
                names=('ID','Ref_file','SNR_B','SNR_V','SNR_R','RV0'))
            for line in lines:
                fit = star.fitline(line, width=wid, tol=tol, func=fun, plot=True)

                RV = round(fit['RV_kms']+star.rv0, 2)
                EW = fit['EW']
                FW = fit['FWHM']
                dep = fit['depth']
                snr = fit['snr']

                if EW != None:
                    if EW < ewcut:
                        EW = FW = np.nan

                if snr > snrcut:
                    snr_min = 100
                else: snr_min = snr

                if 3/snr_min > dep:
                    RV = EW = FW = dep = snr = np.nan

                for par,val in zip(['RV_','EW_','FW_','dep_','snr_'],[RV,EW,FW,dep,snr]):
                    T_source[par+str(round(line))] = val

            next = input("Type 'n' to repeat, hit return to move to the Hb line. ")
            plt.close('all')

        next = 'n'
        while next == 'n':

            if skip == 's': break

            #=======================================================================
            print('\nAnalyzing H beta line...\n')

            star.waveflux(4801,4921); star.cosmic()
            star.plotspec(4821,4901, lines='35-10K')

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

            plt.close()

            fit = star.fitline(4861.325, width=wid, func=fun, iter=1, info=True, outfit=True, plot=True)

            RV = round(fit['RV_kms']+star.rv0, 2)
            EW = fit['EW']
            FW = fit['FWHM']
            dep = fit['depth']
            gam = round(fit['gamma'], 3)

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
                print('Line could not be fitted...')
                T_source['FW14_Hb'] =  T_source['FW34_Hb'] = T_source['gamma_Hb'] = np.nan

            next = input("\nRepeat Hb / continue to the next star / save and exit ['n'/''/'q']: ")
            plt.close()

            if next == 'n':
                continue
            else:
                output = vstack([output,T_source])

                if next == 'q': quit = next

    output.write(maindir+'tables/'+output_table, format='fits', overwrite=True)

    return 'DONE'


def auto_measure(lines, table='IACOB_new_N+M_ToDo.fits', output_table='new_RVEWFW.fits',
    func='g', width=20, ewcut=10, tol=150, txt=False, table_fit_param=False):

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

        star = spec(id, SNR='bestMF', txt=txt)

        fig = plt.figure(figsize=(12,10))
        fig.suptitle(id, fontsize=9)

        row_data = []; nplot = 1
        for line,element in zip(lines,elements):

            print('\nAnalyzing %s...\n' % (element+' - '+str(line)))

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


# %% ===========================================================================
#''' Scrip to show the stars for which new spectra with higher SNR is available '''
#table_old = findtable('IACOB_O9-B7_SNR20.fits')
#table_new = findtable('IACOB_O9BAs_SNR20.fits')
#
#for row in table_new:
#    snr_old = table_old[table_old['ID']==row['ID']]['SNR_best']
#    if row['SNR_best'] > snr_old:
#        print(str(row['ID']).strip(),row['SNR_best'],float(snr_old))
