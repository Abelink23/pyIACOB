from spec import *
from RV import *


def RVEWFW(table='IACOB_O9BAs_SNR20.fits', output_table='O9BAs_RVEWFWs.fits',
    RV0lines='rv_Bs.lst', RV0tol=150, ewcut=10, tol=100, redo='n'):

    '''
    Function to interactively calculate and store radial velocity, equivalent
    width and full width at half maximum of stars for SiIII and Hb lines.
    Note: input table must contain columns "ID" and "SpC".

    Parameters
    ----------
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

    redo : str, optional
        Coma separated string with the list of stars for which next the analysis.

    Other parameters : optional
        See help for see spec and spec.fitline

    Returns
    -------
    Nothing, but the output table with the results is created.
    '''

    table = findtable(table)

    #===========================================================================
    try: # Check if the input table already exist
        output = findtable(output_table)
        # Remove entries to be re-done
        if redo != 'n':
            try: redo = redo.split(',')
            except: print('Bad input for "redo" parameter. Exitting...'); return None
            for id in redo: output = output[output['ID'] != id]
    except:
        format1 = str(table['ID'].info.dtype)[1:]
        format2 = str(table['SpC'].info.dtype)[1:]

        output = Table(names=(\
        'ID','SNR_B','SNR_V','SNR_R',\
        'RVSiIII1','EWSiIII1','FWSiIII1','depSiIII1',\
        'RVSiIII2','EWSiIII2','FWSiIII2','depSiIII2',\
        'RVSiIII3','EWSiIII3','FWSiIII3','depSiIII3',\
        'RVSiII','EWSiII','FWSiII','depSiII',\
        'RVHb','EWHb','FWHb','FW14Hb','FW34Hb','gamma','depHb'),\
        dtype=('S16','float64','float64','float64',\
        'float64','float64','float64','float64',\
        'float64','float64','float64','float64',\
        'float64','float64','float64','float64',\
        'float64','float64','float64','float64',\
        'float64','float64','float64','float64','float64','float64','float64'))


    quit = ''
    for source in table:
        if quit == 'q': break

        id = source['ID']
        try: spt  = source['SpC']
        except: spt=''

        if id in [i.strip() for i in output['ID']]: continue

        skip = input("%s (%s) - Hit return to continue, type 's' to skip: " % (id,spt))

        next = 'n'
        while next == 'n':

            if skip == 's': break

            star = spec(id,SNR='best')

            snr_b = star.snrcalc(zone='B')
            snr_v = star.snrcalc(zone='V')
            snr_r = star.snrcalc(zone='R')

            #===================================================================
            print('\nAnalyzing Si III triplet...\n')
            star.plotspec(4500,4600)

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

            star.rv0 = RV0(RV0lines,star.spectrum,func=fun,ewcut=30,tol=RV0tol)
            star.waveflux(4500,6380); star.cosmic()
            star.plotspec(4500,4600,poslines='OB')

            input(); plt.close()

            fit = star.fitline(6347.11,width=wid,tol=tol,func=fun,plot=True)
            RVSi4,EWSi4,FWSi4,depSi4 = fit['RV_kms'],fit['EW'],fit['FWHM'],fit['depth']
            fit = star.fitline(4574.757,width=wid,tol=tol,func=fun,plot=True)
            RVSi3,EWSi3,FWSi3,depSi3 = fit['RV_kms'],fit['EW'],fit['FWHM'],fit['depth']
            fit = star.fitline(4567.84 ,width=wid,tol=tol,func=fun,plot=True)
            RVSi2,EWSi2,FWSi2,depSi2 = fit['RV_kms'],fit['EW'],fit['FWHM'],fit['depth']
            fit = star.fitline(4552.622,width=wid,tol=tol,func=fun,plot=True)
            RVSi1,EWSi1,FWSi1,depSi1 = fit['RV_kms'],fit['EW'],fit['FWHM'],fit['depth']

            if EWSi4 != None:
                if EWSi4 < ewcut: EWSi4 = FWSi4 = np.nan
            if EWSi3 != None:
                if EWSi3 < ewcut: EWSi3 = FWSi3 = np.nan
            if EWSi2 != None:
                if EWSi2 < ewcut: EWSi2 = FWSi2 = np.nan
            if EWSi1 != None:
                if EWSi1 < ewcut: EWSi1 = FWSi1 = np.nan

            next = input("Type 'n' to repeat, hit return to move to the Hb line. ")
            plt.close('all')

        next = 'n'
        while next == 'n':

            if skip == 's': break

            #=======================================================================
            print('\nAnalyzing H beta line...\n')

            star.waveflux(4801,4921); star.cosmic()
            star.plotspec(4821,4901,poslines='OB')

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

            fit = star.fitline(4861.325,width=wid,func=fun,iter=1,info=True,outfit=True,plot=True)
            RVHb,EWHb,FWHb,depHb = fit['RV_kms'],fit['EW'],fit['FWHM'],fit['depth']

            try:
                gamma = fit['gamma']

                wave = fit['wave']; flux_fit = fit['flux_fit']
                lowval = (max(flux_fit) + 3*min(flux_fit))/4
                uppval = (3*max(flux_fit) + min(flux_fit))/4

                FWs = []
                for val in lowval,uppval:
                    medpos = [np.where(flux_fit <= val)[0][value] for value in (0,-1)]
                    try: l_val = np.interp(val,[flux_fit[medpos[0]],flux_fit[medpos[0]-1]],
                                                   [wave[medpos[0]],wave[medpos[0]-1]])
                    except: l_val = wave[medpos[0]]
                    try: r_val = np.interp(val,[flux_fit[medpos[1]],flux_fit[medpos[1]+1]],
                                                  [wave[medpos[1]],wave[medpos[1]+1]])
                    except: r_val = wave[medpos[1]]
                    FWs.append(round(r_val-l_val,3))
            except:
                print('Line could not be fitted...')
                FWs = [np.nan,np.nan]; gamma = np.nan

            next = input("\nRepeat Hb / continue to the next star / save and exit ['n'/''/'q']: ")
            plt.close()
            if next == 'n': continue

            else:
                output.add_row(([id],[snr_b],[snr_v],[snr_r],\
                    [RVSi1],[EWSi1],[FWSi1],[depSi1],\
                    [RVSi2],[EWSi2],[FWSi2],[depSi2],\
                    [RVSi3],[EWSi3],[FWSi3],[depSi3],\
                    [RVSi4],[EWSi4],[FWSi4],[depSi4],\
                    [RVHb],[EWHb],[FWHb],[FWs[0]],[FWs[1]],[gamma],[depHb]))

                if next == 'q': quit = next

    output.write(maindir+'tables/'+output_table,format='fits',overwrite=True)

    return 'DONE'


def auto_measure(table='emulated_all.txt', output_table='Emul_RVEWFWs.fits',
    lines='emulated.lst', func='vrg_H', width=20, ewcut=10, tol=150, txt=False):

    '''
    Function to automatically calculate and store radial velocity, equivalent
    width and full width at half maximum of stars for input lines.

    Parameters
    ----------
    table : str
        Name of the input table contaning the list of stars to analyze.

    output_table : str
        Name of the output (new) table contaning the results.

    lines : str, list
        Enter the name of .txt/.lst file containing the lines and identifiers.
        (e.g. 4861.325 Hb)

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
    pp = PdfPages(maindir+'tmp_plots/RVEWFW.pdf')

    table = findtable(table)

    IDs = []; A1 = []; x0 = []; sig = []; sig1 = []; gamma = []; vsini = []; A2 = []; sig2 = []

    '''============================ Create table ============================'''
    lines,elements,_ = findlines(lines)

    columns = ['RV','EW','FW','dep','FW14','FW34']

    names = ['ID']+([j+i for i in elements for j in columns])
    dtypes =['S%i' % len(sorted(table['ID'],key=len)[-1]+' ')]
    dtypes += ['float64' for i in elements for j in columns]

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

        star = spec(id,SNR='best',txt=txt)

        fig = plt.figure()
        fig.suptitle(id,fontsize=9)

        row_data = []; nplot = 1
        for line,element in zip(lines,elements):

            print('\nAnalyzing %s...\n' % element)

            plt.subplot(nrows,ncols,nplot)

            #star.cosmic()

            fitting = star.fitline(line,width=width,func=func,iter=iter,output=True)
            row_data += fitting[2:6]

            try:
                wave,flux,flux_norm,flux_fit,popt = fitting[-1]

                IDs.append(id)
                if func in ['vrg_H','vrg_Z']:
                    A1.append(popt[0]);x0.append(popt[1]);sig1.append(popt[2]);gamma.append(popt[3])
                    vsini.append(popt[4]);A2.append(popt[5]);sig2.append(popt[6])
                elif func in ['vr_H','vr_Z']:
                    A1.append(popt[0]);x0.append(popt[1]);sig.append(popt[2]);gamma.append(popt[3])
                    vsini.append(popt[4])
                elif func in ['r']:
                    A1.append(popt[0]);x0.append(popt[1]);sig.append(popt[2]);vsini.append(popt[3])
                elif func in ['v']:
                    A1.append(popt[0]);x0.append(popt[1]);sig.append(popt[2]);gamma.append(popt[3])
                elif func in ['g','l']:
                    A1.append(popt[0]);x0.append(popt[1]);sig.append(popt[2])


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

            except:
                row_data += [np.nan]*2
                wave,flux,flux_norm,flux_fit,popt = [np.nan]*5

            plt.plot(wave,flux,'orange',lw=.5)
            plt.plot(wave,flux_norm,'b',lw=.5)
            plt.plot(wave,flux_fit,'g',lw=.5)

            plt.title(str(line)+' '+element,size=9,pad=4)
            plt.tick_params(direction='in',top='on')
            plt.ylim(ymax=1.03,ymin=0.4)
            plt.tight_layout()

            plt.close('all')
            nplot += 1

        pp.savefig(fig); plt.close(fig)

        output.add_row(([id]+row_data))

        bar.update(i)

    pp.close()

    output.write(maindir+'tables/'+output_table,format='fits',overwrite=True)

    bar.finish()

    # TO REMOVE TIL RETURN
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
    table = join(output,table,keys='ID',join_type='outer')
    table.write(maindir+'tables/'+output_table,format='fits',overwrite=True)
    # # # # # # # #

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
