from spec import *
from RV import *


def RVEWFW(table='IACOB_O9BAs_SNR20.fits', output_table='O9BAs_RVEWFWs.fits',
           lines='rv_Bs.lst', RV0tol=150, ewcut=10, tol=100, redo='n'):
    '''
    Function to iteratively calculate and store EWs and FWHMs of stars using high-
    resolution spectra, adjusted to find BSGs (i.e. using SiIII triplet + Hb line)


    Parameters
    ----------
    table : str
        Name of the input table contaning the list of stars to analyze.

    output_table : str
        Name of the output (new) table contaning the results.

    lines : str, list
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
    '''

    table = findtable(table)
    #table = findtable('Std_09-B8.fits')

    '''===================== Create table if not exist ======================'''
    try: output = findtable(output_table)
    except:
        format1 = str(table['Name'].info.dtype)[1:]
        format2 = str(table['SpC'].info.dtype)[1:]

        output = Table(names=(\
        'Name','mySpC','SpC','SpT_code','LC_code','SNR_B',\
        'RVSiIII1','EWSiIII1','FWSiIII1','depSiIII1',\
        'RVSiIII2','EWSiIII2','FWSiIII2','depSiIII2',\
        'RVSiIII3','EWSiIII3','FWSiIII3','depSiIII3',\
        'RVSiII','EWSiII','FWSiII','depSiII',\
        'RVHb','EWHb','FWHb','FW14Hb','FW34Hb','gamma','depHb',\
        'Comments'),\
        dtype=('S16','S12','S16','float64','float64','float64',\
        'float64','float64','float64','float64',\
        'float64','float64','float64','float64',\
        'float64','float64','float64','float64',\
        'float64','float64','float64','float64',\
        'float64','float64','float64','float64','float64','float64','float64',\
        'S80'))

    if redo != 'n':
        try: redo = redo.split(',')
        except: print('Bad input for "redo" parameter. Exitting...'); return None
        for name in redo:
            output = output[output['Name'] != name]

    quit = ''
    for row in table:

        if quit == 'quit': break

        name = row['Name'].strip()
        spt  = row['SpC'].strip()

        if name in [i.strip() for i in output['Name']]: continue

        skip = input("%s (%s) - Hit return to continue, type 's' to skip: " % (name,spt))

        next = 'n'
        while next == 'n':

            if skip == 's': break

            star = spec(name,SNR='best')

            #=======================================================================
            snr_b = int(round(star.snrcalc(zone='B')))

            #=======================================================================
            print('\nAnalyzing Si III triplet...\n')
            star.plotspec(4500,4600)#4650

            fun = '-'
            while fun not in ['g','l','v','r','vr','dwarf']:
                fun = input('Choose function to fit between g/l/v/r/vr/dwarf (default is g): ')
                if fun == '': fun = 'g'
            wid = '-'
            while type(wid) is not float:
                wid = input('Choose the initial width in angstroms (default is 15): ')
                if wid == '': wid = 15.
                else:
                    try: wid = float(wid)
                    except: wid = '-'

            plt.close()

            star.offset = RV0(lines,star.spectrum,func=fun,ewcut=30,tol=RV0tol)
            star.waveflux() # Applies the offset
            star.cosmic(sigclip=0.002)

            star.plotspec(4500,4600,poslines='OB')

            input(); plt.close()

            RVSi4,EWSi4,FWSi4,depSi4 = \
            star.fitline(6347.11,width=wid,tol=tol,func=fun,plot='y')[2:6]
            RVSi3,EWSi3,FWSi3,depSi3 = \
            star.fitline(4574.757,width=wid,tol=tol,func=fun,plot='y')[2:6]
            RVSi2,EWSi2,FWSi2,depSi2 = \
            star.fitline(4567.84 ,width=wid,tol=tol,func=fun,plot='y')[2:6]
            RVSi1,EWSi1,FWSi1,depSi1 = \
            star.fitline(4552.622,width=wid,tol=tol,func=fun,plot='y')[2:6]

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

            star.waveflux(); star.cosmic()
            star.plotspec(4821,4901,poslines='OB')

            fun = '-'; iter = 3
            while fun not in ['vr','dwarf']:
                fun = input('Choose function to fit between vr/dwarf (default is dwarf): ')
                if fun == '': fun = 'dwarf'; iter = 1

            wid = '-'
            while type(wid) is not float:
                wid = input('Choose the initial width in angstroms (default is 50): ')
                if wid == '': wid = 50.
                else:
                    try: wid = float(wid)
                    except: wid = '-'

            plt.close()

            RVHb,EWHb,FWHb,depHb =\
            star.fitline(4861.325,width=wid,func=fun,iter=iter,output='y',plot='y')[2:6]

            try:
                wave,flux_fit,popt =\
                star.fitline(4861.325,width=wid,func=fun,iter=iter,output='y',outfit=True)

                gamma = round(popt[3],2)

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
                FWs = popt = [np.nan,np.nan,np.nan,np.nan,np.nan]; gamma = np.nan

            next = input("\nRepeat Hb / continue to the next star / save and exit ['n'/''/'q']: ")
            plt.close()
            if next == 'n': continue

            else:
                mySpC = input('Enter visual SpC (default "OB"): ')
                if mySpC == '': mySpC = 'OB'
                comments = input('Add any comment here, otherwise hit return: ')

                output.add_row(([name],
                    [mySpC],[row['SpC']],[row['SpT_code']],[row['LC_code']],[snr_b],\
                    [RVSi1],[EWSi1],[FWSi1],[depSi1],\
                    [RVSi2],[EWSi2],[FWSi2],[depSi2],\
                    [RVSi3],[EWSi3],[FWSi3],[depSi3],\
                    [RVSi4],[EWSi4],[FWSi4],[depSi4],\
                    [RVHb],[EWHb],[FWHb],[FWs[0]],[FWs[1]],[gamma],[depHb],\
                    [comments]))

                if next == 'q': quit = next

    hdu = fits.BinTableHDU(data=output.filled(np.nan))
    hdu.writeto(maindir+'tables/'+output_table,overwrite=True)

    return('DONE')


# %% ===========================================================================
#''' Scrip to show the stars for which new spectra with higher SNR is available '''
#table_old = findtable('IACOB_O9-B7_SNR20.fits')
#table_new = findtable('IACOB_O9BAs_SNR20.fits')
#
#for row in table_new:
#    snr_old = table_old[table_old['Name']==row['Name']]['SNR_best']
#    if row['SNR_best'] > snr_old:
#        print(str(row['Name']).strip(),row['SNR_best'],float(snr_old))
