from tkinter import font
from rv import *

import random
from copy import deepcopy
from matplotlib.backends.backend_pdf import PdfPages


def gen_ascii(id, txt=False, db_table=None, spt='auto', lwl=None, rwl=None, rv_corr=True, RV0tol=200, 
    export_rv=False, cosmetic=False, cosmic=False, lines_cosmic=None, degrade=None, show_plot=False):

    '''
    Function to remove cosmic rays in the spectra by different approaches.

    Parameters
    ----------
    id : str
        Name or filename of the star.

    txt : boolean, optional
        If True, it will search for two-column ascii files matching the input 'id'.
        See spec.txtwaveflux() for more details.

    db_table : str, optional
        Input table containing information of the spectral type for the star.
        The table must contain either 'SpT'/'SpT'/'SpC' columns.

    spt : str, optional
        Input spectral type of the star. If 'auto' (default), it takes it from either the
        db_table or the fits file. Otherwise enter an input type (e.g. B4Ia).

    lwl : float, optional
        Lower wavelength limit for the exported spectrum. Default is None.

    rwl : float, optional
        Upper wavelength limit for the exported spectrum. Default is None.

    rv_corr : boolean, optional
        If True, the output ascii spectrum will be corrected from radial velocity.
        Default is True.

    RV0tol : int, optional
        Enter the input radial velocity tolerance for the radial velocity correction.

    export_rv : boolean, optional
        If True and rv_corr also True, then the calculated radial velocity is exported
        as a table named 'gen_ascii_RVs.txt' with columns 'ID, RV0, eRV0'

    cosmetic : boolean, optional
        If True, the output ascii spectrum will be corrected from cosmetic defects.
        Default is False.

    cosmic : boolean, optional
        If True, the output ascii spectrum will be corrected from cosmic rays.
        Default is False.

    lines_cosmic : str/list, optional
        Enter the wavelenght(s) of the line(s) to show, either in a coma-separated
        string, or in a .txt/.lst file containing the lines to plot in order to review
        the cosmic rays removal.

    degrade : int/float, optional
        If True, the output ascii spectrum will degraded to the input resolution.
        Default is None.

    show_plot : boolean, optional
        If True, a plot showing the modifications to the original spectrum is shown.
        Default is False.

    Returns
    -------
    Nothing, but the ascii file for the input star is created.
    '''

    if db_table is not None:
        db_table = findtable(db_table)
        if not 'ID' in db_table.colnames:
            print('Column ID must be included in the table columns. Exiting...\n')
            return None

        row_id = db_table[db_table['ID']==id.split('_')[0]]

        if spt == 'auto':
            if 'SpT' in db_table.colnames:
                print('Initial spectral type taken from SpT column.')
                spt = row_id['SpT'][0]
            elif 'SpT' in db_table.colnames:
                print('Initial spectral type taken from SpT column.')
                spt = spc_code(row_id['SpT'])[0]
            elif 'SpC' in db_table.colnames:
                print('Initial spectral type taken from SpC column.')
                spt = spc_code(row_id['SpC'])[0]

    if spt == 'auto':
        if txt == False:
            print('Initial spectral type taken from FITS header.')
        else:
            print('Initial spectral type taken from Simbad query.')

        star = spec(id, SNR='bestHF', txt=txt)
        spt = spc_code(star.SpC)[0]

    elif isinstance(spt, str):
        spt = spc_code(spt)[0]

    if spt == '':
        spt = input('SpC keyword not set, please specify an SpC for the star: ')
        spt = spc_code(spt)[0]

    skip = input('%s - Hit return to continue, type "s" to skip: ' % id)
    if skip == 's':
        return None

    finish = 'n'
    while finish == 'n':

        star = spec(id, SNR='bestHF', txt=txt, cut_edges=True)

        if show_plot == True:
            fig,axs = plt.subplots(3, 1, figsize=(14,8))
            fig.suptitle(star.filename)
            axs[0].tick_params(direction='in', top='on')
            axs[0].set_xlim(3950, 6850)
            axs[0].set_ylim(0.5, 1.1)
            axs[0].plot(star.wave_0, star.flux_0, c='orange', lw=.2, label='RV corrected')

        # Correct the spectrum form cosmetic defects:
        if cosmetic == True:

            if '_F_' in star.filename:
                i = 0
                std = np.std(star.flux[(star.wave > 4968) & (star.wave < 4985)])
                for win in [
                [4089.9,4090.3], [4505.,4508.], [4693.0,4695.7],
                [4794.3,4796.5], [4900.7,4902.2], [6636.,6638.],
                ]:
                    mask = (star.wave > win[0]) & (star.wave < win[1])
                    gap = star.flux[mask]
                    if np.std(gap) < 5*std:
                        continue # To skip if the artifact is not there
                    length = len(gap)
                    cont = np.mean(np.concatenate((gap[:5], gap[-5:])))
                    star.flux[mask] = [random.uniform(cont-2*std, cont+2*std) for j in range(length)]
                    axs[0].plot(star.wave[mask], star.flux[mask],
                        c='r', lw=.2, label='Cosmetic fix' if i==0 else '')
                    i += 1

        # Correct the spectrum form radial velocity:
        if rv_corr == True:

            if spt <= 2.:
                spt_list = 'rv_Os.lst'
            elif spt > 2. and spt < 2.7:
                spt_list = 'rv_Bs.lst'
            elif spt >=2.7:
                spt_list = 'rv_As.lst'

            next_rv0 = 'n'; fun = 'g'; wid = 15; tmp_wave = star.wave
            while next_rv0 == 'n':

                print('Current input for RV0 correction is list({})'.format(spt_list)
                    +' / function({})'.format(fun)+' / width({})'.format(wid))
                change = input('To change type list/funcion/width, otherwise hit return: ')

                plt.close('all')

                if change == 'list':
                    spt_list = '-'
                    while spt_list not in ['rv_Os.lst','rv_Bs.lst','rv_As.lst']:
                        spt_list = input('Choose list of lines between O/B/A: ')
                        if spt_list == '' or spt_list in ['B','b']: SpT = 'B'; spt_list = 'rv_Bs.lst'
                        elif spt_list in ['A','a']: SpT = 'A'; spt_list = 'rv_As.lst'
                        elif spt_list in ['O','o']: SpT = 'O'; spt_list = 'rv_Os.lst'
                elif change in ['fun','function','func']:
                    fun = '-'
                    while fun not in ['g','l','v','r','vr']:
                        fun = input('Choose function to fit between g/l/v/r/vr: ')
                        if fun == '': fun = 'g'
                elif change in ['width','wid']:
                    wid = '-'
                    while type(wid) is not float:
                        wid = input('Choose the initial width in angstroms: ')
                        if wid == '': wid = 15.
                        else: wid = float(wid)

                star.rv0, erv0 = RV0(spt_list, star.filename, txt=txt, ewcut=30, width=wid, tol=RV0tol, func=fun)
                star.wave = tmp_wave*(1 - 1000*star.rv0/cte.c)

                star.plotspec(4821,4901, lines='35-10K')
                plt.figure()
                if spt <= 2.5:
                    star.plotspec(4530, 4590, lines='35-10K')
                else:
                    star.plotspec(6361.37, 6381.37, lines='35-10K')

                next_rv0 = input("Type 'n' to repeat, hit return to continue. ")
                if next_rv0 not in ['n','']:
                    next_rv0 = 'n'

            plt.close('all')

            if export_rv == True:
                rv_table = open(maindir+'tables/gen_ascii_RVs.txt', 'a+')
                rv_table.write('{}, {}, {:.3f}, {:.3f}\n'.format(star.id_star, star.filename, star.rv0, erv0))
                rv_table.close()

        elif rv_corr not in [True,False]:
            try:
                star.rv0 = float(rv_corr)
                star.wave = star.wave*(1 - 1000*star.rv0/cte.c)
                print('Correcting spectrum using the input RV correction.')
            except:
                print('Spectrum could not be corrected using the input RV correction.')

        # Correct the spectrum form cosmic rays:
        if cosmic == True:

            next_cosm = 'n'; dmin = 0.05; zs_cut = 3; blue_cut = '4000'
            while next_cosm == 'n':

                print('Current input for cosmic rays correction is zs_cut({})'.format(zs_cut)
                    +' / dmin({})'.format(dmin))
                change = input('To change type cut/dmin, otherwise hit return: ')

                if change in ['cut','zs_cut']:
                    zs_cut = '-'
                    while type(zs_cut) is not float:
                        zs_cut = input('Choose a new zs_cut value: ')
                        if zs_cut != '': zs_cut = float(zs_cut)
                elif change == 'dmin':
                    dmin = '-'
                    while type(dmin) is not float:
                        dmin = input('Choose a new dmin value: ')
                        if dmin != '': dmin = float(dmin)

                tmp_star = deepcopy(star)
                tmp_star.cosmic(method='zscore', dmin=dmin, zs_cut=zs_cut)

                # To prevent the noisier blue part of the spectrum to be taken for cosmic removal:
                try:
                    blue_cut = float(input('Set initial wavelength from where to apply the removal (now is %s): ' % str(blue_cut)))
                except:
                    print('Only int/float are accepted. Choosing %s... ' % str(blue_cut))
                    blue_cut = float(blue_cut)

                tmp_star.flux = np.concatenate([star.flux[star.wave<blue_cut],tmp_star.flux[star.wave>=blue_cut]])

                fig_cosm,ax_cosm = plt.subplots(figsize=(14,4))
                ax_cosm.plot(star.wave, star.flux, c='orange', lw=.7, label='RV corrected')
                ax_cosm.plot(star.wave, tmp_star.flux, c='b', lw=.5, label='Cosmic corrected')

                if ax_cosm.get_ylim()[1] > 4:
                    ax_cosm.set_ylim(top=4)
                if ax_cosm.get_ylim()[0] <0:
                    ax_cosm.set_ylim(bottom=0)

                ax_cosm.legend()
                fig_cosm.tight_layout()
                fig_cosm.show()

                if lines_cosmic is not None:

                    lines,elems,_ = findlines(lines_cosmic)
                    nrows = int(np.ceil(np.sqrt(len(lines))))
                    ncols = round(len(lines)/np.ceil(np.sqrt(len(lines)))+0.4)

                    fig_lines,ax_lines = plt.subplots(nrows, ncols, figsize=(13,8.5))
                    fig_lines.suptitle(star.filename, y=0.97, fontsize=8)

                    ax_lines = ax_lines.flatten()
                    for ax_i,line,elem in zip(ax_lines,lines,elems):

                        if elem in ['Hdelta','Hgamma','Hbeta','Halpha']:
                            width = 60
                        else:
                            width = 15

                        mask = (star.wave >= line-width/2.) & (star.wave <= line+width/2.)
                        ax_i.plot(star.wave[mask], star.flux[mask], c='orange', lw=.7, label='RV corrected')
                        ax_i.plot(star.wave[mask], tmp_star.flux[mask], c='b', lw=.5, label='Cosmic corrected')
                        ax_i.set_title(elem, fontsize=6, pad=0.55)
                        ax_i.tick_params(direction='in', top='on', labelsize=5)
                        ax_i.set_yticks([])
                        if ax_i.get_ylim()[0] > 0.9:
                            ax_i.set_ylim(bottom=0.9)

                    [fig_lines.delaxes(ax_lines[i]) for i in np.arange(len(lines), len(ax_lines), 1)]

                    fig_lines.tight_layout()
                    fig_lines.subplots_adjust(wspace=.05, hspace=0.25)

                    fig_lines.show()

                next_cosm = input("Type 'n' to repeat, hit return to continue. ")
                if next_cosm not in ['n','']:
                    next_cosm = 'n'

                if next_cosm == '':
                    star.flux = tmp_star.flux

                plt.close('all')

            if show_plot == True:
                axs[1].tick_params(direction='in', top='on')
                axs[1].set_xlim(3950, 6850)
                axs[1].set_ylim(0.5, 1.1)
                axs[1].plot(star.wave, star.flux, c='g', lw=.2, label='Cosmic corrected')

        # Degrade the spectra to a different resolution
        if degrade is not None:
            while type(degrade) is not float:
                try:
                    degrade = float(degrade)
                except:
                    degrade = float(input('Choose a valid degrading resolution: '))

            star.degrade(resol=degrade)

        # Show the final plot:
        if show_plot == True:
            axs[2].tick_params(direction='in', top='on')
            axs[2].set_xlim(3950, 6850)
            axs[2].set_ylim(0.5, 1.1)
            axs[2].plot(star.wave, star.flux, c='b', lw=.2, label='Final')

            fig.tight_layout()
            fig.legend(ncol=3)
            fig.show()

        finish = input('Type "n" to repeat, hit return to move to the next star. ')
        if finish not in ['n','']:
            finish = 'n'

    plt.close(fig)
    plt.close('all')
    plt.close()

    # Cut the spectra to the selected limits:
    if lwl is not None and rwl is not None:
        mask = (star.wave >= lwl) & (star.wave <= rwl)
        star.wave = star.wave[mask]
        star.flux = star.flux[mask]

    star.export(tail='_RV', extension='.ascii')

    return None


def gen_ascii_ML(input_table='OBAs_ML_raw.fits', not_do=None, cosmic_manual=False, txt=False):

    '''
    IN DEVELOPMENT

    Similar function as 'gen_ascii' but only to quickly generate ascii spectra to be used
    in Matchine Learning projects, as the spectra is degraded and resampled.

    Parameters
    ----------
    input_table : str, optional
        Input table from where to take the stars ID and spectral type.
        Temporal default is 'OBAs_ML_raw.fits'

    not_do : list, optional
        List of IDs to skip.

    cosmic_manual : boolean, optional
        If True, the user manually supervises the cosmic ray removal step.

    Returns
    -------
    Nothing, output ascii files with the spectra are generated.
    '''

    if type(table) is type(Table()): pass # In case the input table is already a table
    else: table = findtable(table) # file where star names and quality flags are

    output = open(maindir + 'tmp/results_ML.txt', 'a')
    pp = PdfPages(maindir + 'tmp_plots/ML_results.pdf')

    for row in table[:3]:

        # Skip all sources in the 'not_do' input list :
        if not_do is not None and row['ID'] in not_do: continue

        if row['SNR_best'] < 60: continue
        else: print('Analysing %s' % row['ID'])

        # Determines best list for RV calculation and line for sanity check based on SpT
        if row['SpT'] < 2:
            rv_list = 'rv_Os.lst'
            line = 5411.52 # 5592.252

        elif 2 <= row['SpT'] < 2.5:
            rv_list = 'rv_Bs.lst'
            line = 4552.622

        elif 2.5 <= row['SpT'] < 2.9:
            rv_list = 'rv_Bs.lst'
            line = 6371.37

        elif row['SpT'] >= 2.9:
            rv_list = 'rv_As.lst'
            line = 4233.129

        # Determines the RV with a default fitting function and width
        fun = 'g'; wid = 15
        best_star = spec(row['ID'], SNR='bestHF', txt=txt)
        best_star.rv0, erv0 = RV0(rv_list, best_star.filename, txt=txt, ewcut=50, func=fun, width=wid, tol=150)
        best_star.waveflux(3950,6850, cut_edges=True) # Applies the rv0 correction

        # If the line is >.1A from where should be, you determine new best function, width and line
        RV_A = abs(best_star.fitline(line, func=fun, width=wid, tol=20)['RV_A'])
        if np.isnan(RV_A) or RV_A > 0.1: # 5km/s at 5400A
            print('Fitted line offset is %.3f (tolerance is 0.1A)' % RV_A)

            next_rv0 = 'n'; tmp_wave = best_star.wave
            while next_rv0 == 'n':

                print('Current input for RV0 correction is function({})'.format(fun)+' / width({})'.format(wid))
                change = input('To change type funcion/width, otherwise hit return: ')

                plt.close('all')

                if change in ['fun','function','func']:
                    fun = '-'
                    while fun not in ['g','l','v','r','vr']:
                        fun = input('Choose function to fit between g/l/v/r/vr: ')
                        if fun == '': fun = 'g'
                elif change in ['width','wid']:
                    wid = '-'
                    while type(wid) is not float:
                        wid = input('Choose the initial width in angstroms: ')
                        if wid == '': wid = 15.
                        else: wid = float(wid)

                best_star.rv0, erv0 = RV0(rv_list, best_star.filename, ewcut=30, func=fun, width=wid, tol=150)
                best_star.wave = tmp_wave*(1 - 1000*best_star.rv0/cte.c) # Applies the rv0 correction

                best_star.plotspec(4821,4901, lines='35-10K')
                plt.figure()
                if row['SpT'] <= 2.5:
                    best_star.plotspec(4530, 4590, lines='35-10K')
                else:
                    best_star.plotspec(6361.37, 6381.37, lines='35-10K')

                next_rv0 = input("Type 'n' to repeat, hit return to continue with last chosen parameters. ")

        # For all available spectra withing a limit, create the final output ascii
        goodspec = findstar(row['ID'], SNR=60)
        if len(goodspec) > 10:
            goodspec = random.sample(goodspec, 10)

        for j,n in zip(goodspec, range(len(goodspec))):
            star = spec(j.split(os.sep)[-1], txt=txt)
            star.waveflux(3950,6850, cut_edges=True)

            # Create the plot:
            fig,axs = plt.subplots(3, 1, figsize=(14,8))
            fig.suptitle(star.filename)
            axs[0].tick_params(direction='in', top='on')
            axs[0].set_xlim(3950, 6850)
            axs[0].set_ylim(0.5, 1.1)
            axs[1].tick_params(direction='in', top='on')
            axs[1].set_xlim(3950, 6850)
            axs[1].set_ylim(0.5, 1.1)

            axs[0].plot(star.wave, star.flux, c='orange', lw=.2, label='Original')

            # Remove artifacts in FEROS spectra
            if '_F_' in star.filename:
                i = 0
                std = np.std(star.flux[(star.wave > 4968)&(star.wave < 4985)])
                for win in [
                [4089.9,4090.3], [4505.,4508.], [4693.0,4695.7],
                [4794.3,4796.5], [4900.7,4902.2], [6636.,6638.],
                ]:
                    mask = (star.wave > win[0]) & (star.wave < win[1])
                    gap = star.flux[mask]
                    if np.std(gap) < 5*std: continue # To skip if the artifact is not there
                    length = len(gap)
                    cont = np.mean(np.concatenate((gap[:5],gap[-5:])))
                    star.flux[mask] = [random.uniform(cont-2*std,cont+2*std) for j in range(length)]
                    axs[0].plot(star.wave[mask], star.flux[mask],
                        c='r', lw=.2, label='FEROS fixed issues' if i==0 else "")
                    i += 1

            star.rv0, erv0 = RV0(rv_list, star.filename, ewcut=50, width=wid, tol=150, func=fun)
            star.wave = star.wave*(1 - 1000*star.rv0/cte.c) # Applies the rv0 correction

            axs[1].plot(star.wave, star.flux, c='orange', lw=.2, label='Original*')

            # Correct the spectrum form cosmic rays:
            if cosmic_manual == True:
                next_cosm = 'n'; dmin = 0.05; zs_cut = 4
                while next_cosm == 'n':

                    print('Current input for cosmic rays correction is zs_cut({})'.format(zs_cut)
                        +' / dmin({})'.format(dmin))
                    change = input('To change type cut/dmin, otherwise hit return: ')

                    if change in ['cut','zs_cut']:
                        zs_cut = '-'
                        while type(zs_cut) is not float:
                            zs_cut = input('Choose a new zs_cut value: ')
                            if zs_cut != '': zs_cut = float(zs_cut)
                    elif change == 'dmin':
                        dmin = '-'
                        while type(dmin) is not float:
                            dmin = input('Choose a new dmin value: ')
                            if dmin != '': dmin = float(dmin)

                    tmp_star = deepcopy(star)
                    tmp_star.cosmic(method='zscore', dmin=dmin, zs_cut=zs_cut, protect_em_lines=True)

                    fig_cosm,ax_cosm = plt.subplots(figsize=(14,4))
                    ax_cosm.plot(star.wave, star.flux, c='orange', lw=1, label='RV corrected')
                    ax_cosm.plot(star.wave, tmp_star.flux, c='g', lw=.5, label='Cosmic corrected')
                    fig_cosm.tight_layout()
                    fig_cosm.show()

                    next_cosm = input("Type 'n' to repeat, hit return to continue. ")

            else:
                tmp_star = deepcopy(star)
                tmp_star.cosmic(method='zscore', zs_cut=6, dmin=0.05, protect_em_lines=True)

            i = 0
            for wl_em in [3967.79,4958.911,5006.843,6300.304,6548.04,6583.46,6716.44,6730.82]:
                mask = (star.wave > wl_em-0.8) & (star.wave < wl_em+0.8)
                axs[1].plot(star.wave[mask], tmp_star.flux[mask], c='g', lw=.2, label='Em. lines' if i==0 else "", zorder=10)
                i += 1

            if tmp_star.flux.min() < 0:
                tmp_star.flux = np.where(tmp_star.flux < 0, 0.0, tmp_star.flux)

            axs[1].plot(star.wave, tmp_star.flux, c='b', lw=.2, label='Final')

            star.flux = tmp_star.flux

            # Degrade the spectra to a different resolution
            star.degrade(resol=5000)

            mask = (star.wave > line-10) & (star.wave < line+10)
            axs[2].plot(star.wave[mask], star.flux[mask], c='b', lw=.2, label='Final spectra')
            axs[2].plot([line,line], [min(star.flux[mask]),max(star.flux[mask])], c='k', label='At. line position')
            medval = (max(star.flux[mask]) + min(star.flux[mask]))/2
            medpos = [np.where(star.flux[mask] <= medval)[0][value] for value in (0,-1)]
            center = round((star.wave[mask][medpos[1]]+star.wave[mask][medpos[0]])/2,3)
            axs[2].plot([center,center], [min(star.flux[mask]),max(star.flux[mask])], c='g', label='Aprox. center')
            axs[2].set_title('Difference in angstroms is: %.5f' % abs(center - line))
            axs[2].set_ylim()

            # Resample the spectra
            star.resamp(10*0.02564975, 3950, 6850)

            star.export(tail='_RV_ML_%i' % n, extension='.ascii')
            if max(star.flux) > 1.5:
                output.write(star.filename + ' | Max flux: %.3f\n' % max(star.flux))
            if abs(center - line) > 0.2:
                output.write(star.filename + ' | RV difference (A): %.3f\n' % (center - line))

            axs[0].legend(); axs[1].legend(); axs[2].legend()

            fig.tight_layout()
            pp.savefig(fig)
            plt.close('all')

    output.close()
    pp.close()
    np.savetxt(maindir + 'tmp/wavelenghtML.ascii', np.c_[star.wave], fmt=('%.4f'))


def remove_wave(path=maindir+'tmp/', only_list='to_correct.txt'):

    '''
    To remove the wavelenght columns from the ascii files containing the spectrum.
    '''

    if only_list != None:
        table = findtable(only_list, delimiter=' ')

    for file in os.listdir(path):
        if file.endswith('.ascii'):

            if only_list != None and not file in table['File']:
                continue

            data = Table.read(path + file, format='ascii', delimiter=' ')
            try:
                if max(data['col2']) > 1.5:
                    plt.plot(data['col1'], data['col2'], lw=.5)
                    plt.plot([3950,6850], [2,2], c='k', lw=.5)
                    plt.show(block=False)
                    inp = input('Wanna print data for %s it? [y/ ]: ' % file)
                    if inp == 'y':
                        print(file,'>1.5',max(data['col2']))
                    plt.close()
                if min(data['col2']) < 0: print(file,'<0',min(data['col2']))
                data.remove_columns(['col1'])
                data.write(maindir + 'tmp/new/' + file, format='ascii.no_header')
            except:
                print(file)

        else: continue


#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# Update tables with new data | / IMPLEMENT IN DB AT SOME POINT!!
#table1,table2 = findtable('OBAs_ML_raw.fits'),findtable('MAUI_results.fits')
#columns_to_update = table2.colnames[1:-9]+[table2.colnames[-1]]
#for row,i in zip(table1,range(len(table1))):
#    updated = False
#    if row['ID'] not in table2['ID']: continue
#    if row['Teff'] == table2[table2['ID'] == row['ID']]['Teff'].data[0]: continue
#    for col_name in columns_to_update:
#        table1[i][col_name] = table2[table2['ID'] == row['ID']][col_name].data[0]
#        if updated == False: print(row['ID']); updated = True
#table1.write(maindir+'tables/table1_updated.fits',format='fits',overwrite=True)

# This one is to empty bad data
#table1 = findtable('OBAs_ML_raw.fits')
#columns_to_update = [i for i in table1.columns[36:-1]]
#for row,i in zip(table1,range(len(table1))):
#    for col_name in columns_to_update:
#        if row[col_name] in ['d','<','>']:
#            for j in columns_to_update[columns_to_update.index(col_name)+1:columns_to_update.index(col_name)+4]:
#                table1[i][j] = np.nan
#table1.write(maindir+'tables/OBAs_ML_ver1.fits',format='fits',overwrite=True)

#table1 = findtable('OBAs_ML_ver1b.fits')
#for row,i in zip(table1,range(len(table1))):
#    if row['QSiIII'] < 3:
#        table1[i]['EWSiIII1'] = table1[i]['FWSiIII1'] = table1[i]['depSiIII1'] = np.nan
#        table1[i]['EWSiIII2'] = table1[i]['FWSiIII2'] = table1[i]['depSiIII2'] = np.nan
#        table1[i]['EWSiIII3'] = table1[i]['FWSiIII3'] = table1[i]['depSiIII3'] = np.nan
#    if row['QSiII'] < 3:
#        table1[i]['EWSiII'] = table1[i]['FWSiII'] = table1[i]['depSiII'] = np.nan
#    if row['QHb'] < 3:
#        table1[i]['EWHb'] = table1[i]['FWHb'] = table1[i]['FW14Hb'] = table1[i]['FW34Hb'] = table1[i]['depHb'] = table1[i]['gamma'] = np.nan
#table1.write(maindir+'tables/OBAs_ML_ver1.fits',format='fits',overwrite=True)

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# Replace values in table
#table = findtable('OBAs_ML_raw.fits')
#for row,i in zip(table,range(len(table))):
#    for column,j in zip(row,range(len(row))):
#        #if column == 'N': table[i][j] = '='
#        if str(column) == str(1e+20): table[i][j] = np.nan
#table.write(maindir+'tables/OBAs_ML_raw_.fits',format='fits',overwrite=True)
