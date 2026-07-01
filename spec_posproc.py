from rv import *

import random
from copy import deepcopy
from matplotlib.backends.backend_pdf import PdfPages


def gen_ascii(id, orig='IACOB', rv_corr=True, rv_method='fitting', rv_tol=200, export_rv=False, spt_table=None, spt='auto',
              folder_ccf=None, lwl=None, rwl=None, cosmetic=False, cosmic=False, lines_cosmic=None, degrade=None,
              show_plot=False):

    '''
    Function to generate ascii spectra from the FITS files of the IACOB database, with
    the possibility to apply different corrections to the spectra.

    Parameters
    ----------
    id : str
        Name or filename of the star. If the input is a name, the function will search
        for the best available spectrum in the IACOB database.

    orig : str, optional
        See spec() function for more information. Default is 'IACOB'.

    rv_corr : boolean, optional
        If True, the output ascii spectrum will be corrected from radial velocity.
        Default is True.

    rv_method : str, optional
        If rv_corr is True, choose the method to calculate the radial velocity correction:
        - 'fitting' : Uses a pre-defined list of lines depending on the spectral type.
                    --> See spt_table and spt input parameters too.
        - 'CCF' : Uses a cross-correlation function with a template spectrum.
                    --> See folder_ccf input parameter too.
                    NOTE: lwl and rwl are parsed to the RV0_cc() function, and thus
                    are used to limit the wavelength range for the CCF calculation.

    rv_tol: int, optional
        If rv_corr is True, and rv_method is 'fitting':
        Enter the input radial velocity tolerance use for discarding bad-fitted lines.
        Default is 200 km/s. See also "tol" input parameter in rv.RV0() function.

    export_rv : boolean, optional
        If True and rv_corr also True, then the calculated radial velocity is exported
        as a table named 'gen_ascii_RVs.txt' with columns 'ID, RV0, eRV0'

    spt : str, optional
        Input spectral type of the star. If 'auto' (default), it takes it from either the
        spt_table or the fits file. Otherwise enter an input type (e.g. B4Ia).

    spt_table : str, optional
        Input table containing information of the spectral type for the star.
        The table must contain either 'SpT'/'SpT'/'SpC' columns.

    folder_ccf : str, optional
        Folder containing the CCF templates for the radial velocity correction.
        The function will search for a template with the same ID as the input star.

    lwl : float, optional
        Lower wavelength limit for the exported spectrum. Default is None.

    rwl : float, optional
        Upper wavelength limit for the exported spectrum. Default is None.

    cosmetic : boolean, optional
        If True, the output ascii spectrum will be corrected from cosmetic defects.
        Default is False.

    cosmic : boolean, optional
        If True, the output ascii spectrum will be corrected from cosmic rays.
        Default is False.

    lines_cosmic : str/list, optional
        Enter the wavelenght(s) of the line(s) to show, either in a coma-separated
        string, or in a .txt/.txt file containing the lines to plot in order to review
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

    msg.info('Generating ascii spectrum for %s' % id)

    if rv_corr == True and rv_method == 'fitting' and spt_table is not None:
        spt_table = findtable(spt_table)
        if not 'ID' in spt_table.colnames:
            msg.error('Column ID must be included in the table columns. Exiting...\n')
            return None

        row_id = spt_table[spt_table['ID']==id.split('_')[0]]

        if spt == 'auto':
            if 'SpT' in spt_table.colnames:
                msg.info('Initial spectral type taken from SpT column.')
                spt = row_id['SpT'][0]
            elif 'SpC' in spt_table.colnames:
                msg.info('Initial spectral type taken from SpC column.')
                spt = spc_code(row_id['SpC'][0])[0]
            else:
                star = spec(id, snr='bestHF', orig=orig)
                if orig == 'IACOB':
                    msg.info('Initial spectral type taken from FITS header.')
                else:
                    msg.info('Initial spectral type taken from Simbad query.')
                    star.get_spc()
                spt = spc_code(star.SpC)[0]
        elif isinstance(spt, str) and spt not in ['auto','']:
            spt = spc_code(spt)[0]
        else:
            msg.warn('Input spectral is not set, or not valid. Exiting...\n')
            spt = input('SpC keyword not set, please specify an SpC for the star: ')
            spt = spc_code(spt)[0]

        msg.info('Initial spectral type is set to: %s' % spt)

    elif rv_corr == True and rv_method.lower() == 'ccf':
        if folder_ccf is None:
            msg.error('Input folder for CCF templates is not set. Exiting...\n')
            return None
        elif folder_ccf is not None and not os.path.exists(folder_ccf):
            msg.error('Input folder for CCF templates does not exist. Exiting...\n')
            return None

        # find a template file for the CCF with the ID in the filename:
        ccf_file = [ccf for ccf in os.listdir(folder_ccf) if id.split('_')[0] in ccf]
        while len(ccf_file) == 0 or len(ccf_file) > 1:
            for ccf in ccf_file:
                print('  - %s' % ccf)
            ccf_input = input('Please type a file to use as template for the CCF: ')
            ccf_file = [ccf for ccf in os.listdir(folder_ccf) if ccf_input in ccf]

        msg.info('Using %s as template for the CCF.' % ccf_file[0])

    skip = input('%s - Hit return to continue, type "s" to skip: ' % id)
    if skip == 's':
        return None

    finish = 'n'
    while finish == 'n':

        if type(rv_corr) is float:
            rv0 = rv_corr
        else:
            rv0 = 0

        star = spec(id, snr='bestHF', orig=orig, cut_edges=True, rv0=rv0)

        if show_plot == True:
            fig,axs = plt.subplots(3, 1, figsize=(14,8))
            fig.suptitle(star.filename)
            axs[0].tick_params(direction='in', top='on')
            axs[0].set_xlim(3950, 6850)
            axs[0].set_ylim(0.5, 1.1)
            label = 'RV corrected' if rv_corr == True else 'Original'
            axs[0].plot(star.wave_0, star.flux_0, c='orange', lw=.2, label=label)

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

        # Correct the spectrum form radial velocity using the fitting method:
        if rv_corr == True and rv_method.lower() == 'fitting':
            if spt <= 2.:
                spt_list = 'rv_Os.txt'
            elif spt > 2. and spt < 2.7:
                spt_list = 'rv_Bs.txt'
            elif spt >=2.7:
                spt_list = 'rv_As.txt'

            next_rv0 = 'n'; fun = 'g'; wid = 15; tmp_wave = star.wave
            while next_rv0 == 'n':

                print('Current input for RV0 correction is list({})'.format(spt_list)
                    +' / function({})'.format(fun)+' / width({})'.format(wid))
                change = input('To change type list/funcion/width, otherwise hit return: ')

                plt.close('all')

                if change == 'list':
                    spt_list = '-'
                    while spt_list not in ['rv_Os.txt','rv_Bs.txt','rv_As.txt']:
                        spt_list = input('Choose list of lines between O/B/A: ')
                        if spt_list == '' or spt_list in ['B','b']: SpT = 'B'; spt_list = 'rv_Bs.txt'
                        elif spt_list in ['A','a']: SpT = 'A'; spt_list = 'rv_As.txt'
                        elif spt_list in ['O','o']: SpT = 'O'; spt_list = 'rv_Os.txt'
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

                star.rv0, erv0 = RV0(spt_list, star.filename, orig=orig, ewcut=30, width=wid, tol=rv_tol, func=fun)
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

        if rv_corr == True and rv_method.lower() == 'ccf':
            star.rv0, erv0, _ , _ = RV0_cc(star.filename, ccf_file[0], orig1=orig, orig2='ascii',
                                    outside_dir=folder_ccf, method='windows', lwl=lwl, rwl=rwl)
            star.wave = star.wave*(1 - 1000*star.rv0/cte.c)

        msg.info('RV0 correction is {:.3f} +/- {:.3f} km/s'.format(star.rv0, erv0))
        plt.close('all')

        if export_rv == True:
            rv_table = open(maindir+'tables/gen_ascii_RVs.txt', 'a+')
            rv_table.write('{}, {}, {:.3f}, {:.3f}\n'.format(star.id_star, star.filename, star.rv0, erv0))
            rv_table.close()

        # Correct the spectrum form cosmic rays:
        if cosmic == True:
            next_cosm = 'n'; dmin = 0.05; zs_cut = 5; niter=3; blue_cut = 4000
            while next_cosm == 'n':

                print('Current input for cosmic rays correction is zs_cut({})'.format(zs_cut)
                    +' / dmin({})'.format(dmin) +' / niter({})'.format(niter))
                change = input('To change type cut/dmin/niter, otherwise hit return: ')

                if change in ['cut','zs_cut']:
                    zs_cut = '-'
                    while type(zs_cut) is not float:
                        zs_cut = input('Choose a new zs_cut value: ')
                        if zs_cut != '': zs_cut = float(str(zs_cut).replace(',','.'))
                elif change == 'dmin':
                    dmin = '-'
                    while type(dmin) is not float:
                        dmin = input('Choose a new dmin value: ')
                        if dmin != '': dmin = float(str(dmin).replace(',','.'))
                elif change == 'niter':
                    niter = '-'
                    while type(niter) is not int:
                        niter = input('Choose a new number of iterations: ')
                        if niter != '': niter = int(niter)

                tmp_star = deepcopy(star)
                tmp_star.cosmic(method='zscore', dmin=dmin, zs_cut=zs_cut, iter=niter)

                # To prevent the noisier blue part of the spectrum to be taken for cosmic removal:
                new_blue_cut = input('Set initial wavelength from where to apply the removal (now is %s): ' % str(blue_cut))
                if new_blue_cut.isnumeric() == False and new_blue_cut != '':
                    print('Only int/float are accepted. Choosing %s... ' % str(blue_cut))
                elif new_blue_cut != '':
                    blue_cut = float(new_blue_cut)

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

                    fig_lines,ax_lines = plt.subplots(nrows, ncols, figsize=(ncols*3,nrows*1.5))
                    fig_lines.suptitle(star.filename, y=0.97, fontsize=8)

                    ax_lines = ax_lines.flatten()
                    for ax_i,line,elem in zip(ax_lines,lines,elems):
                        width = 60 if elem in ['Hdelta','Hgamma','Hbeta','Halpha'] else 15

                        mask = (star.wave >= line-width/2.) & (star.wave <= line+width/2.)
                        ax_i.plot(star.wave[mask], star.flux[mask], c='orange', lw=.7, label='RV corrected')
                        ax_i.plot(star.wave[mask], tmp_star.flux[mask], c='b', lw=.5, label='Cosmic corrected')
                        ax_i.set_title(elem, pad=0.55)
                        ax_i.tick_params(direction='in', top='on')
                        ax_i.set_yticks([])
                        if ax_i.get_ylim()[0] > 0.9:
                            ax_i.set_ylim(bottom=0.9)

                    [fig_lines.delaxes(ax_lines[i]) for i in np.arange(len(lines), len(ax_lines), 1)]

                    fig_lines.tight_layout()
                    fig_lines.subplots_adjust(wspace=.05, hspace=0.25)

                    fig_lines.show()

                next_cosm = input('Type "n" to repeat the cosmic ray removal, hit return to continue. ')
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

        finish = input('Type "n" to repeat re-do the processing, hit return to move to the next star. ')
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

    # Create the output ascii file:
    if not os.path.exists(datadir+'ASCII/POSPROC/NEW/'):
        os.makedirs(datadir+'ASCII/POSPROC/NEW/')

    star.export(output_dir=datadir+'ASCII/POSPROC/NEW/', tail='_RV', extension='.ascii')

    return None


def gen_ascii_ML(input_table='OBAs_ML_raw.fits', not_do=None, cosmic_manual=False, orig='txt'):

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

    l_wl = 3950
    r_wl = 6850

    for row in table:

        # Skip all sources in the 'not_do' input list :
        if not_do is not None and row['ID'] in not_do: continue

        if row['SNR_best'] < 60: continue
        else: print('Analysing %s' % row['ID'])

        # Determines best list for RV calculation and line for sanity check based on SpT
        if row['SpT'] < 2:
            rv_list = 'rv_Os.txt'
            line = 5411.52 # 5592.252

        elif 2 <= row['SpT'] < 2.5:
            rv_list = 'rv_Bs.txt'
            line = 4552.622

        elif 2.5 <= row['SpT'] < 2.9:
            rv_list = 'rv_Bs.txt'
            line = 6371.37

        elif row['SpT'] >= 2.9:
            rv_list = 'rv_As.txt'
            line = 4233.129

        # Determines the RV with a default fitting function and width
        fun = 'g'; wid = 15
        best_star = spec(row['ID'], snr='bestHF', orig=orig)
        best_star.rv0, erv0 = RV0(rv_list, best_star.filename, orig=orig, ewcut=50, func=fun, width=wid, tol=150)
        best_star.waveflux(l_wl-1, r_wl+1, cut_edges=True) # Applies the rv0 correction

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
        goodspec = findstar(row['ID'], snr=60)
        if len(goodspec) > 10:
            goodspec = random.sample(goodspec, 10)

        for j,n in zip(goodspec, range(len(goodspec))):
            star = spec(j.split(os.sep)[-1], orig=orig)
            star.waveflux(l_wl-1, r_wl+1, cut_edges=True)

            # Create the plot:
            fig,axs = plt.subplots(3, 1, figsize=(14,8))
            fig.suptitle(star.filename)
            axs[0].tick_params(direction='in', top='on')
            axs[0].set_xlim(l_wl-1, r_wl+1)
            axs[0].set_ylim(0.5, 1.1)
            axs[1].tick_params(direction='in', top='on')
            axs[1].set_xlim(l_wl-1, r_wl+1)
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
            star.resamp(10*0.02564975, l_wl, r_wl)

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

    l_wl = 3950
    r_wl = 6850

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
                    plt.plot([l_wl,r_wl], [2,2], c='k', lw=.5)
                    plt.show(block=False)
                    inp = input('Do you want to print data for %s it? [y/ ]: ' % file)
                    if inp == 'y':
                        print(file,'>1.5',max(data['col2']))
                    plt.close()
                if min(data['col2']) < 0: print(file,'<0',min(data['col2']))
                data.remove_columns(['col1'])
                data.write(maindir + 'tmp/new/' + file, format='ascii.no_header')
            except:
                print(file)

        else: continue

