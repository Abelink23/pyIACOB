'''=============================================================================
Program to calculate the radial velocities for a given spectra/spectrum via
lists of individual lines or via cross correlation using synthetic spectra.
============================================================================='''

from spec import *

from scipy.signal import correlate,correlation_lags

import random


def RV0_cc(spec1, spec2, orig1='IACOB', orig2='synthetic', method='windows', 
           lwl=3800, rwl=8000, show_plot=False,
           windows=[(3950,4160),(4310,4360),(4370,4490),(4540,4690),(4840,4950)]):

    '''
    Function to obtain the radial velocity of a spectrum using a cross correlation technique
    relative to a second reference spectrum of the same star.

    If you compare a IACOB spectrum with a synthetic one, use the 'IACOB' as the origin of the
    first spectrum and 'synthetic' as the origin of the second one.

    NOTE 1: This function is very sensitive if you include a part of the spectrum with much
            lower S/N ratio (e.g., below 4000A). It is recommended to cut this part.

    NOTE 2: The 2*1/snr1 is used to mask the continuum and to vary the flux of the spectra
            in the Monte Carlo simulation. The factor *2* has been chosen based on testing.

    NOTE 3: As it is now, it is spec1 the spectrum which is being degraded to the resolution
            of spec2. This might be changed in the future.

    Parameters
    ----------
    spec1 : str
        Original spectrum.

    spec2 : str
        Spectrum to compare (reference).

    orig1 : str, optional
        Select the origin of the reference spectrum (see spec() for more information).
        Default is 'IACOB' for reference and'synthetic' for the comparison.

    orig2 : str, optional
        Select the origin of the spectrum to compare (see spec() for more information).
        Default is 'synthetic'.

    method : str, optional
        Options are 'windows', 'simple' and 'mcmc'. Default is 'windows'.

    lwl, rwl : float, optional
        Left and right wavelength limits of the spectra.
        Default is 3800 and 8000, respectively.
    
    show_plot : boolean, optional
        True if you want to see the parts of the spectra used for the cross correlation.

    windows : list, optional
        List of pairs of wavelengths where the cross correlation is going to be
        computed. If more than one pair, the result is averaged.
        Default is a set of default windows.

    Returns
    -------
    Radial velocity in km/s and Angstroms, and their respective uncertainties.
    '''

    spec1 = spec(spec1, orig=orig1)
    spec2 = spec(spec2, orig=orig2)

    spec1.waveflux(lwl, rwl)
    spec2.waveflux(lwl, rwl)
    
    snr1 = spec1.snrcalc()
    snr2 = spec2.snrcalc()

    # Change the resolution of the spectra if needed
    if orig1 == 'IACOB' and orig2 in ['synthetic','syn']:
        print('Changing the resolution of the synthetic spectrum to the one of spectrum 1.')
        spec2.degrade(resol=spec1.resolution)
        #resol = 1/np.sqrt((1/spec1.resolution)**2-(1/85000)**2)
        #if not resol == np.inf:
        #    spec2.degrade(resol=resol)
    if orig1 == 'IACOB' and orig2 == 'IACOB':
        if spec1.resolution > spec2.resolution:
            print('The spectrum 1 will be degraded to the resolution of spectrum 2.')
            spec1.degrade(resol=spec2.resolution)
        elif spec1.resolution < spec2.resolution:
            print('The spectrum 2 will be degraded to the resolution of spectrum 1.')
            spec2.degrade(resol=spec1.resolution)

    # Resample the spectra to have the same number of points
    if spec1.wave[-1]-spec1.wave[0] > spec2.wave[-1]-spec2.wave[0]:
        spec1.waveflux(spec2.wave[0], spec2.wave[-1])
        spec2.resamp(dlam=spec1.dlam, lwl=spec1.wave[0], rwl=spec1.wave[-1], method='linear')
        spec1.resamp(dlam=spec1.dlam) # To have the same number of points in the CC

    else:
        spec2.waveflux(spec1.wave.min(), spec1.wave.max())
        spec1.resamp(dlam=spec1.dlam, lwl=spec2.wave[0], rwl=spec2.wave[-1])
        spec2.resamp(dlam=spec1.dlam)

    # use the S/N ratio to define the mask where to evaluate the cross-correlation
    if orig1 == 'IACOB' and orig2 != 'IACOB':
        print('The S/N ratio of spectrum 1 is used in spectrum 2 to mask the continuum.')
        mask = spec2.flux < 1-2*1/snr1
    elif orig1 == 'IACOB' and orig2 == 'IACOB':
        if snr1 < snr2:
            print('The S/N ratio of spectrum 1 is used in spectrum 1 to mask the continuum.')
            mask = spec1.flux < 1-2*1/snr1
        else:
            print('The S/N ratio of spectrum 2 is used in spectrum 2 to mask the continuum.')
            mask = spec2.flux < 1-2*1/snr2
    else:
        print('A maximum flux of 0.998 is used in spectrum 2 to mask the continuum.')
        mask = spec2.flux < .998

    # mask out all the points between 5885 and 5900
    mask = mask & ((spec1.wave < 5885) | (spec1.wave > 5900))
    # mask out all the points between 6265 and 6330
    mask = mask & ((spec1.wave < 6275) | (spec1.wave > 6330))
    # mask out all the points between 6532.8 and 6592.8 (Halpha)
    mask = mask & ((spec1.wave < 6532.8) | (spec1.wave > 6592.8))
    # mask out all the points between 6850 and 7430 in all cases
    mask = mask & ((spec1.wave < 6850) | (spec1.wave > 7430))
    # mask out all the points between 7585 and 7720
    mask = mask & ((spec1.wave < 7585) | (spec1.wave > 7720))

    if show_plot == True:
        fig, ax = plt.subplots(1, 2, width_ratios=[3, 1], figsize=(12, 4), constrained_layout=True)
        ax[0].plot(spec1.wave[~mask], spec1.flux[~mask], 'r', lw=2, alpha=.4)
        ax[0].plot(spec1.wave, spec1.flux, 'b', lw=.5)
        ax[0].plot(spec2.wave[~mask], spec2.flux[~mask], 'r', lw=2, alpha=.4)
        ax[0].plot(spec2.wave, spec2.flux, 'g', lw=.5)
        ax[0].set_xlabel('Wavelength [$\AA$]', size=10)

    spec1.wave = spec1.wave[mask]
    spec1.flux = spec1.flux[mask]

    spec2.wave = spec2.wave[mask]
    spec2.flux = spec2.flux[mask]

    if method == 'windows':
        RVs_A = []; RVs_kms = []
        for win in windows:
            mask_widw = (spec1.wave >= win[0]) & (spec1.wave <= win[1])

            flux1 = spec1.flux[mask_widw]
            wave1 = spec1.wave[mask_widw]

            if len(wave1) == 0:
                print('No points within the window %s-%s. Skipping...' % (win[0],win[1]))
                continue

            flux2 = spec2.flux[mask_widw]

            corr = correlate(flux1-1, flux2-1)
            corr /= np.max(corr)

            lags = correlation_lags(len(flux1), len(flux2))*spec1.dlam
            RV_A_i = lags[np.argmax(corr)]

            # Obtain the radial velocity in Angstroms
            RVs_A.append(RV_A_i)
            RVs_kms.append(RV_A_i/np.mean(wave1)*cte.c/1000)

            if show_plot == True:
                # plot a shaded region with the windows used for the cross-correlation
                ax[0].axvspan(win[0], win[1], alpha=0.2, color='gray')
                # plot the cross-correlation function
                ax[1].plot(lags, corr, lw=.5)

        RV_A = round(np.mean(RVs_A), 6)
        e_RV_A = round(np.std(RVs_A), 6)
        RV_kms = round(np.mean(RVs_kms), 4)
        e_RV_kms = round(np.std(RVs_kms), 4)

    elif method == 'simple':
        corr = correlate(spec1.flux-1, spec2.flux-1)
        corr /= np.max(corr)

        lags = correlation_lags(len(spec1.flux), len(spec2.flux))*spec1.dlam

        # Obtain the radial velocity in Angstroms and km/s
        RV_A = lags[np.argmax(corr)]
        # define the uncertainty as the half-width of the correlation function at 0.9 from its maximum
        e_RV_A = 1/2*np.abs(lags[np.where(corr >= 0.99)[0][-1]] - lags[np.where(corr >= 0.99)[0][0]])
        RV_kms = RV_A/np.mean(spec1.wave)*cte.c/1000
        e_RV_kms = e_RV_A/np.mean(spec1.wave)*cte.c/1000

        if show_plot == True:
            # plot the uncertainty of the cross-correlation
            ax[1].plot([RV_A-e_RV_A, RV_A+e_RV_A], [0.99, 0.99], c='k', marker='x', lw=0.5)
            # plot the cross-correlation function
            ax[1].plot(lags,corr,lw=.5)

    elif method == 'mcmc':

        RVs_A = []; RVs_kms = []; mcmc_n = 500

        print('Running a Monte Carlo simulation with %d iterations...' % mcmc_n)

        for i in range(mcmc_n):
            # use snr2 to vary spec1.flux (with noise and continuum) for the Monte Carlo simulation
            if orig1 in ['synthetic','syn']:
                flux1_i = spec1.flux-1 + np.random.normal(0, 1/snr2, len(spec1.flux))
                flux1_i += np.random.normal(0, 1/snr2, 1)
            else:
                flux1_i = spec1.flux-1

            # use snr1 to vary spec2.flux (with noise and continuum) for the Monte Carlo simulation
            if orig2 in ['synthetic','syn']:
                flux2_i = spec2.flux-1 + np.random.normal(0, 1/snr1, len(spec2.flux))
                flux2_i += np.random.normal(0, 1/snr1, 1)
            else:
                flux2_i = spec2.flux-1

            if orig1 in ['synthetic','syn'] and orig2 in ['synthetic','syn']:
                print('Both spectra are synthetic. This is not possible.')
                return None

            #// plot the individual modified spectra (FOR DEVELOPMENT)
            #//if show_plot == True:
            #//    ax[0].plot(spec1.wave, flux1_i+1, 'b', lw=.5)
            #//    ax[0].plot(spec2.wave, flux2_i+1, 'g', lw=.5)

            corr = correlate(flux1_i, flux2_i)
            corr /= np.max(corr)

            lags = correlation_lags(len(flux1_i), len(flux2_i))*spec1.dlam

            # Obtain the radial velocity in Angstroms and km/s
            RV_A_i = lags[np.argmax(corr)]
            RVs_A.append(RV_A_i)
            RVs_kms.append(RV_A_i/np.mean(spec1.wave)*cte.c/1000)

        if show_plot == True:
            # plot the cross-correlation function in units of Angstroms
            ax[1].plot(lags, corr, lw=.5)

        # Obtain the mean and standard deviation of the radial velocity in Angstroms and km/s
        RV_A = round(np.mean(RVs_A),6)
        e_RV_A = round(np.std(RVs_A),6)
        RV_kms = round(np.mean(RVs_kms),4)
        e_RV_kms = round(np.std(RVs_kms),4)

    # If the uncertainty is lower than the step in lambda, set it to that value
    if e_RV_A < spec1.dlam:
        e_RV_A = round(spec1.dlam, 6)
        e_RV_kms = round(spec1.dlam/np.mean(spec1.wave)*cte.c/1000, 4)

    if show_plot == True:
        # draw the zero line and the maximum of the cross-correlation function
        ax[1].axvline(0, color='gray', lw=.5, ls='--')
        ax[1].axvline(RV_A, color='k', lw=0.5)
        xlim = 10*RV_A if abs(RV_A) < 0.2 else 4
        ax[1].set_xlim(-xlim, xlim)
        ax[1].set_ylim(bottom=0.5)
        ax[1].set_xlabel('RV [$\AA$]', size=10)
        ax[1].text(0.05, 0.15, 'RV=%.3f [$\AA$]' % RV_A, transform=ax[1].transAxes, fontsize=8)
        ax[1].text(0.05, 0.10, 'RV=%.2f [km/s]' % RV_kms, transform=ax[1].transAxes, fontsize=8)

        plt.show(block=False)

    return RV_kms, e_RV_kms, RV_A, e_RV_A


def RV_cc(id_star, snr=0, n_max=50, orig='IACOB', method='windows', lwl=3800, rwl=8000, 
          windows=[(3950,4160),(4310,4360),(4370,4490),(4540,4690),(4840,4950)]):

    '''
    Function to obtain the radial velocities of a star using a cross correlation technique
    relative to a reference spectrum of the same star.

    The program will try to find a synthetic spectrum for the given star when 'IACOB' is used
    in the orig parameter.

    Parameters
    ----------
    id_star : str
        Star ID to analyze the radial velocity.

    snr : float, optional
        Enter the minimum S/N ratio to search for the spectra. Default is 0.

    n_max : int, optional
        Enter the maximum number of spectra to be used. Default is 50.

    orig : str, optional
        Select the origin of the reference spectrum [IACOB/synthetic]. Default is 'IACOB'.

    For the rest of the parameters, see RV0_cc() function.

    Returns
    -------
    Radial velocity curve.
    '''

    spectra = findstar(id_star, snr=snr)
    spectra = [i.split('/')[-1] for i in spectra]

    if len(spectra) > n_max:
        n = input('Number of spectra is more then %i, do you want to take a random number of them? [#/no/n]: ' % n_max)
        if n not in ['no','n','']:
            spectra = random.sample(spectra, int(n))
        elif n == '':
            print('Plotting %i random spectra...' % n_max)
            spectra = random.sample(spectra, n_max)

    synthetic = []
    for root, dirs, files in os.walk(datadir+'ASCII/Synthetic_MAUI/'):
        for file in files:
            if id_star+'_' in file:
                synthetic.append(os.path.join(root, file))

    if len(synthetic) == 0:
        print('No synthetic files found for %s.\n' % id_star)
        return None
    elif len(synthetic) == 1:
        synthetic = synthetic[0]
    else:
        for name,i in zip(synthetic,range(len(synthetic))):
            print(name,i)
        which = input('Enter the number of the synthetic spectra you want to use: ')
        synthetic = synthetic[int(which)]

    synthetic = synthetic.split('/')[-1]

    if not os.path.isfile(maindir+'radial_velocity/summary.txt'):
        out_f0 = open(maindir+'radial_velocity/summary.txt', 'a')
        out_f0.write('ID, RV_all, RV_all_std, RV_p2p, RV_p2p_err, Tspan, Nspec\n')
    else:
        out_f0 = open(maindir+'radial_velocity/summary.txt', 'a')

    out_f1 = open(maindir+'radial_velocity/%s.csv' % (id_star+'_RV'), 'a')
    out_f1.write('rv, rv_error, mbjd, spectrum\n')

    fig, ax = plt.subplots(constrained_layout=True)

    i = 0
    RVs_kms = []; e_RVs_kms = []
    for spectrum in spectra:
        spectrum = spectrum.split('/')[-1]

        print('\n##########################################################')
        print('Analyzing spectrum: ' + spectrum)

        print(spectrum, synthetic)
        RV_kms_i, e_RV_kms_i, _,_ = RV0_cc(spectrum, synthetic, orig1=orig, orig2='syn', 
                                    method=method, lwl=lwl, rwl=rwl, windows=windows)

        spectrum = spec(spectrum, rv0=RV_kms_i)

        date_obs = spectrum.hjd - 2400000.5
        if i == 0:
            date_obs_0 = date_obs

        '''=================== Getting the important data ==================='''
        RVs_kms.append(RV_kms_i)
        e_RVs_kms.append(e_RV_kms_i)

        out_f1.write('%.4f, %.4f, %.4f, %s\n' %
            (RV_kms_i,e_RV_kms_i,date_obs,spectrum.filename.split('.')[0]))

        '''============================== Plot =============================='''
        ax.errorbar(date_obs, RV_kms_i, yerr=e_RV_kms_i, elinewidth=.4, marker='o',
                    color='b',capsize=2,markersize=3)

        i = i + 1

    out_f1.close()

    if i == 0:
        print('The program did not work, check the output lines...')
        plt.close()
        return None

    '''==================== Mean of RVs and peak-to-peak ===================='''
    RVs_mean_all = np.mean(RVs_kms)    # Mean of all spectra
    std_RVs_mean_all = np.std(RVs_kms) # Std of the mean of all spectra
    if i == 1:
        peak2peak = peak2peak_err = tspan = 0
    else:
        ax.plot([date_obs_0, date_obs], [RVs_mean_all,RVs_mean_all], '-k', lw=.5)
        peak2peak = abs(max(RVs_kms) - min(RVs_kms))
        peak2peak_err = np.sqrt(e_RVs_kms[RVs_kms.index(max(RVs_kms))]**2
            + e_RVs_kms[RVs_kms.index(min(RVs_kms))]**2)
        tspan = date_obs-date_obs_0

    '''============================== Output ================================'''
    print('====================================================')
    print('Results for '+ id_star)
    print('The mean radial velocity is: %.4f +/- %.4f' % (RVs_mean_all, std_RVs_mean_all))
    print('The peak to peak value is: %.4f +/- %.4f' % (peak2peak, peak2peak_err))
    print('The time span of the spectra is: %d' % tspan)
    print('The number of spectra used is: %d' % i)
    print('====================================================')

    out_f0.write('%s, %.4f, %.4f, %.4f, %.4f, %.2f, %d\n' %
        (id_star, RVs_mean_all, std_RVs_mean_all, peak2peak, peak2peak_err,
        tspan, i))

    out_f0.close()

    '''================================ Plot ================================'''
    ax.set_title(id_star, size=10)
    ax.set_xlabel('MBJD',size=10)
    ax.set_ylabel('V$_{r}$ [km/s]',size=10)
    ax.tick_params(direction='in',top='on')

    completeName = os.path.join(maindir+'tmp_plots/', 'RVcc_'+id_star+'.png')
    fig.savefig(completeName, format='png', dpi=200, bbox_inches='tight')

    plt.show(block=False)


def RV0(lines, spectrum, orig='IACOB', ewcut=50, width=20, tol=150, func='g', check_fit=False, plot=False):

    '''
    Function to calculate the radial velocity of a given spectrum using a set of input
    lines where the individual RV is measured and then a sigma-clipping and final average
    is used to compute it.

    Parameters
    ----------
    lines : str, list
        Enter the wavelength(s) of the line(s) to fit, either in a coma-separated
        string, or in a .txt/.lst file containing the lines.

    spectrum : str
        Enter the filename of the spectrum.

    orig : str, optional
        See spec() for more information. Default is 'IACOB'.

    ewcut : float, optional
        Enter the EW threshold value for a line to be used for RV. Default is 30.

    check_fit : boolean, optional
        True if you want to see the individual information of each fitting and discard
        potential bad fittings from a plot. Note: this set the plot option to True.
        Default is False.

    plot : boolean, optional
        True if you want to see a plot with the individual line fittings. Default is False.

    Other parameters : optional
        See help for spec and spec.fitline

    Returns
    -------
    Mean radial velocity in km/s.
    '''

    lines = findlines(lines)[0]

    star = spec(spectrum, orig=orig)

    RVs = []; i = 0
    for line in lines:

        fit = star.fitline(line, width=width, tol=tol, func=func, info=check_fit, outfit=True)

        if np.isnan(fit['RV_kms']):
            continue
        elif fit['EW'] < ewcut:
            continue
        else:
            RVs.append(fit['RV_kms'])

        if plot == True or check_fit == True:

            if i == 0:
                fig, axs = plt.subplots(1,len(lines), tight_layout=True, figsize=(16,2))
                fig.subplots_adjust(wspace=0, hspace=0)
                axs = axs.flatten()

            axs[i].plot(fit['wave'], fit['flux_norm'], c='b', lw=.5)
            axs[i].plot(fit['wave'], fit['flux_fit'], c='g', lw=.5)
            axs[i].set_title('Line: ' + str(fit['line'])
                + ' - RV[A/Kms]= ' + str(fit['RV_A']) + '/' + str(fit['RV_kms']), fontsize=5)
            axs[i].set_yticks([]); axs[i].set_xticks([])

            i += 1

    if len(RVs) == 0:
        print('\tWARNING: No lines were fitted for RV0 calculation.\n')
        return 0,0

    try:
        RVs_f = sigma_clip(RVs, sigma_lower=1.7, sigma_upper=1.7, masked=False)
        idx_bad = [j for j in range(len(RVs)) if RVs[j] not in RVs_f]
        #print(RVs)
    except:
        print('Not enought values for sigma clipping. Skipping... ')
        return 0,0

    if plot == True or check_fit == True:
        # Remove plots with failed fittings
        [fig.delaxes(axs[j]) for j in np.arange(i, len(axs), 1)]
        # Remove plots with distarded fittings after the sigma clipping
        [fig.delaxes(axs[j]) for j in idx_bad if idx_bad != []]

        # Set the minimum y-value of all the plots based on the global minimum
        ymin = np.min([axs[j].get_ylim()[0] for j in range(i) if i != 0 and not j in idx_bad])
        [axs[j].set_ylim(bottom=ymin) for j in range(i) if i != 0 and not j in idx_bad]

        plt.show(block=False)

        if check_fit == True:
            remove = input('Which lines from the plotted ones you want to remove'\
            '(e.g. 0,3). Hit return to continue.\n').split(',')
            remove = [int(j) for j in remove if not remove == ['']]

            RVs_f = [RVs_f[j] for j in range(len(RVs_f)) if not j in remove]

    RV_0 = np.mean(RVs_f)
    e_RV_0 = np.std(RVs_f)

    print('\nRV0=%s (%s/%s lines used with std=%s [km/s])' %
        (round(RV_0,2),len(RVs_f),len(lines),round(e_RV_0,2)))

    return RV_0, e_RV_0


def RV(lines, id_star, snr=0, linesRV0=None, n_max=50, linecut=1, ewcut=25, width=None, tol=50,\
       func='g', info=False, plot=False):

    '''
    Function to obtain RV measurements for a given multi-epoch data.

    Parameters
    ----------
    lines : str
        Enter the filename (either txt or lst), with the lines to calculate the RV.

    id_star : str, list
        Enter the name of the star, or list of fits-files to analyze.

    snr : float, optional
        Enter the minimum S/N ratio to search for the spectra. Default is 0.

    linesRV0 : str, optional
        Enter the filename (either txt or lst), with the lines to calculate the RV0.

    n_max : int, optional
        Enter the maximum number of spectra to be used. Default is 50.

    linecut : int, optional
        Enter the minimum number of lines to take the RV measurements into account.
        Default is 1.

    ewcut : float, optional
        Enter the EW threshold value for a line to be used for RV. Default is 25.

    width : float, optional
        Enter the width of the fitting window. If None, it will be set automatically.

    Other parameters : optional
        See help for see spec and spec.fitline and spec() class.

    Returns
    -------
    Nothing, but output text files are created with the individual and global results.
    Terminal output is also prompted.
    
    Notes
    -----
    The tolerance of the initial RV0 calculation is set to 3 times the tolerance
    of the individual line fitting.
    '''


    '''============================ PARAMETERS =============================='''
    if width == None:
        if   lines.startswith('O'):
            width = 20
            color = 'purple'
        elif lines.startswith('B'):
            width = 15
            color = 'b'
        elif lines.startswith('A'):
            width = 10
            color = 'teal'
        elif lines.startswith('M'):
            width = 10
            color = 'r'
        else:
            width = 15
            color = 'g'
    else:
        width = 15
        color = 'g'


    '''=============================== SPECTRA =============================='''
    spectra = findstar(spectra=id_star, snr=snr)

    if len(spectra) > n_max:
        n = input('Number of spectra is more then %i, do you want to take a random number of them? [#/no/n]: ' % n_max)
        if n not in ['no','n','']:
            spectra = random.sample(spectra, int(n))
        elif n == '':
            print('Plotting %i random spectra...' % n_max)
            spectra = random.sample(spectra, n_max)

    lines,elements,_ = findlines(lines)

    if not os.path.isfile(maindir+'radial_velocity/summary.txt'):
        out_f0 = open(maindir+'radial_velocity/summary.txt', 'a')
        out_f0.write('ID, RV_all, RV_all_std, RV_p2p, RV_p2p_err, Tspan, Nlines, Nspec\n')
    else:
        out_f0 = open(maindir+'radial_velocity/summary.txt', 'a')

    out_f1 = open(maindir+'radial_velocity/%s.csv' % (id_star+'_RV'), 'a')
    out_f1.write('rv, rv_error, mbjd, num_lines, spectrum\n')

    fig, ax = plt.subplots(constrained_layout=True)

    i = 0
    RVs_all_means = []
    std_RVs_all = []
    num_lines = []
    for spectrum in spectra:
        spectrum = spectrum.split('/')[-1]

        print('\n##########################################################')
        print('Analyzing spectrum: ' + spectrum)

        if linesRV0 == None:
            RV_0 = 0
        else:
            RV_0, eRV_0 = RV0(linesRV0, spectrum, ewcut=ewcut, width=width, func=func, tol=3*tol)

        spectrum = spec(spectrum, rv0=RV_0)

        date_obs = spectrum.hjd - 2400000.5
        if i == 0:
            date_obs_0 = date_obs

        out_f2 = open(maindir+'radial_velocity/%s.csv' % (spectrum.id_star+'_lines'), 'a')
        out_f2.write(spectrum.id_star +' | '+ spectrum.filename +'\n')
        out_f2.write('target_line, fitted_line, RV, EW, FWHM, q_fit\n')

        RVs = []
        for line in lines:

            fit = spectrum.fitline(line, width=width, tol=tol, func=func, info=info, plot=plot)

            if np.isnan(fit['RV_kms']):
                continue
            elif fit['EW'] < ewcut:
                continue
            else:
                RV_i = fit['RV_kms'] + RV_0
                RVs.append(RV_i)
                out_f2.write('%.3f, %.3f, %.2f, %d, %.2f, %.3f\n' %
                    (line,fit['line'],RV_i,fit['EW'],fit['FWHM'],fit['q_fit']))

        if RVs == []:
            print('No lines were found for spectrum: %s\n' % spectrum.filename)
            continue

        if len(RVs) < linecut:
            print('Only %d line found for: %s\n' % (len(RVs),spectrum.filename))
            continue

        '''========================= Sigma Clipping ========================='''
        # Enable/disable histogram (also enable also the lines after sigma clip)
        #ax.hist(RVs, bins = 60, alpha = 0.8)
        RVs = sigma_clip(RVs, sigma_lower=2, sigma_upper=2, masked=False, \
            return_bounds=True, cenfunc='median')
        #ax.plot([np.median(RVs[0]),np.median(RVs[0])],[0,10],'-k') # Prints central value
        #ax.plot([RVs[1],RVs[1]],[0,10],'--r') # Prints low  clipping value
        #ax.plot([RVs[2],RVs[2]],[0,10],'--r') # Prints high clipping value

        '''=================== Exporting the important data ================='''
        num_lines.append(len(RVs[0]))
        RVs_mean = np.mean(RVs[0])                      # Mean for a spectrum
        std = np.std(RVs[0])/np.sqrt(num_lines[i])      # std  for a spectrum
        RVs_all_means.append(RVs_mean)
        std_RVs_all.append(std)

        out_f1.write('%.4f, %.4f, %.4f, %d, %s\n' %
            (RVs_mean,std,date_obs,len(RVs[0]),spectrum.filename.split('.')[0]))

        '''============================== Plot =============================='''
        ax.errorbar(date_obs, RVs_mean, yerr=std, elinewidth=.4, marker='o', 
                    color=color, capsize=2, markersize=3)

        i = i + 1

    out_f1.close()
    out_f2.close()

    if i == 0:
        print('The program did not work, check the output lines...')
        plt.close()
        return None

    '''==================== Mean of RVs and peak-to-peak ===================='''
    RVs_mean_all = np.mean(RVs_all_means)    # Mean of all spectra
    std_RVs_mean_all = np.std(RVs_all_means) # Std of the mean of all spectra
    if i == 1:
        peak2peak = peak2peak_err = tspan = 0
    else:
        ax.plot([date_obs_0, date_obs], [RVs_mean_all,RVs_mean_all], '-k', lw=.5)
        peak2peak = abs(max(RVs_all_means) - min(RVs_all_means))
        peak2peak_err = np.sqrt(std_RVs_all[RVs_all_means.index(max(RVs_all_means))]**2
            + std_RVs_all[RVs_all_means.index(min(RVs_all_means))]**2)
        tspan = date_obs-date_obs_0

    '''============================== Output ================================'''
    print('====================================================')
    print('Results for '+ id_star)
    print('The mean radial velocity is: %.4f +/- %.4f' % (RVs_mean_all, std_RVs_mean_all))
    print('The peak to peak value is: %.4f +/- %.4f' % (peak2peak, peak2peak_err))
    print('The time span of the spectra is: %d' % tspan)
    print('The average number of lines used is: %.1f' % np.mean(num_lines))
    print('The number of spectra used is: %d' % i)
    print('====================================================')

    out_f0.write('%s, %.4f, %.4f, %.4f, %.4f, %.2f, %.1f, %d\n' %
        (id_star, RVs_mean_all, std_RVs_mean_all, peak2peak, peak2peak_err,
        tspan, np.mean(num_lines), i))

    out_f0.close()

    '''================================ Plot ================================'''
    ax.set_title(id_star, size=10)
    ax.set_xlabel('MBJD',size=10)
    ax.set_ylabel('V$_{r}$ [km/s]',size=10)
    ax.tick_params(direction='in',top='on')

    completeName = os.path.join(maindir+'tmp_plots/', 'RV_'+id_star+'.png')
    fig.savefig(completeName, format='png', dpi=200, bbox_inches='tight')

    plt.show(block=False)
