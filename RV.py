from spec import *
from scipy.signal import correlate,correlation_lags


'''=============================================================================
Program to calculate the radial velocities for a given spectra/spectrum via
lists of individual lines or via cross correlation using synthetic spectra.
============================================================================='''

def RV_cc(id_star, windows=[(3950,4160),(4310,4360),(4370,4490),(4540,4690),(4840,4950)]):
    spectra = findstar(id_star)

    synthetic = []
    for file in os.listdir(datadir+'ASCII/Synthetic_MAUI/'):
        if id_star+'_' in file:
            synthetic.append(file)

    if len(synthetic) == 0:
        print('No files found for %s.\n' % id_star)
    elif len(synthetic) == 1:
        synthetic = synthetic[0]
    else:
        for name,i in zip(synthetic,range(len(synthetic))):
            print(name,i)
        which = input('Enter the number of the synthetic spectra you want to use: ')
        synthetic = synthetic[int(which)]

    fig, ax = plt.subplots()

    RVs_angs = []; RVs_kms = []
    for spectrum in spectra:

        print('Analyzing spectrum: ' + spectrum.split('/')[-1])

        RV_A_i,RV_kms_i = RV1_cc(spectrum,synthetic,windows=windows)
        RVs_angs.append(RV_A_i); RVs_kms.append(RV_kms_i)

        date_obs = spec(spectrum).hjd - 2400000.5
        if spectrum == spectra[0]: date_obs_0 = date_obs

        ax.scatter(date_obs,RV_kms_i,s=10,c='b')
        #ax.errorbar(date_obs,RVs_mean,yerr=std,elinewidth=.4,marker='o',
        #            color=color,capsize=2,markersize=3)

    peak_to_peak = max(RVs_kms) - min(RVs_kms)

    ax.plot([date_obs_0, date_obs],[np.mean(RVs_kms),np.mean(RVs_kms)],'-k',lw=.5)


    '''============================== Output ================================'''
    print('====================================================')
    print('Results for '+ id_star)
    print('The mean radial velocity is: ' + str(round(np.mean(RVs_kms), 4)))
    print('The peak to peak value is: ' + str(round(peak_to_peak, 4)))
    print('The time span of the spectra is: ' + str(round(date_obs-date_obs_0, 2)))
    print('The number of spectra used is: ' + str(len(RVs_kms)))
    print('====================================================')

    ax.set_title(id_star)
    ax.set_xlabel('MBJD',size=13)
    ax.set_ylabel('V$_{r}$ [km/s]',size=13)
    ax.tick_params(direction='in',top='on')
    ax.figure.subplots_adjust(top=.9,bottom=.1,right=.95,left=.15)

    completeName = os.path.join(maindir+'tmp_plots/','RVCC_%s.eps' % id_star)
    ax.figure.savefig(completeName,dpi=300)

    plt.show(block=False)


def RV1_cc(spectra1, spectra2, windows=[(3950,4160),(4310,4360),(4370,4490),(4540,4690),(4840,4950)]):
    '''

    Parameters
    ----------
    spectra1 : str
        Original spectra.

    spectra2 : str
        Synthetic spectra.

    windows : list, optional
        List of pairs of wavelengths where the cross correlation is going to be
        computed. If more than one pair, the result is averaged.
        Default is a set of defautl windows.

    Returns: XXXXXXX
    '''

    spec1 = spec(spectra1)
    spec2 = spec(spectra2,txt=True)

    resol = 1/np.sqrt((1/spec1.resolution)**2-(1/85000)**2)
    if not resol == np.inf: spec2.degrade(resol=resol)

    spec2.resamp(dx=spec1.dx,lwl=spec1.wave[0],rwl=spec1.wave[-1])
    spec1.resamp(dx=spec1.dx)

    mask = [spec2.flux < .998]

    spec1.wave = spec1.wave[mask]
    spec1.flux = spec1.flux[mask]
    spec2.wave = spec2.wave[mask]
    spec2.flux = spec2.flux[mask]

    RVs_angs = []; RVs_kms = []
    for win in windows:
        flux2 = spec2.flux[(spec2.wave >= win[0]) & (spec2.wave <= win[1])]
        flux1 = spec1.flux[(spec1.wave >= win[0]) & (spec1.wave <= win[1])]
        wave1 = spec1.wave[(spec1.wave >= win[0]) & (spec1.wave <= win[1])]

        corr = correlate(flux2-1,flux1-1)
        corr /= np.max(corr)

        # Requires scipy version 1.6.0 or above to run
        lags = correlation_lags(len(flux1),len(flux2))
        corr_shift = -lags[np.argmax(corr)]

        #lags = np.arange(-len(flux1)+1,len(flux1),1)
        #corr_shift = -corr.argmax()+len(flux1)-1

        RVs_angs.append(corr_shift*spec1.dx)
        RVs_kms.append(RVs_angs[-1]/np.mean(wave1)*cte.c/1000)

        #plt.plot(lags,corr,lw=.5)

    RV_A = round(np.mean(RVs_angs),8)
    RV_kms = round(np.mean(RVs_kms),4)

    return RV_A,RV_kms


def RV0(lines, spectrum, ewcut=50, width=20, tol=150, func='g', info=False, plot=False):
    '''
    Parameters
    ----------
    lines : str, list
        Enter the wavelenght(s) of the line(s) to fit, either in a coma-separated
        string, or in a .txt/.lst file containing the lines.

    ewcut : float, optional
        Enter the EW threshold value for a line to be used for RV. Default is 30.

    Other parameters : optional
        See help for spec and spec.fitline

    Returns: Mean radial velocity in angstroms.
    '''

    lines = findlines(lines)[0]

    RVs = []
    for line in lines:

        fit = spec(spectrum).fitline(line,width=width,tol=tol,func=func,info=info,plot=plot)

        if np.isnan(fit['RV_A']): continue
        elif fit['EW'] < ewcut: continue
        else: RVs.append(fit['RV_kms'])

    if len(RVs) == 0:
        print('\tWARNING: No lines were fitted for RV0 calculation.\n')
        return 0

    try:
        RVs = sigma_clip(RVs, sigma_lower=1.7, sigma_upper=1.7, masked=False)
        #print(RVs)
    except:
        print('Not enought values for sigma clipping. Skipping... ')
        return 0

    print('\nRV0: %s/%s lines used (after clipping) with std=%s [km/s]\n' % (len(RVs),len(lines),round(np.std(RVs),3)))

    RV_0 = np.mean(RVs)

    return RV_0


def RV(lines, spectra, SNR=None, linesRV0=None, linecut=1, ewcut=25, width=None, tol=50,\
       func='g', info=False, plot=False):
    '''

    Parameters
    ----------
    lines : str
        Enter the filename (either txt or lst), with the lines to calculate the RV.

    linesRV0 : str, optional
        Enter the filename (either txt or lst), with the lines to calculate the RV0.

    linecut : int, optional
        Enter the minimum number of lines to take the RV measurements into account.
        Default is 1.

    ewcut : float, optional
        Enter the EW threshold value for a line to be used for RV. Default is 25.

    Other parameters : optional
        See help for see spec and spec.fitline and spec() class.
    '''

    spectra = findstar(spectra=spectra,SNR=SNR)
    flines = lines
    lines,elements,_ = findlines(lines)

    fig, ax = plt.subplots()


    '''============================ PARAMETERS =============================='''
    if width == None:
        try:
            if   lines.startswith('O'):
                width = 12
                color = 'purple'
            elif lines.startswith('B'):
                width = 10
                color = 'b'
            elif lines.startswith('A'):
                width = 10
                color = 'teal'
            elif lines.startswith('M'):
                width = 12
                color = 'r'
            else:
                width = 10
                color = 'g'
        except:
            width = 10
            color = 'g'


    '''=============================== SPECTRA =============================='''
    i = 0
    RVs_all_means = []
    RVs_all_std = []
    num_lines = []
    for spectrum in spectra:

        print('\n##########################################################')
        print('Analyzing spectrum: ' + spectrum.split('/')[-1])

        if linesRV0 == None:
            RV_0 = 0
        else:
            RV_0 = RV0(linesRV0,spectrum,ewcut=ewcut,width=width,func=func)

        spectrum = spec(spectrum,offset=RV_0)

        date_obs = spectrum.hjd - 2400000.5
        if i == 0:
            date_obs_0 = date_obs
            out_f1 = open(maindir+'radial_velocity/%s.csv' % (spectrum.id_star+'_RV'),'a')
            out_f1.write('mean_rv, mean_rv_error, mbjd, num_lines, spectrum\n')

        out_f2 = open(maindir+'radial_velocity/%s.csv' % (spectrum.id_star+'_lines'),'a')
        out_f2.write(spectrum.id_star +' | '+ spectrum.filename +'\n')
        out_f2.write('target_line, fitted_line, RV, EW, FWHM, q_fit\n')

        RVs = []
        for line in lines:

            fit = spectrum.fitline(line,width=width,tol=tol,func=func,info=info,plot=plot)

            if np.isnan(fit['RV_kms']):
                continue

            elif fit['EW'] < ewcut:
                continue

            else:
                RVs.append(fit['RV_kms'])
                out_f2.write('%.3f, %.3f, %d, %d, %.2f, %.3f\n' %
                    line,fit['line'],fit['RV_kms'],fit['EW'],fit['FWHM'],fit['q_fit'])

        if RVs == []:
            print('No lines were found for spectrum: %s\n' % spectrum.filename)
            continue
        if len(RVs) < linecut:
            print('Only %d line found for: %s\n' % (len(RVs),spectrum.filename))
            continue


        '''========================= Sigma Clipping ========================='''
        # Enable/disable histogram (also enable also the lines after sigma clip)
        #ax.hist(RVs, bins = 60, alpha = 0.8)
        RVs = sigma_clip(RVs,sigma_lower=2,sigma_upper=2,masked=False, \
            return_bounds=True,cenfunc='median')
        #ax.plot([np.median(RVs[0]),np.median(RVs[0])],[0,10],'-k') # Prints central value
        #ax.plot([RVs[1],RVs[1]],[0,10],'--r') # Prints low  clipping value
        #ax.plot([RVs[2],RVs[2]],[0,10],'--r') # Prints high clipping value

        '''=================== Getting the important data ==================='''
        num_lines.append(len(RVs[0]))
        RVs_mean = np.mean(RVs[0])                      # Mean for a spectrum
        std = np.std(RVs[0])/np.sqrt(num_lines[i])      # std  for a spectrum
        RVs_all_means.append(RVs_mean)
        RVs_all_std.append(std)

        out_f1.write(
            str(round(RVs_mean, 4)) + ', ' + str(round(std, 4)) + ', ' \
            +  str(round(date_obs,4)) + ', ' + str(len(RVs[0])) + ', ' \
            +  spectrum.filename.split('.')[0] + '\n')


        '''============================== Plot =============================='''
        ax.errorbar(date_obs,RVs_mean,yerr=std,elinewidth=.4,marker='o',
            color=color,capsize=2,markersize=3)

        i = i + 1

    out_f1.close()
    out_f2.close()

    if i == 0:
        print('The program did not work, check the output lines...')
        plt.close()
        return None


    '''==================== Mean of RVs and peak-to-peak ===================='''
    RVs_mean_all_means = np.mean(RVs_all_means)    # Mean of all spectra
    RVs_mean_all_means_std = np.std(RVs_all_means) # Std of the mean of all spectra
    ax.plot([date_obs_0, date_obs],[RVs_mean_all_means,RVs_mean_all_means],'-k',lw=.5)
    peak_to_peak = abs(max(RVs_all_means) - min(RVs_all_means))
    peak_to_peak_err = np.sqrt(RVs_all_std[RVs_all_means.index(max(RVs_all_means))]**2
        + RVs_all_std[RVs_all_means.index(min(RVs_all_means))]**2)


    '''============================== Output ================================'''
    print('====================================================')
    print('Results for '+ spectrum.id_star)
    if i == 1:
        print('The mean RV of the spectrum is: ' + str(round(RVs_mean,4)) \
        + ' +/- ' + str(round(std,4)))
    else:
        print('The mean radial velocity is: ' + str(round(RVs_mean_all_means, 4)) \
        + ' +/- ' + str(round(RVs_mean_all_means_std, 4)))
        print('The peak to peak value is: ' + str(round(peak_to_peak, 4)) \
        + ' +/- ' + str(round(peak_to_peak_err,4)))
        print('The time span of the spectra is: ' + str(round(date_obs-date_obs_0, 2)))
    print('The average number of lines used is: ' + str(round(np.mean(num_lines), 1)))
    print('The number of spectra used is: ' + str(len(RVs_all_means)))
    print('====================================================')


    '''================================ Plot ================================'''
    ax.set_title(spectrum.id_star)
    ax.set_xlabel('MBJD',size=13)
    ax.set_ylabel('V$_{r}$ [km/s]',size=13)
    ax.tick_params(direction='in',top='on')
    ax.figure.subplots_adjust(top=.9,bottom=.1,right=.95,left=.15)

    completeName = os.path.join(maindir+'tmp_plots/','%s.eps' % ('RV_'+spectrum.id_star))
    ax.figure.savefig(completeName,dpi=300)

    plt.show(block=False)
