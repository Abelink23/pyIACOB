# Python Standard Library packages
import os

# Main packages
import numpy as np
import matplotlib.pyplot as plt

# Astro-packages
import astropy.units as u
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation

# Scientific packages:
import scipy.constants as cte
from scipy.signal import convolve
from scipy.interpolate import interp1d


def raw_to_IACOB(path_to_spectra, table_with_spc=None, plot=False):

    '''
    Function to turn raw files from the telescope to the IACOB format.

    Parameters
    ----------
    path_to_spectra : str
        Enter the path of the folder containing the files to convert.

    Returns
    -------
    Nothing, but the files are reduced and saved in the IACOB format.
    '''

    spectra = [os.path.join(path_to_spectra, i) for i in os.listdir(path_to_spectra)
                if i.endswith('.fits') and not i.startswith('.')]

    for spectrum in spectra:
        # Retrieve the key values from the fits header
        hdu = fits.open(spectrum) # Open the fits image file
        hdu0 = hdu[0]           # Load the header list of primary header
        header = hdu0.header    # Read the values of the headers

        object = header['OBJECT'].replace(' ','').upper()
        DATE = header['DATE-AVG']
        date, time = DATE.split('T')
        date = date.replace('-','')
        time = time.replace(':','').split('.')[0]

        if header['INSTRUME'].strip() == 'FIES':
            telcode = '_N_'
            TELES = 'NOT' ; TELESINFO = 'Telescope (www.not.iac.es)'
            INSTR = 'FIES'
            #loc = EarthLocation.of_site('LaPalma')
            loc = EarthLocation.from_geodetic(lat=28.7624*u.deg, lon=-17.8792*u.deg, height=2350*u.m)
            if header['FIFMSKNM'].strip() == 'F4 HiRes':
                RESOL = 'V67000'
            elif header['FIFMSKNM'].strip() == 'F3 MedRes':
                RESOL = 'V46000'
            elif header['FIFMSKNM'].strip() == 'F1 LowRes':
                RESOL = 'V25000'
        elif header['INSTRUME'].strip() == 'HERMES':
            telcode = '_M_'
            TELES = 'MERCATOR'; TELESINFO = 'Telescope (www.mercator.iac.es)'
            INSTR = 'HERMES'
            RESOL = 'V85000'
            loc = EarthLocation.from_geodetic(lat=28.7624*u.deg, lon=-17.8792*u.deg, height=2350*u.m)
        elif header['INSTRUME'].strip() == 'FEROS':
            telcode = '_F_'
            TELES = 'MPI-2.2' ; TELESINFO = 'Telescope (www.mpia-hd.mpg.de/FEROS)'
            INSTR = 'FEROS'
            RESOL = 'V48000'
            #loc = EarthLocation.of_site('LaSilla')
            loc = EarthLocation.from_geodetic(lat=-29.2567*u.deg, lon=-70.7366*u.deg, height=2350*u.m)

        sc = SkyCoord(ra=header['RA']*u.deg, dec=header['DEC']*u.deg)
        VBAR = sc.radial_velocity_correction(obstime=Time(DATE), location=loc)

        # Add some keywords with comments to the header
        header['I-TELES'] = (TELES, TELESINFO)
        header['I-INSTR'] = (INSTR, 'Instrument')
        header['I-RESOL'] = (int(RESOL[1:]), 'Nominal resolving power (R)')
        header['I-DATE'] = (DATE[:21], '[YYYY-MM-DDTHH:MM:SS.S] Adq.date at midpoint')
        header['I-VBAR'] = (round(VBAR.to(u.km/u.s).value, 3), 'Barycentric velocity correction (km/s)')
        header['I-HJD'] = (round(Time(DATE, scale='utc').jd, 5), 'Heliocentric JD at midpoint')
        header['I-RA'] = (round(header['RA'], 5), 'RA [deg]')
        header['I-DEC'] = (round(header['DEC'], 5), 'DEC [deg]')
        header['I-TEXP'] = (round(header['EXPTIME'], 1), '[s]')

        # Get the wavelength of the spectrum
        ctype = header['CTYPE1']       # Type of wavelength calibration (linear, log, etc.)
        lam0 = header['CRVAL1']        # Get the wavelength of the first pixel
        dlam = header['CDELT1']        # Step of increase in wavelength
        pix0 = header['CRPIX1']        # Reference pixel (generally 1, FEROS -49)
        spec_length = header['NAXIS1'] # Length of the spectrum
        # Calculate the wavelength and flux arrays
        flux = hdu0.data
        wave = lam0 + dlam*(np.arange(spec_length) - pix0 + 1)
        if 'log' in ctype:
            wave = np.exp(wave)

        # Normalize the spectrum
        ft, ynt, fkt, lam, snr, snr4500 = normalize(TELES, wave, flux)
        # Add the normalize flux to the fits data
        hdu0.data = [ft, hdu0.data]

        if plot:
            fig, ax = plt.subplots(figsize=(10, 5))
            ax.plot(wave, flux/np.nanmean(flux), label='Original', lw=0.5, alpha=0.5)
            ax.plot(wave, ft, label='Normalized', lw=0.5, alpha=0.5)
            ax.plot(wave, ynt/np.nanmean(ynt), label='Continuum fit', lw=0.5, alpha=0.5)
            ax.plot(wave, fkt, label='Mask', lw=0.5, alpha=0.5, drawstyle='steps-mid')
            ax.set_title(f"{object}_{DATE}_{telcode}_{RESOL}")
            ax.legend()
            fig.tight_layout()
            plt.show(block=False)

        header['I-SNR'] = (int(np.median(snr)), 'Median SNR of the normalized regions')
        header['I-SN450'] = (int(snr4500), 'SNR (4500 A)')

        # Add the spectral classification to the header
        if table_with_spc is not None:
            spc, spc_ref = retrieve_spc(object, table_with_spc)

        header['I-SPC'] = (spc, 'Spectral classification')
        header['I-SPCREF'] = (spc_ref, 'Spectral classification reference')

        # Add the data release version
        header['I-VERS'] = ('DR3.0', 'Data release version')

        # Add the comments
        header['COMMENT'] = ' --- The IACOB spectroscopic database of Galactic OB stars'
        if telcode == '_N_':
            header['COMMENT'] = ' --- FEROS spectrum downloaded from the ESO archive'
        header['COMMENT'] = ' --- Some extra information about the spectra has been added above'
        if telcode == '_N_' or (telcode == '_M_' and not '_log' in spectrum):
            header['COMMENT'] = ' --- The spectrum is NOT corrected from barycentric velocity'
            header['COMMENT'] = ' --- Use the keyword I-VBAR to correct from barycentric velocity'
        elif telcode == '_F_' or (telcode == '_M_' and '_log' in spectrum):
            header['COMMENT'] = ' --- The spectrum is corrected from barycentric velocity'
        header['COMMENT'] = ' --- A normalized and original merged spectrum are included'
        header['COMMENT'] = ' --- Be aware of the reference for the I-SPC keyword'
        header.add_blank('', before='COMMENT')

        end = '.fits' if not '_log' in ctype else '_log.fits'
        new_name = object + '_' + date + '_' + time + telcode + RESOL + end

        # save the fits file with the new name in the same folder
        print('Renaming',spectrum,'to',new_name)
        hdu.writeto(os.path.join(path_to_spectra, new_name), overwrite=True)
        hdu.close()

    return None


def calculate_snr(flux):
    # Calculate the SNR of the spectrum with an iterative sigma clipping method
    me0 = np.nanmean(flux)
    st0 = np.nanstd(flux)
    snr0 = me0 / st0
    snr = snr0 + 100
    while abs(snr0 - snr) >= 1.0:
        snr = snr0
        tt = np.where(np.abs(flux - me0) <= 2.0 * st0)
        me0 = np.nanmean(flux[tt])
        st0 = np.nanstd(flux[tt])
        snr0 = me0 / st0

    return int(snr0)


def normalize(telescope, w, f):
    '''
    Normalize a complete spectrum by regions.
    Parameters
    ----------
    telescope : str
        Name of the telescope (e.g., 'NOT', 'MERCATOR', 'MPI-2.2').
    w : array-like
        Wavelength array of the spectrum.
    f : array-like
        Flux array of the spectrum.
    Returns
    -------
    ft : array-like
        Normalized flux array of the spectrum.
    ynt : array-like
        Continuum fit of the spectrum.
    fkt : array-like
        Mask of the points used for the normalization (1 for good points, 0 for bad points).
    lam : array-like
        Wavelengths of the strong lines used for normalization.
    '''

    lam0 = np.min(w.astype(float))

    # Selection of the individual regions for normalization
    if telescope == 'MERCATOR':
        #! TO BE DEFINED
        lam = np.array([3815,3920,4060,4150,4290,4605,4830,4950,5070,5320,5500,5620,5750,5840,
                        5980,6200,6460,6700,7100,7620,8000,8400,8710,8930,9020,9300,9400], float)
    elif telescope == 'NOT':
        lam = np.array([3815,3920,4060,4150,4290,4605,4830,4950,5070,5320,5500,5620,5750,5840,
                        5980,6200,6460,6700,7100,7620,8000,8400,8710,8930], float)
    elif telescope == 'MPI-2.2':
        #! TO BE DEFINED
        lam = np.array([3815,3920,4060,4150,4290,4605,4830,4950,5070,5320,5500,5620,5750,5840,
                        5980,6200,6460,6700,7100,7620,8000,8400,8710,8930,9020,9300,9400], float)
    else:
        print(f"Warning: Telescope '{telescope}' not recognized. Returning original spectrum.")
        return f, np.zeros_like(w), np.ones_like(w, float), lam, np.array([np.nan]), np.nan

    # Take regions within the range of the spectrum
    mask = (lam >= w.min()) & (lam <= w.max())
    lam = np.concatenate([lam[mask], [w.max()]]).astype(float)

    ws_list, fs_list, fn_list, yn_list, fk_list, snr_list = [], [], [], [], [], []
    # First region
    ws0, fs0, fn0, yn0, fk0, snr0, wlim1, wlim2, lam1 = normalize_slice(w, f, lam0, lam[0])
    ws0, fs0, fn0, yn0, fk0, snr0, wlim1, wlim2, lam1 = normalize_slice(
        w, f, lam0, lam[0], wlim1=wlim1, wlim2=wlim2, lam1=lam1, iter=1)
    ws_list.append(ws0); fs_list.append(fs0); fn_list.append(fn0); yn_list.append(yn0); fk_list.append(fk0)
    snr_list.append(snr0)

    # Rest of the regions
    for i in range(1, len(lam)):
        ws0, fs0, fn0, yn0, fk0, snr0, wlim1, wlim2, lam1 = normalize_slice(
            w, f, lam[i-1], lam[i], wlim1=wlim1, wlim2=wlim2, lam1=lam1)
        ws0, fs0, fn0, yn0, fk0, snr0, wlim1, wlim2, lam1 = normalize_slice(
            w, f, lam[i-1], lam[i], wlim1=wlim1, wlim2=wlim2, lam1=lam1, iter=1)
        ws_list.append(ws0); fs_list.append(fs0); fn_list.append(fn0); yn_list.append(yn0); fk_list.append(fk0)
        snr_list.append(snr0)

    ws = np.concatenate(ws_list)
    fs = np.concatenate(fs_list)
    fn = np.concatenate(fn_list)
    yn = np.concatenate(yn_list)
    fk = np.concatenate(fk_list)

    snr = np.array(snr_list)
    # find the SNR for the index of lam where lam is closest to 4500
    idx_4500 = np.argmin(np.abs(lam - 4500))
    snr4500 = snr[idx_4500] if idx_4500 < len(snr) else np.nan

    if ws is None or fn is None or len(ws) == 0 or len(fn) == 0:
        ft = np.zeros_like(w)
        ynt = np.zeros_like(w)
    else:
        ft = np.interp(w, ws, fn)
        ynt = np.interp(w, ws, yn)
    # global continuum interpolated
    # global mask passing from ws to w, trusting ws with w
    fkt = np.zeros_like(w, float)
    fkt[np.searchsorted(w, ws)] = fk

    return ft, ynt, fkt, lam, snr, snr4500


def normalize_slice(w, f, w0, w1, wlim1=None, wlim2=None, lam1=None, iter=None):

    '''
    Function to normalize the spectrum by fitting a polynomial to the continuum
    regions. It removes the strong lines from the spectrum and iteratively fits a
    polynomial to the remaining points, applying sigma-clipping to remove points
    affected by weak lines or noise.

    Parameters
    ----------
    w : array-like
        Wavelength array of the spectrum.
    f : array-like
        Flux array of the spectrum.
    w0 : float
        Lower limit of the wavelength region to be normalized.
    w1 : float
        Upper limit of the wavelength region to be normalized.
    wlim1 : array-like, optional
        Lower limits of the wavelength regions to be normalized.
    wlim2 : array-like, optional
        Upper limits of the wavelength regions to be normalized.
    lam1 : array-like, optional
        Wavelengths of the strong lines to be removed.
    iter : int, optional
        If not None, it indicates that this is an iterative call to the function, and the strong lines should not be added to lam1 again.

    Returns
    -------
    ws0 : array-like
        Wavelength array of the normalized spectrum in the given region.
    fs0 : array-like
        Flux array of the normalized spectrum in the given region.
    fn : array-like
        Normalized flux array of the spectrum in the given region.
    yn : array-like
        Continuum fit of the spectrum in the given region.
    fk : array-like
        Mask of the points used for the normalization (1 for good points, 0 for bad points) in the given region.
    wlim1 : array-like
        Updated lower limits of the wavelength regions to be normalized.
    wlim2 : array-like
        Updated upper limits of the wavelength regions to be normalized.
    lam1 : array-like
        Updated wavelengths of the strong lines to be removed.
    '''

    # Strong lines
    lamH = np.array([3712,3722,3735,3750,3771,3797,3835,3890,3970,
                     4102,4340,4860,6563,8438,8465,8500,8545,8600,8665,8748,8860], dtype=float)
    lamHeI  = np.array([3820,3926,4009.,4026.,4143.,4387.,4471.5,4713,4922,5875,6678])
    lamHeII = np.array([4200,4542,4686,5411])
    lamISM  = np.array([4430])

    lam0  = np.concatenate([lamH, lamHeI, lamHeII, lamISM])
    dlam0 = np.concatenate([lamH*0+5., lamHeI*0+2., lamHeII*0+2., lamISM*0+5.])
    tt = np.argsort(lam0)
    lam0 = lam0[tt]
    dlam0 = dlam0[tt]
    snr = np.nan

    # Cut region
    mask = (w > w0) & (w <= w1) & (f > 0)
    if np.sum(mask) < 50:
        print(f"Warning: Region has too few points to be normalized.")
        # Return neutral arrays
        ws0 = w[mask] if mask.any() else np.array([])
        fs0 = f[mask] if mask.any() else np.array([])
        return ws0, fs0, np.ones_like(ws0), np.zeros_like(ws0), np.zeros_like(ws0), snr, wlim1, wlim2, lam1

    ws0 = w[mask]
    fs0 = f[mask]
    ws00 = ws0.copy()
    fs00 = fs0.copy()

    # Keywords initialization
    lam1 = [] if lam1 is None else list(lam1)
    if wlim1 is None or not isinstance(wlim1, np.ndarray) or len(wlim1) != len(lam0):
        wlim1 = np.zeros(len(lam0))
    if wlim2 is None or not isinstance(wlim2, np.ndarray) or len(wlim2) != len(lam0):
        wlim2 = np.zeros(len(lam0))

    # Calculates the FWHM of the line as average size of the line wings
    # The * 0.9 factor reduces that width by 10%.
    if iter is not None:
        dlam0 = 0.5*(np.abs(wlim2)+np.abs(wlim1))*0.9

    # Remove strong lines
    for j in range(len(lam0)):
        if ws0.size == 0:
            break
        if (ws0.min()-10 <= lam0[j]) and (ws0.max()+10 >= lam0[j]):
            good = (ws0 <= lam0[j]-dlam0[j]) | (ws0 >= lam0[j]+dlam0[j])
            ws0 = ws0[good]
            fs0 = fs0[good]
            if iter is None:
                lam1.append(lam0[j])

    ws = ws0
    fs = fs0

    # Initial adjustment of the continuum. If there are not enough points, return the original spectrum.
    ord = 1
    if len(ws) < ord+1:
        # Not enough points to fit a polynomial of the given order.
        print(f"Warning: not enough points to fit a polynomial of order {ord}")
        fn = np.ones_like(ws00)
        yn = np.zeros_like(ws00)
        fk = np.zeros_like(ws00)
        if np.shape(fn)!=np.shape(yn):
            print('fn return de normaliza region sin puntos: ', np.shape(fn))
            print('yn return de normaliza region sin puntos: ', np.shape(yn))
        return ws00, fs00, fn, yn, fk, snr, wlim1, wlim2, np.array(lam1)

    # Note: Simón-Díaz used to add an extra smoothing here (not implemented)
    dat = np.polyfit(ws, fs, ord)
    ys  = np.polyval(dat, ws)

    fn = fs / ys
    fk = np.ones_like(fn)
    mfn = np.mean(fn)
    sfn = np.std(fn)

    # Iterative sigma-clipping
    error0, error = 40., 20.
    value, step = 2., 0
    while abs(error-error0) >= 0.05 and ws.size > 0:
        for i in range(len(ws)):
            if fn[i] < mfn-sfn or fn[i] > mfn+value*sfn:
                bad = (fn < mfn - sfn) | (fn > mfn + value * sfn)
                fn[bad] = mfn + np.random.randn(np.sum(bad)) * sfn
                fk[bad] = 0.

        nn = np.where(fk > 0)[0]
        error0 = error
        if len(nn) < 2:
            break
        dat2 = np.polyfit(ws[nn], fn[nn], 2)
        ys2  = np.polyval(dat2, ws)
        error2 = np.std(fn - ys2)

        dat1 = np.polyfit(ws[nn], fn[nn], 1)
        ys1  = np.polyval(dat1, ws)
        error1 = np.std(fn - ys1)

        if error2 <= error1:
            ys = ys2
            error = abs(error2)
        else:
            ys = ys1
            error = abs(error1)

        if w0 >= 8300:
            dat = np.polyfit(ws[nn], fn[nn], 1)
            ys  = np.polyval(dat, ws)
            error = np.std(fn - ys)

        fn = fn / ys
        mfn = np.median(fn)
        sfn = np.std(fn)
        step += 1

    # Final normalization
    nn = np.where(fk > 0)[0]
    wp = ws[nn] if nn.size > 0 else np.array([])
    fp = fs[nn] if nn.size > 0 else np.array([])

    ord = 2 if w0 < 8300 else 1
    for _ in range(3):
        if wp.size < ord+1:
            break
        dat = np.polyfit(wp, fp, ord)
        yfit = np.polyval(dat, wp) / fp
        keep = np.abs(yfit - np.median(yfit)) <= 2*np.std(yfit)
        wp = wp[keep]
        fp = fp[keep]

    dat = np.polyfit(wp, fp, ord) if wp.size >= ord+1 else np.array([1.])
    yn     = np.polyval(dat, ws00) if wp.size >= ord+1 else np.ones_like(ws00)
    ynoise = np.polyval(dat, ws0)   if wp.size >= ord+1 else np.ones_like(ws0)

    fn = fs00 / yn # This is the key to why fk has fewer points. fn is redefined here, but fk is not.

    # Calculate the SNR - NOT ENABLED
    fnoise = fs0 / ynoise
    snr = calculate_snr(fnoise) if ws0.size > 0 else 0

    # Save the limits of the§ strong lines for the next iteration.
    if iter is None and len(lam1) > 0:
        lam1 = np.array(lam1)
        for line in lam1:
            tt = np.where(lam0 == line)[0]
            if wp.size > 0:
                left = np.where(wp <= line)[0]
                wlim1[tt] = wp[left.max()] - line if left.size > 0 else wp.min() - line
                right = np.where(wp >= line)[0]
                wlim2[tt] = wp[right.min()] - line if right.size > 0 else wp.max() - line
            else:
                wlim1[tt] = 0
                wlim2[tt] = 0

    # Fix small mismatch as fk was defined in ws, and everything else in ws00. It was coming out smaller
    fkt = np.zeros_like(ws00)       #fkt = np.zeros_like(ws00, dtype=float)
    idx = np.searchsorted(ws00, ws) #idx = np.searchsorted(ws00, ws)
    # ensure that we only assign valid indices: valid = (idx >= 0) & (idx < len(ws00)) [possible alternative]
    fkt[idx] = fk                   #fkt[idx[valid]] = fk[valid]

    return ws00, fs00, fn, yn, fkt, snr, wlim1, wlim2, lam1


def retrieve_spc(star_id, spc_table):
    """
    Function to ingest the spectral classification of the stars into the FITS headers.

    The user must provide a table containing the names of the stars in the first column,
    the spectral classification in the second column, and reference in the third column.

    Parameters
    ----------
    star_id : str
        Name of the star to update the spectral classification for.

    spc_table : pandas.DataFrame
        Path to the table containing the spectral classification information for the stars.
        The table must be a txt or csv file with the following format:

            ID, spc, spc_ref
            HD12345, B8Ia, Negueruela et al. (2018)

    Returns
    -------
    spc : str
        Spectral classification of the star.
    spc_ref : str
        Reference for the spectral classification of the star.
    """

    # Open the spectral classification table and create a dictionary with the information
    table = {}
    for line in open(spc_table, 'r').readlines():
        id, spc, spc_ref = line.split(',')
        table[id.upper()] = (spc, spc_ref)

    if star_id.upper() not in table:
        print(f"Warning: Spectral classification not found for star ID '{star_id}'")
        return 'Unknown', 'Unknown'
    else:
        spc = table[star_id.upper()][0].strip()
        spc_ref = table[star_id.upper()][1].strip()

    return spc, spc_ref


def degrade_spec(path_to_spectra, output_dir, resol=5000, lwl=None, rwl=None, step=None,
                 delimiter=';', extension='csv', output_name='fullname'):
    """
    Function to degrade the spectra in a given directory to a given resolution and/or delta-lambda.

    Parameters
    ----------
    path_to_spectra : str
        Path to the directory containing the spectra in FITS format.

    output_dir : str
        Path to the directory where the degraded spectra will be stored.

    resol : int, optional
        Resolution to which the spectra will be degraded. Default is 5000.

    lwl : float, optional
        Left wavelength limit of the spectra. Default is None (no change).
        Only applied if dlam is not None.

    rwl : float, optional
        Right wavelength limit of the spectra. Default is None (no change).
        Only applied if dlam is not None.

    step : float, optional
        The step size or delta-lambda to which the spectra will be spaced.
        Default is None (no change).

    delimiter : str, optional
        Delimiter of the output files. Default is ;.

    extension : str, optional
        Extension of the output files. Only 'ascii','txt','csv','fits' are accepted.
        Default is csv.

    output_name : str, optional
        Name of the output files. Default is fullname (the full name of the input file).
        Alternatively, the user can choose 'name' (only the name of the input file).
        Default is fullname.

    Returns
    -------
    None, but the degraded and/or resampled input spectra are stored in the output_dir.
    """

    spectra = [os.path.join(path_to_spectra, i) for i in os.listdir(path_to_spectra) if i.endswith('.fits')]

    for spectrum in spectra:

        # Retrieve the key values fron the fits header
        hdu = fits.open(spectrum)       # Open the fits image file
        hdu.verify('fix')               # Fix possible issues with the keywords
        header0 = hdu[0].header         # Read header of primary extension

        instrum = header0['INSTRUME']   # Instrument

        lam0 = header0['CRVAL1']        # Get the wavelenght of the first pixel
        dlam = header0['CDELT1']        # Step of increase in wavelength
        pix0 = header0['CRPIX1']        # Reference pixel (generally 1, FEROS -49)
        spec_length = header0['NAXIS1'] # Length of the spectrum

        if 'I-VBAR' in header0:
            vbar = header0['I-VBAR']    # Barycentric velocity correction
        else:
            print('No barycentric correction applied to' + spectrum.split(os.sep)[-1])
            vbar = 0

        # Correct Mercator CRVAL1 20101018-19:
        if any(bad in spectrum for bad in ['_20101018_','_20101019_']) and lam0 == 3763.9375:
            lam0 = 3763.61

        # Make lists with wavelenght and flux for each spectrum
        # Those with log or from FEROS are already corrected from heliocentric velocity
        wave = lam0 + dlam*(np.arange(spec_length) - pix0 + 1)
        if '_log' in spectrum:
            wave = np.exp(wave)
        elif not instrum == 'FEROS':
            wave = wave*(1 + 1000*vbar/cte.c)

        # Get the flux of the spectrum
        try:
            flux = hdu[0].data[0] # IACOB
        except:
            flux = hdu[0].data # Other

        # Close the fits image file
        hdu.close()

        # Cut the spectrum to the desired wavelenght range given by lwl and rwl
        if lwl != None and rwl != None:
            if wave[0] > lwl+dlam or wave[-1] < rwl-dlam:
                print('WARNING: Wavelength limits outside spectrum wavelength range.')
            flux = flux[(wave >= lwl) & (wave <= rwl)]
            wave = wave[(wave >= lwl) & (wave <= rwl)]

        # Recalculate the delta lambda if the spectrum is in log
        if '_log' in spectrum.split(os.sep)[-1]:
            dlam = (wave[-1]-wave[0])/(len(wave)-1)

        # Degradation of the spectum
        lambda0 = np.mean(wave)

        sigma = lambda0/(2.35482*float(resol))

        x = np.arange(-10*sigma, 10*sigma+dlam, dlam)
        gauss = f_gaussian(x, sigma)
        kernel = gauss/np.trapz(gauss)

        # Remove the nans from the flux
        mask = np.where(np.isnan(flux) == False)[0]

        # Convolve the flux with the kernel keeping the nans of the original flux
        flux[mask] = 1 + convolve(flux[mask] - 1, kernel, mode='same')

        # Resample the spectrum to the new delta lambda
        if step is not None:

            if not isinstance(step, float) and not isinstance(step, int):
                print('ERROR: The input delta lambda is not a float or an integrer.')
                return None

            if step > np.mean(wave)/resol/3:
                # It is divided by 3 to at least have 3 pixels in a gaussian
                print('WARNING: The new delta lambda implies lossing information...')

            if lwl is None or lwl < wave[0]:
                lwl = wave[0]
            if rwl is None or rwl > wave[-1]:
                rwl = wave[-1]

            # Interpolate the spectrum to the new delta lambda / step size
            f = interp1d(wave, flux, kind='linear', fill_value='extrapolate')
            wave = np.arange(lwl, rwl, step)
            flux = f(wave)

        # Remove the nan values from the flux and waveleght arrays
        mask = np.where(np.isnan(flux) == False)[0]
        flux = flux[mask]
        wave = wave[mask]

        # Save the spectrum
        if output_name == 'fullname':
            filename = spectrum.split(os.sep)[-1].replace('.fits', '')
        elif output_name == 'name':
            filename = spectrum.split(os.sep)[-1].split('_')[0]

        # If directory "degraded" inside output_dir does not exist, it is created
        if not os.path.exists(output_dir+'degraded/'):
            os.makedirs(output_dir+'degraded/')

        # Save the spectrum
        if extension in ['ascii', 'txt', 'csv']:
            np.savetxt(output_dir+'degraded/%s' % (filename + '.' + extension),
                np.transpose([wave, flux]), delimiter=delimiter, fmt=('%.4f','%.6f'),
                header='lambda%sflux' % delimiter, comments='')
        elif extension == 'fits':
            hdu = fits.PrimaryHDU(np.transpose([flux, wave]))
            hdu.writeto(output_dir+'degraded/%s' % (filename + '.' + extension), overwrite=True)

    return None


def fix_edges(path_to_spectra):
    """
    Function to correct the edges of the normalized spectra of the files in the given directory.

    Parameters
    ----------
    path_to_spectra : str
        Path to the directory containing the spectra in FITS format.

    Returns
    -------
    None
    """

    spectra = [os.path.join(path_to_spectra, i) for i in os.listdir(path_to_spectra) if i.endswith('.fits')]

    for spectrum in spectra:

        hdu = fits.open(spectrum) # Open the fits image file
        hdu.verify('fix') # Fix possible issues with the keywords

        flux = hdu[0].data[0]
        wave = hdu[0].data[1]

        if flux[1] < flux[0]:
            l_cut = next(pix for pix in range(len(wave)) if flux[pix+1] > flux[pix])
        else:
            l_cut = next(pix for pix in range(len(wave)) if flux[pix+1] < flux[pix])

        if flux[-2] < flux[-1]:
            r_cut = next(pix for pix in reversed(range(len(wave))) if flux[pix-1] > flux[pix])
        else:
            r_cut = next(pix for pix in reversed(range(len(wave))) if flux[pix-1] < flux[pix])

        flux[:l_cut] = np.nan
        flux[r_cut:] = np.nan

        hdu[0].data[0] = flux

        hdu.writeto(spectrum, output_verify='ignore', overwrite=True)

    return None


def f_gaussian(x, sigma):
    return np.exp(-(x/sigma)**2/2)