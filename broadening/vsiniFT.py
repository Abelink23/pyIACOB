from __future__ import annotations

# =======================================================================
# ============================== VSINI(FT) ==============================
# =======================================================================

# ==============================
def _vsiniFT(
             lines     : list[float],
             w         : np.ndarray[float],
             f         : np.ndarray[float],
             res       : float,
             xr0       : np.ndarray[float] | None,
             xl        : np.ndarray[float] | None = None,
             xyc       : np.ndarray[float] | None = None,
             xsnr      : np.ndarray[float] | None = None,
             xrec      : np.ndarray[float] | None = None,
             rerec     : bool = False,
             vmin      : float | None = None,
             beta      : float | None = None,
             sigma_FT0 : float | None = None,
             w_intprof : np.ndarray[float] | None = None,
             f_intprof : np.ndarray[float] | None = None
    ) -> dict:
    """
    Pre-process a spectral line and measure the projected rotational
    velocity from the first zero of its Fourier transform.

    The observed line is prepared for the Fourier transform (FT) analysis:
    the spectral window, the local continuum (optional renormalisation), the
    signal-to-noise ratio and the limits of the line are defined, the line is
    interpolated to a fixed wavelength step and normalised to the continuum
    level, and part of the line can be clipped. The FT of the line is then
    shown and the user marks its first zero, which gives vsini(FT).

    Every quantity given as a keyword (xl, xyc, xsnr, xrec) is used
    directly; otherwise it is selected interactively on the plots.

    ----------
    Parameters
    ----------
        w : np.ndarray[float]
            Wavelength values of the spectrum [A].

        f : np.ndarray[float]
            Normalised flux values of the spectrum.

        res : float
            Resolving power of the spectrum. Only used to broaden the
            intrinsic profile of the user before its FT.

        xr0 : np.ndarray[float] | None
            Initial wavelength window around the line [A]. If None, it is
            selected interactively with two mouse clicks.

        lines : list[float]
            Wavelengths of the lines of pyiacob_broad_lines.dat [A], drawn as
            vertical dashed lines in the interactive plots.

        xl : np.ndarray[float] | None, optional
            Limits of the line used for the FT and GOF [A]. If None, they
            are selected interactively (and clipping is offered).

        xyc : np.ndarray[float] | None, optional
            Wavelength of the line core [A] and flux of the adjacent
            continuum. If None, selected with one click.

        xsnr : np.ndarray[float] | None, optional
            Continuum window used to compute the SNR [A]. If None, selected
            with two clicks.

        xrec : np.ndarray[float] | None, optional
            Two continuum windows [x1, x2, x3, x4] for the local
            renormalisation [A] (only with rerec). If None, selected with
            four clicks.

        rerec : bool, optional
            Local renormalisation of the continuum with a straight line
            fitted to the xrec windows. Defaults to False.

        vr : float, optional
            Radial velocity correction applied to the wavelengths [km/s].
            Defaults to 0.

        vmin : float | None, optional
            Minimum vsini and vmac of the GOF grid [km/s]. If None, one pixel
            of the line in velocity units.

        beta : float | None, optional
            Limb-darkening coefficient, beta = epsilon/(1 - epsilon). Only
            used to broaden the intrinsic profile of the user.

        sigma_FT0 : float | None, optional
            Position of the first zero of the rotational FT, in units of
            lambda*vsini/c. If None, 0.660.

        w_intprof : np.ndarray[float] | None, optional
            Wavelength [A] of the intrinsic profile of the user.

        f_intprof : np.ndarray[float] | None, optional
            Normalised flux of the intrinsic profile of the user.

    -------
    Returns
    -------
        prep : dict
            Pre-processed line (saved in <star>_<line>.idl):

            - 'vmin' : float
                Minimum vsini and vmac of the GOF grid [km/s].

            - 'w0' : np.ndarray[float]
                Wavelengths of the spectrum in the initial window,
                normalised to the continuum level.

            - 'f0' : np.ndarray[float]
                Fluxes of the spectrum in the initial window, normalised to
                the continuum level.

            - 'wx0' : np.ndarray[float]
                Wavelengths of the line profile used for the FT and the GOF,
                on a fixed wavelength step (clipped points removed).

            - 'fx0' : np.ndarray[float]
                Fluxes of the line profile used for the FT and the GOF, on a
                fixed wavelength step (clipped points removed).

            - 'xc' : float
                Wavelength of the line core [A].

            - 'ew' : float
                Equivalent width of the line from its FT [mA].

            - 'snr' : float
                Signal-to-noise ratio of the adjacent continuum.

            - 'vft' : float
                vsini from the first zero of the FT [km/s].

            - 'xfou_obs' : np.ndarray[float]
                FT of the line: vsini scale [km/s].

            - 'yfou_obs' : np.ndarray[float]
                FT of the line: log10 of the normalised amplitude.

            - 'max_sigma_obs' : float
                Minimum measurable vsini set by the Nyquist frequency
                [km/s].

            - 'clip' : int
                1 if part of the line was clipped, 0 otherwise.

            - 'w_intprof' : Only with an intrinsic profile of the user:
                wavelengths of the profile moved to the line core and
                scaled to the observed EW.

            - 'f_intprof' : Only with an intrinsic profile of the user:
                fluxes of the profile moved to the line core and scaled to
                the observed EW.

            - 'ew_intprof' : Only with an intrinsic profile of the user:
                its original EW [mA].

    -----
    Notes
    -----
        - The line is interpolated to a fixed wavelength step equal to the IDL
          MEDIAN of the pixel differences, as required by the FT.

        - The fluxes are normalised to the continuum level ycont marked by the
          user: F' = F - ycont + 1.

        - If a region is clipped, its flux is replaced by a straight line
          between the closest points outside it for the FT, and the points are
          removed afterwards for the GOF.

        - With an intrinsic profile of the user, its FT (broadened with
          vsini = 5 km/s and the instrument) is also shown, and, as in IDL,
          max_sigma_obs is the Nyquist value of that FT.
    """

    import numpy as np

    # ########## PANORAMA ##########
    from .fourier import vfft
    from .interactive import _COLORS, _DASH_DOT3, _clip_curves, _clip_line, _inside, _limits, _show, _window, _xhair
    from .numerics import _interpol, _median_upper
    from .profiles import _convol
    from .variables import C_KMS


    # ########## LINE SELECTION ##########
    """
    - When no line identifier is provided, the user is prompted to select
      the region of the line interactively.
    """

    # ##### INTERACTIVE FIGURE #####
    fig = _window()

    # ##### SELECTING THE SPECTRAL WINDOW #####
    if xr0 is None:
        if xl is not None:
            raise ValueError("A 'line' is needed to define the initial window when 'xl' is given.")

        ax = _show(fig, 'Select the region where the line is located', [(w, f, {})], lines)
        xr0 = _limits(None, ax) # Two clicks around the line

    # ##### SPECTRAL WINDOW #####
    m = _inside(w, xr0)
    w0, f0 = w[m], f[m] # Spectrum in the initial window

    y_sh = min(f0.max(), 1.5) - f0.min() # Flux range of the plots (top limited to 1.5)
    y_range = [f0.min() - 0.1 * y_sh, f0.max() + 0.1 * y_sh]


    # #################### PREPROCESS ####################
    # ########## NORMALIZATION ##########
    """
    - Optional local renormalisation: a straight line is fitted to the
      points of the two continuum windows of xrec and the flux is divided
      by it.
    - While the windows are selected, the points of each window are drawn
      in red. In the next plots, the original spectrum is drawn in grey
      and the renormalised one in black.
    """
    raw = [] # Original spectrum (grey) drawn below the renormalised one in the next plots
    if rerec:
        if xrec is None:
            ax = _show(fig,
                       'Normalization: Select two continuum regions',
                       [(w0, f0, {})],
                       lines,
                       xr0,
                       y_range
                      )
            xrec = []
            for _ in range(2): # Two clicks per continuum window; its points are drawn in red
                xw = [x for x, _ in _xhair(ax, 2)]
                mw = _inside(w0, np.sort(xw))
                ax.plot(w0[mw], f0[mw], color=_COLORS[2], lw=0.6)
                fig.canvas.draw_idle()
                xrec += xw
            print(f"{'xrec':>13} = [{', '.join(f'{x:.3f}' for x in xrec)}] A")
        m = _inside(w0, np.sort(xrec[:2])) | _inside(w0, np.sort(xrec[2:4])) # Points of both continuum windows
        coef = np.polynomial.polynomial.polyfit(w0[m], f0[m], 1) # Straight line fitted to the continuum
        raw = [(w0, f0, dict(color=_COLORS[11]))]
        f0 = f0 / np.polynomial.polynomial.polyval(w0, coef) # Renormalised flux
        y_range = [0.9 * f0.min(), min(1.1 * f0.max(), 1.5)]


    # ########## SNR COMPUTATION ##########
    """
    - SNR = 1 / standard deviation of the flux in the continuum window xsnr.
    """
    ax = None if xsnr is not None else _show(
                                             fig,
                                             'Mark a region to compute the SNR',
                                             raw + [(w0, f0, {})],
                                             lines,
                                             xr0,
                                             y_range
                                            )
    xs = _limits(xsnr, ax) # Limits of the continuum window (keyword or two clicks)
    msnr = _inside(w0, xs)
    snr = 1. / np.std(f0[msnr], ddof=1) # Sample standard deviation (N-1), as IDL STDDEV

    print(f"{'xsnr':>13} = [{xs[0]:.3f}, {xs[1]:.3f}] A")
    print(f"{'snr':>13} = {int(snr)}\n")


    # ########## LINE LIMIT ##########
    ax = None if xl is not None else _show(
                                           fig,
                                           'Mark the limits of the line',
                                           raw + [(w0, f0, {}),
                                            (w0[msnr], f0[msnr], dict(color=_COLORS[2]))],
                                           lines,
                                           xr0,
                                           y_range)
    x1, x2 = xlim = _limits(xl, ax) # Limits of the line (keyword or two clicks)
    m = _inside(w0, xlim)
    wx0, fx0 = w0[m], f0[m] # Line profile used for the FT and the GOF


    # ########## INTERPOLATION TO A FIXED WAVELENGTH STEP ##########
    """
    - The FT needs an equidistant grid. The step is the IDL MEDIAN of
      |wx0 - SHIFT(wx0, 1)| (the first difference wraps around), and only
      the points within the original range are kept.
    """
    w_sh = _median_upper(np.abs(wx0 - np.roll(wx0, 1))) # Fixed wavelength step
    wxa = wx0.min() + np.arange(int((wx0.max() - wx0.min()) / w_sh)) * w_sh # Equidistant grid
    keep = _inside(wxa, [wx0.min(), wx0.max()]) # Points within the original range
    wx0, fx0 = wxa[keep], _interpol(fx0, wx0, wxa)[keep] # Line on the equidistant grid

    span = wx0.max() - wx0.min() # Plot ranges
    x_range = [wx0.min() - 0.3 * span, wx0.max() + 0.3 * span]
    y_range = [0.98 * fx0.min(), min(1.02 * fx0.max(), 1.5)]


    # ########## LINE CORE & CONTINUUM ##########
    # ##### INTERACTIVE SELECTION #####
    if xyc is None:
        ax = _show(fig,
                   'Mark the continuum & the core of the line (just 1 click!)',
                   raw + [(w0, f0, {}), (wx0, fx0, dict(color=_COLORS[2]))],
                   lines,
                   xr0,
                   y_range)
        xyc = _xhair(ax)[0] # One click: x = line core, y = continuum level
        ax.axvline(xyc[0], color=_COLORS[3], ls='--') # Mark the selected line core
        ax.axhline(xyc[1], color=_COLORS[3], ls='--') # Mark the selected continuum level
        fig.canvas.draw_idle()

    # ##### LINE CORE & CONTINUUM #####
    xc, ycont = xyc[0], xyc[1] # Wavelength of the line core [A] and continuum level

    # ##### VSINI RANGE #####
    vft_max = np.mean(np.abs(np.array([x1, x2]) / xc - 1)) * C_KMS # Maximum vsini of the FT plot: half-width of the window in km/s
    print(f"{'vft_max':>13} = {vft_max:.1f} km/s")

    fx0, f0 = fx0 - ycont + 1., f0 - ycont + 1. # Normalisation to the continuum level: F' = F - ycont + 1
    print(f"{'xl':>13} = [{x1:.3f}, {x2:.3f}] A")
    print(f"{'xyc':>13} = [{xc:.3f}, {ycont:.4f}]")

    vmin = vmin or _median_upper(np.abs(wx0 - np.roll(wx0, 1))) * C_KMS / xc # Minimum vsini and vmac: one pixel in km/s (unless given by the user)
    print(f"{'vmin':>13} = {vmin:.2f} km/s")

    # ##### INTRINSIC PROFILE #####
    if w_intprof is not None:
        w_intprof = w_intprof - w_intprof[np.argmin(f_intprof)] + xc # Minimum of the user profile moved to the line core


    # ########## CLIPPING ##########
    """
    - Only in interactive mode (xl not given), as in IDL.
    """
    cut = None # Mask of clipped points (None -> no clipping)
    if xl is None:
        fx0, cut = _clip_line(fig, w0, f0, wx0, fx0, lines, xr0, xc, 1.) # fx0 is already normalised to continuum = 1


    # #################### FOURIER TRANSFORM ####################
    # ########## PLOTTING ##########
    fig = _window() # Left panel: line profile; right panel: its FT
    ax_line = _show(fig,
                    '',
                    _clip_curves(w0, f0, wx0, fx0, cut, shift=ycont - 1.), # Clipped regions: interpolation in blue, real data hidden
                    xlim = x_range,
                    ylim = y_range,
                    ax = fig.add_subplot(121))
    ax_line.axvline(xc, ls='--', color=_COLORS[3]) # Mark the line core
    ax_line.axhline(ycont, ls='--', color=_COLORS[3]) # Mark the continuum level


    # ########## FOURIER TRANSFORM OF THE LINE ##########
    ft = vfft(wx0, fx0, xc, sigma_FT0, nyq=True)
    ew = 1.e3 * ft.ew # Equivalent width from the zero-frequency amplitude [mA]
    print(f"{'ew':>13} = {ew:.1f} mA")
    ax = _show(fig, 'Fourier transform, mark the position of the first zero',
               [(ft.vsini, ft.log_amp, dict(color=_COLORS[2]))], xlim=[vft_max, 0.5], ylim=[-4., 0.],
               ax=fig.add_subplot(122))
    ax.axvline(ft.vsini_nyq, color='black', ls=_DASH_DOT3) # Nyquist limit
    max_sigma_obs, factor, ew_intprof = ft.vsini_nyq, ft.factor, None


    # ########## FOURIER TRANSFORM OF THE INTRINSIC PROFILE ##########
    """
    - The intrinsic profile of the user is broadened with vsini = 5 km/s and
      the instrument, and its FT is drawn for comparison. As in IDL,
      max_sigma_obs is overwritten by the Nyquist value of this FT.
    - The profile is then scaled to the EW of the observed line.
    """
    if w_intprof is not None:
        ft_int = vfft(*_convol(w_intprof, f_intprof, vsini=5., res0=0.01, resol=res, beta=beta),
                      xc, sigma_FT0, nyq=True)
        ax.plot(ft_int.vsini, ft_int.log_amp, color=_COLORS[4])
        max_sigma_obs, factor, ew_intprof = ft_int.vsini_nyq, ft_int.factor, 1.e3 * ft_int.ew
        print(f"{'ew_intprof':>13} = {ew_intprof:.1f} mA")
        f_intprof = (f_intprof - 1.) * ew / ew_intprof + 1. # Depth scaled to the observed EW


    # ########## VSINI(FT) ##########
    vft = _xhair(ax)[0][0] # The user marks the first zero of the FT
    ax.axvline(vft, color=_COLORS[3], ls='--') # Mark the selected zero of the FT
    fig.canvas.draw_idle()
    print(f"{'vft':>13} = {vft:.1f} km/s")
    print(f"{'sigma(FT)':>13} = {1. / (factor * vft):.3f} s-1\n") # Frequency of the marked zero


    # ########## CLIPPED POINTS ##########
    if cut is not None: # Clipped points removed for the GOF (IDL: fx0 set to 0, keep fx0 > 0)
        fx0[cut] = 0
        wx0, fx0 = wx0[fx0 > 0], fx0[fx0 > 0]


    # #################### OUTPUT ####################
    prep = dict(vmin=vmin, w0=w0, f0=f0, wx0=wx0, fx0=fx0, xc=xc, ew=ew, snr=snr, vft=vft,
                xfou_obs=ft.vsini, yfou_obs=ft.log_amp, max_sigma_obs=max_sigma_obs,
                clip=int(cut is not None))
    if w_intprof is not None:
        prep.update(w_intprof=w_intprof, f_intprof=f_intprof, ew_intprof=ew_intprof)

    return prep
