from __future__ import annotations

# =======================================================================
# ============================== CONSTANTS ==============================
# =======================================================================

# ==============================
"""
    - RSMPL_FAC: resampling factor. If the power-of-two grid would be too
                 coarse (nn/nn_x > RSMPL_FAC), the number of points is doubled.
    - EPS: tolerance, relative to the central wavelength, used to decide
           whether the input spectrum is equidistant.
"""
RSMPL_FAC = 0.7 # Resampling factor
EPS = 2.e-3 / 5000. # Relative tolerance for equidistant wavelength grids




# ===============================================================================
# ============================== INTRINSIC PROFILE ==============================
# ===============================================================================

# ==============================
def _ewgauss(
             delta : np.ndarray[float],
             ew    : float,
             resol : float,
             lam0  : float,
             ycont : float = 1.
    ) -> tuple[np.ndarray[float], np.ndarray[float]]:
    """
    Gaussian absorption profile with a given equivalent width.

    The profile is centred at lam0, has a full width at half maximum
    FWHM = lam0/resol and an integrated depth (equivalent width) equal
    to ew. It is evaluated on a new, regular wavelength grid covering the
    range of delta. IACOB-BROAD uses it with R = 200000 as a delta-like
    intrinsic profile carrying the EW of the observed line.

    ----------
    Parameters
    ----------
        delta : np.ndarray[float]
            Wavelength array [A]. Only its minimum and maximum are used, to
            define the range of the new grid.

        ew : float
            Equivalent width of the profile [A].

        resol : float
            Resolving power defining the width of the Gaussian,
            FWHM = lam0/resol.

        lam0 : float
            Central wavelength of the profile [A].

        ycont : float, optional
            Continuum level. Defaults to 1.

    -------
    Returns
    -------
        wave : np.ndarray[float]
            New wavelength grid of 10010 points: 10000 steps over
            [min(delta), max(delta)] plus 10 extra steps beyond the maximum.

        flux : np.ndarray[float]
            Gaussian profile evaluated on wave.

    -----
    Notes
    -----
        - The Gaussian is written as
            F(lambda) = ycont - A * exp(-b * (lambda - lam0)^2),
        with
            b = ln(2) * (2/FWHM)^2    and    A = ew * sqrt(b/pi),
        so that the integral of the depth, A*sqrt(pi/b), is equal to ew.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # #################### WAVELENGTH GRID ####################
    delta = np.asarray(delta) # Input wavelengths (only the range is used)
    paso = (delta.max() - delta.min()) / 10000. # Step of the new grid: 10000 steps over the input range
    wave = np.arange(10010) * paso + delta.min() # 10010 points, i.e. 10 extra steps beyond the maximum (as in IDL)


    # #################### GAUSSIAN PROFILE ####################
    b = np.log(2.) * (2. / (lam0 / resol))**2 # b = ln2 * (2/FWHM)^2, with FWHM = lam0/resol
    amp = ew * np.sqrt(b / np.pi) # Amplitude such that the integrated depth equals ew
    flux = ycont - amp * np.exp(-b * (wave - lam0)**2) # Absorption profile below the continuum level


    # #################### OUTPUT ####################
    return wave, flux




# =================================================================================
# ============================== BROADENING PROFILES ==============================
# =================================================================================

# ==============================
def _rotation_profile(
                      px    : np.ndarray[float],
                      lam0  : float,
                      vsini : float,
                      beta  : float
    ) -> np.ndarray[float]:
    """
    Rotational broadening profile with linear limb darkening.

    ----------
    Parameters
    ----------
        px : np.ndarray[float]
            Wavelength offsets from the line centre [A], centred on element
            n//2.

        lam0 : float
            Central wavelength [A].

        vsini : float
            Projected rotational velocity [km/s].

        beta : float
            Limb-darkening coefficient, beta = epsilon/(1 - epsilon).

    -------
    Returns
    -------
        profile : np.ndarray[float]
            Normalised rotational profile on px.

    -----
    Notes
    -----
        - With xi = px / (lam0 * vsini / c), the profile is
            G(xi) = [2*sqrt(1 - xi^2)/pi + beta*(1 - xi^2)/2] * normf / (1 + 2*beta/3),
        where normf = c/(vsini*lam0). Outside the stellar disc (|xi| > 1)
        xi^2 is set to 1, so the profile is zero there, as in IDL.

        - The factor 2/3 of the normalisation is written as 0.6666666, as in
          the IDL code.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # ########## PANORAMA ##########
    from .variables import C_KMS_CONVOL


    # #################### ROTATIONAL PROFILE ####################
    normf = C_KMS_CONVOL / (vsini * lam0) # Inverse of the maximum Doppler shift, 1/(lam0*vsini/c) [1/A]
    xi2 = np.minimum((normf * px)**2, 1.) # Normalised coordinate squared; |xi| > 1 (outside the disc) -> xi^2 = 1 -> profile = 0
    profile = ((2. * np.sqrt(1 - xi2) / np.pi + (1 - xi2) * beta / 2.)
               * normf / (1. + 0.6666666 * beta)) # Linear limb-darkening law, normalised to unit area


    # #################### OUTPUT ####################
    return profile




# ==============================
def _rt_component(
                  px2  : np.ndarray[float],
                  lam0 : float,
                  zeta : float,
                  frac : float
    ) -> np.ndarray[float]:
    """
    One component (radial or tangential) of the radial-tangential
    macroturbulence profile, on the positive half of the x-scale.

    ----------
    Parameters
    ----------
        px2 : np.ndarray[float]
            Positive wavelength offsets from the line centre [A].

        lam0 : float
            Central wavelength [A].

        zeta : float
            Macroturbulent velocity of the component [km/s].

        frac : float
            Fraction of the surface covered by this component.

    -------
    Returns
    -------
        component : np.ndarray[float]
            Profile of the component on px2.

    -----
    Notes
    -----
        - With m = zeta*lam0/c and u = px2/m, the component is
            M(u) = 2*frac/(sqrt(pi)*m) * u * [-sqrt(pi) + exp(-u^2)/u + sqrt(pi)*erf(u)],
        i.e. the radial-tangential profile of Gray (1975).
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np
    from scipy.special import erf


    # ########## PANORAMA ##########
    from .variables import C_KMS_CONVOL


    # #################### RADIAL-TANGENTIAL COMPONENT ####################
    sqrt_pi = np.sqrt(np.pi)
    m = zeta * lam0 / C_KMS_CONVOL # Characteristic width of the component [A]
    u = px2 / m # Offsets in units of the characteristic width
    component = 2. * frac / sqrt_pi / m * u * (-sqrt_pi + np.exp(-u**2) / u + sqrt_pi * erf(u))


    # #################### OUTPUT ####################
    return component




# ==============================
def _macroturbulence_profile(
                             px     : np.ndarray[float],
                             lam0   : float,
                             zeta_r : float | None,
                             zeta_t : float | None
    ) -> np.ndarray[float]:
    """
    Radial-tangential macroturbulence profile.

    ----------
    Parameters
    ----------
        px : np.ndarray[float]
            Wavelength offsets from the line centre [A], centred on element
            n//2.

        lam0 : float
            Central wavelength [A].

        zeta_r : float | None
            Radial macroturbulent velocity [km/s]. A value of None (or 0)
            means that the component is not used.

        zeta_t : float | None
            Tangential macroturbulent velocity [km/s]. A value of None (or 0)
            means that the component is not used.

    -------
    Returns
    -------
        profile : np.ndarray[float]
            Symmetric macroturbulence profile on px.

    -----
    Notes
    -----
        - The tangential fraction is 0.5 when both components are given, and
          1 (0) when only zeta_t (zeta_r) is given. IACOB-BROAD always passes
          only zeta_t, so its profile is purely tangential, as in IDL.

        - The profile is computed on the positive half of px and mirrored as
            [reverse(p), p[0], p[0], p],
        i.e. the value p[0] is repeated at the two central points, as in the
        IDL code.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # #################### COMPONENTS ####################
    n = px.size
    aat = 0.5 if (zeta_r and zeta_t) else (1. if zeta_t else 0.) # Tangential fraction of the surface
    px2 = px[n // 2 + 1:] # Positive half of the x-scale (element n//2 is the centre)


    # ########## SUM OF COMPONENTS ##########
    py2 = 0.
    if zeta_r: # Radial component first (IDL: PY2 = PY2R + PY2T)
        py2 = py2 + _rt_component(px2, lam0, zeta_r, 1. - aat)
    if zeta_t: # Tangential component
        py2 = py2 + _rt_component(px2, lam0, zeta_t, aat)


    # #################### SYMMETRIC PROFILE ####################
    profile = np.concatenate([py2[::-1], py2[:1], py2[:1], py2]) # [reverse, p[0], p[0], p]: p[0] repeated at the centre, as in IDL


    # #################### OUTPUT ####################
    return profile




# ==============================
def _gaussian_profile(
                      px     : np.ndarray[float],
                      sigma  : float,
                      centre : float = 0.
    ) -> np.ndarray[float]:
    """
    Normalised Gaussian profile.

    ----------
    Parameters
    ----------
        px : np.ndarray[float]
            Wavelength offsets from the line centre [A].

        sigma : float
            1/e half-width of the Gaussian [A].

        centre : float, optional
            Centre of the Gaussian [A]. Defaults to 0.

    -------
    Returns
    -------
        profile : np.ndarray[float]
            exp(-((px - centre)/sigma)^2) / (sqrt(pi)*sigma), which has unit
            area.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # #################### GAUSSIAN ####################
    profile = np.exp(-((px - centre) / sigma)**2) / (np.sqrt(np.pi) * sigma)


    # #################### OUTPUT ####################
    return profile




# =========================================================================
# ============================== CONVOLUTION ==============================
# =========================================================================

# ==============================
def _conv_spec(
               sy : np.ndarray[float],
               py : np.ndarray[float],
               dx : float
    ) -> np.ndarray[float]:
    """
    FFT convolution of a spectrum with a profile sampled on the same
    equidistant grid (translation of CONV_SPEC in convol.pro).

    ----------
    Parameters
    ----------
        sy : np.ndarray[float]
            Flux on the equidistant wavelength grid.

        py : np.ndarray[float]
            Broadening profile on the same number of points, centred on
            element n//2.

        dx : float
            Step of the wavelength grid [A] (integration element of the
            convolution).

    -------
    Returns
    -------
        conv : np.ndarray[float]
            Convolved flux.

    -----
    Notes
    -----
        - The mean of the two end values of the spectrum is subtracted before
          the (circular) convolution and added back afterwards, so that the
          continuum does not wrap around the edges.

        - The convolution is computed with the convolution theorem,
            conv = dx * IFFT( FFT(sy - level) * FFT(py) ),
        and shifted by n//2 to undo the centring of the profile. The IDL
        expression FFT(FFT(a,-1)*FFT(b,-1),1)*N*dx reduces to this one with
        the NumPy normalisation of the FFT. Real FFTs are used because both
        inputs are real.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # #################### CONVOLUTION ####################
    n = sy.size
    level = (sy[0] + sy[-1]) / 2. # Mean level of the two edges, removed before the convolution
    conv = np.fft.irfft(np.fft.rfft(sy - level) * np.fft.rfft(py), n=n) # Circular convolution through the convolution theorem
    conv = dx * np.roll(conv, n // 2) + level # Undo the centring of the profile, scale by the step and restore the level


    # #################### OUTPUT ####################
    return conv




# ==============================
def _convol(
            sx     : np.ndarray[float],
            sy     : np.ndarray[float],
            resol  : float | None = None,
            vsini  : float | None = None,
            beta   : float | None = None,
            zeta_r : float | None = None,
            zeta_t : float | None = None,
            vgauss : float | None = None,
            vg0    : float = 0.,
            res0   : float | None = None
    ) -> tuple[np.ndarray[float], np.ndarray[float]]:
    """
    Convolve a spectrum with the rotational, radial-tangential
    macroturbulence, instrumental and Gaussian profiles (translation of
    CONVOL in convol.pro).

    Each profile is applied only when its parameter is given, and always in
    the order rotation -> macroturbulence -> instrument -> Gaussian. The
    spectrum is first resampled to an equidistant grid with a power-of-two
    number of points, on which all the profiles are built.

    ----------
    Parameters
    ----------
        sx : np.ndarray[float]
            Input spectrum wavelength [A].

        sy : np.ndarray[float]
            Input spectrum flux.

        resol : float | None, optional
            Resolving power of the instrumental profile.

        vsini : float | None, optional
            Projected rotational velocity [km/s].

        beta : float | None, optional
            Limb-darkening coefficient, beta = epsilon/(1 - epsilon).
            None (or 0, as in IDL) -> BETA = 1.5.

        zeta_r : float | None, optional
            Radial macroturbulent velocity [km/s].

        zeta_t : float | None, optional
            Tangential macroturbulent velocity [km/s].

        vgauss : float | None, optional
            Gaussian broadening velocity (1/e half-width) [km/s].

        vg0 : float, optional
            Centre of the Gaussian profile [km/s]. Defaults to 0.

        res0 : float | None, optional
            Step of the resampled and of the output spectrum [A]. IACOB-BROAD
            always sets it.

    -------
    Returns
    -------
        wave : np.ndarray[float]
            Wavelengths of the convolved spectrum [A], with step res0,
            starting at the first point of the input.

        flux : np.ndarray[float]
            Convolved flux.

    -----
    Notes
    -----
        - Resampling rules of the IDL code:
            >> The new step is res0 (or a fraction RSMPL_FAC of the minimum
            step, for non-equidistant data without res0).
            >> The number of points is raised to the next power of two, and
            doubled if the sampling would still be too coarse
            (nn/nn_x > RSMPL_FAC) for equidistant data.
            >> The new grid starts at the first ORIGINAL wavelength.

        - All the profiles are built on the x-scale px = (i - n//2)*dx [A],
          centred on element n//2, and applied with _conv_spec.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # ########## PANORAMA ##########
    from .numerics import _interpol
    from .variables import BETA, C_KMS_CONVOL


    # #################### INPUT SAMPLING ####################
    # ########## INPUT SPECTRUM ##########
    sx, sy = np.asarray(sx, dtype=float), np.asarray(sy, dtype=float)
    sxx, syy = sx, sy # Working copies (resampled below)
    n = sxx.size


    # ########## EQUIDISTANT GRID? ##########
    step = np.abs(np.diff(sxx)) # Distances between consecutive points
    not_eqdist = (step.max() - step.min()) > EPS * ((sxx[-1] + sxx[0]) / 2.) # Non-equidistant if the steps differ by more than EPS*lambda


    # ########## NEW STEP ##########
    """
    - By default the input sampling is kept (nn points, mean step rdd).
    - With res0 (or non-equidistant data) the new step is res0 (or a fraction
      of the minimum step), and nn becomes the non-integer number of points
      that this step gives over the input range.
    """
    nn, rdd, intp = n, (sxx[-1] - sxx[0]) / (n - 1), False
    if not_eqdist or res0:
        rdd = res0 if res0 else max(step.min() * RSMPL_FAC, 2.7 * EPS * (sxx[-1] + sxx[0]) / 2.)
        nn, intp = (sxx[-1] - sxx[0]) / rdd, True


    # #################### RESAMPLING TO A POWER OF TWO ####################
    xpow = np.log2(nn) # Exponent of the number of points
    ipow = int(xpow)
    if ((2 < ipow <= 15) or intp) and xpow != ipow:
        nn_x = 2**(ipow + 1) # Next power of two
        nn = nn_x * 2 if (nn / nn_x > RSMPL_FAC and not not_eqdist) else nn_x # Doubled if the sampling would be too coarse
        rdd = (sxx[-1] - sxx[0]) / (nn - 1) # Step of the resampled grid
        sxx = np.arange(nn) * rdd + sx[0] # New equidistant grid, starting at the first ORIGINAL wavelength
        syy = _interpol(sy, sx, sxx) # Flux on the new grid (IDL INTERPOL, with extrapolation)


    # #################### PROFILES X-SCALE ####################
    nn = int(nn)
    px = (np.arange(nn) - nn // 2) * rdd # Offsets from the centre [A], centred on element nn//2
    lam0 = (sxx[-1] + sxx[0]) / 2. # Central wavelength of the spectrum
    dx = (sxx[-1] - sxx[0]) / (nn - 1) # Integration step of the convolutions


    # #################### CONVOLUTIONS ####################
    """
    - Order of the convolutions, as in IDL:
        >> rotation -> macroturbulence -> instrument -> Gaussian
    """

    # ########## ROTATION ##########
    if vsini:
        syy = _conv_spec(syy, _rotation_profile(px, lam0, vsini, beta or BETA), dx)


    # ########## MACROTURBULENCE ##########
    if zeta_r or zeta_t:
        syy = _conv_spec(syy, _macroturbulence_profile(px, lam0, zeta_r, zeta_t), dx)


    # ########## INSTRUMENT ##########
    if resol:
        width = lam0 / resol / 2. / np.sqrt(np.log(2.)) # Gaussian with FWHM = lam0/R: 1/e half-width FWHM/(2*sqrt(ln2))
        syy = _conv_spec(syy, _gaussian_profile(px, width), dx)


    # ########## GAUSSIAN ##########
    if vgauss:
        syy = _conv_spec(syy, _gaussian_profile(px, vgauss * lam0 / C_KMS_CONVOL, vg0 * lam0 / C_KMS_CONVOL), dx)


    # #################### OUTPUT SPECTRUM ####################
    nnx = (sxx.max() - sxx.min()) / res0 - 1. # Number of output points with step res0
    wave = sxx.min() + np.arange(int(nnx)) * res0 # Output grid, starting at the first point
    flux = _interpol(syy, sxx, wave) # Convolved flux on the output grid


    # #################### OUTPUT ####################
    return wave, flux
