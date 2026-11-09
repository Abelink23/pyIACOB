from __future__ import annotations
from typing import NamedTuple

# =====================================================================
# ============================== RESULTS ==============================
# =====================================================================

# ==============================
class Envelope(NamedTuple):
    """
    Spline fit to the lower envelope of a chi2 distribution (env_dist).
    """
    x: np.ndarray | None # Variable grid of the spline (None if fewer than 3 distinct values)
    y: np.ndarray | None # Fitted lower envelope of chi2
    resul: dict[str, float] # IDL structure 'resul': bfm, min/max/bfm/err for 1x and 2x, var_min/max




# ==============================
class Estimates(NamedTuple):
    """
    Final estimates of IACOB-BROAD (final_estimates). Errors are [lower, upper].
    """
    vs0: float # vsini (GOF) assuming vmac = 0 [km/s]
    vmft: float # vmac (GOF) for vsini = vsini(FT) [km/s]
    evmft: np.ndarray[float] # 1-sigma errors of vmft [km/s]
    vsgof: float # vsini (GOF) [km/s]
    evsgof: np.ndarray[float] # 1-sigma errors of vsgof [km/s]
    vmgof: float # vmac (GOF) [km/s]
    evmgof: np.ndarray[float] # 1-sigma errors of vmgof [km/s]
    ew: float # Equivalent width (GOF) [mA]




# =====================================================================
# ============================== HELPERS ==============================
# =====================================================================

# ==============================
def _model_window(
                  w0  : np.ndarray[float],
                  xc  : float,
                  vsf : float
    ) -> tuple[np.ndarray[float], float]:
    """
    Wavelengths used to build the synthetic profiles and output step of the
    convolutions.

    ----------
    Parameters
    ----------
        w0 : np.ndarray[float]
            Wavelengths of the observed spectrum around the line [A].

        xc : float
            Wavelength of the line core [A].

        vsf : float
            Maximum vsini of the GOF grid [km/s].

    -------
    Returns
    -------
        w_inst : np.ndarray[float]
            Observed wavelengths within xc +/- 3*vsf/c*xc [A], wide enough to
            contain the broadest synthetic profile of the grid.

        res0 : float
            Output step of the convolutions [A]: 0.1 times the step of the
            observed spectrum.

    -----
    Notes
    -----
        - The step of the observed spectrum is the IDL MEDIAN of
          |w - SHIFT(w, 1)|, i.e. the upper median of the differences, where
          the first difference wraps around (last - first point).
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # ########## PANORAMA ##########
    from .numerics import _median_upper
    from .variables import C_KMS


    # #################### WINDOW & STEP ####################
    w_inst = w0[np.abs(w0 - xc) <= xc * 3. * vsf / C_KMS] # Window of +/- 3 times the maximum Doppler shift of the grid
    t = _median_upper(np.abs(w_inst - np.roll(w_inst, 1))) # Step of the observed spectrum (IDL MEDIAN of the differences)
    res0 = 0.1 * t # Output step of the convolutions


    # #################### OUTPUT ####################
    return w_inst, res0




# ==============================
def _rv_adjust(
               w_inst : np.ndarray[float],
               f_inst : np.ndarray[float],
               wx0    : np.ndarray[float],
               fx0    : np.ndarray[float],
               vft    : float,
               res0   : float,
               resol  : float,
               beta   : float | None = None
    ) -> np.ndarray[float]:
    """
    Shift the intrinsic profile of the user to the radial velocity that
    maximises its correlation with the observed line.

    ----------
    Parameters
    ----------
        w_inst : np.ndarray[float]
            Wavelengths of the intrinsic profile of the user, on the model wavelengths.

        f_inst : np.ndarray[float]
            Fluxes of the intrinsic profile of the user, on the model wavelengths.

        wx0 : np.ndarray[float]
            Wavelengths of the observed line profile.

        fx0 : np.ndarray[float]
            Fluxes of the observed line profile.

        vft : float
            vsini from the first zero of the FT [km/s].

        res0 : float
            Output step of the convolutions [A].

        resol : float
            Resolving power of the spectrum.

        beta : float | None, optional
            Limb-darkening coefficient (None -> 1.5).

    -------
    Returns
    -------
        f_shift : np.ndarray[float]
            Intrinsic profile shifted to the best radial velocity, on w_inst.

    -----
    Notes
    -----
        - The intrinsic profile is broadened with vsini(FT) and the instrument,
          and cross-correlated (Pearson coefficient at lag 0, IDL C_CORRELATE)
          with the observed profile for radial-velocity shifts from -RV_RANGE
          to +RV_RANGE km/s in steps of RV_STEP km/s.

        - The shift with the maximum correlation (the first one, if repeated)
          is applied to the unbroadened intrinsic profile:
            lambda' = lambda * (1 - rv/c).
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # ########## PANORAMA ##########
    from .numerics import _interpol
    from .profiles import _convol
    from .variables import C_KMS, RV_RANGE, RV_STEP


    # #################### BROADENED INTRINSIC PROFILE ####################
    w_rot, f_rot = _convol(w_inst, f_inst, vsini=vft, res0=res0, resol=resol, beta=beta) # Rotation (vsini(FT)) + instrument
    model = _interpol(f_rot, w_rot, wx0) # On the observed wavelengths


    # #################### CROSS-CORRELATION ####################
    nrv = int(round(RV_RANGE / RV_STEP))
    rv = np.arange(-nrv, nrv + 1) * RV_STEP # Radial-velocity shifts: -RV_RANGE, ..., +RV_RANGE [km/s]
    corr = [np.corrcoef(_interpol(model, wx0 * (1. - r / C_KMS), wx0), fx0)[0, 1] for r in rv] # Correlation of the shifted model with the observed line
    best = rv[int(np.argmax(corr))] # First maximum of the correlation


    # #################### SHIFTED PROFILE ####################
    f_shift = _interpol(f_inst, w_inst * (1. - best / C_KMS), w_inst) # Best shift applied to the unbroadened intrinsic profile


    # #################### OUTPUT ####################
    return f_shift




# ======================================================================
# ============================== GOF GRID ==============================
# ======================================================================

# ==============================
def _gof_grid(
              prep     : dict,
              resol    : float,
              fact_vrt : float,
              vmin     : float | None = None,
              em       : bool = False,
              beta     : float | None = None
    ) -> dict[str, float | np.ndarray[float]]:
    """
    chi2 of a grid of synthetic profiles (vsini x vmac x EW) with respect to
    the observed line profile.

    ----------
    Parameters
    ----------
        prep : dict
            Pre-processed line (content of <star>_<line>.idl, see _vsiniFT).

        resol : float
            Resolving power of the spectrum.

        fact_vrt : float
            Maximum vmac of the grid relative to vsini(FT).

        vmin : float | None, optional
            Minimum vsini and vmac of the grid [km/s]. Only used if prep does
            not contain it (files of old IDL versions).

        em : bool, optional
            Emission line: the delta-like intrinsic profile is flipped around
            the continuum. Defaults to False.

        beta : float | None, optional
            Limb-darkening coefficient (None -> 1.5).

    -------
    Returns
    -------
        gof : dict
            Results of the grid, with the models sorted by increasing chi2:

            - 'sigma' : float
                Noise of the normalised flux, 1/SNR.

            - 'sigma_chi' : np.ndarray[float]
                Standard deviation of the residuals of each model.

            - 'chi2' : np.ndarray[float]
                chi2 per point of each model.

            - 'chi2_vs' : np.ndarray[float]
                vsini [km/s] of each model.

            - 'chi2_vrt' : np.ndarray[float]
                vmac [km/s] of each model.

            - 'chi2_ew' : np.ndarray[float]
                EW [mA] of each model.

            - 'avs' : float
                Step of the vsini and vmac grids [km/s].

            - 'vrtf' : Only with an intrinsic profile of the user (as in IDL).

            - 'vsf' : Only with an intrinsic profile of the user (as in IDL).

            - 'w_intprof' : Only with an intrinsic profile of the user (as in IDL).

            - 'f_intprof' : Only with an intrinsic profile of the user (as in IDL).

    -----
    Notes
    -----
        - Grids:
            >> vsini and vmac from vmin to FACT_VS*vft and fact_vrt*vft, in
               steps of STEP_V*vft.
            >> EW: EW0*(1 + k*EW_STEP), with k from -n_down to n_up (N_EW for
               normal lines, N_EW_CUT for clipped lines).

        - Intrinsic profile:
            >> Delta-like Gaussian with the EW of each step of the EW grid
               (_ewgauss with R = R_DELTA). As in IDL, _ewgauss replaces the
               wavelength grid, so the grid changes at every EW step.
            >> Or the intrinsic profile of the user, shifted in radial velocity
               to match the observed line (_rv_adjust).

        - For each model, the intrinsic profile is convolved with rotation,
          macroturbulence and the instrument, interpolated on the observed
          wavelengths, and
            chi2 = sum(((F_model - F_obs)/sigma)^2) / N.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    from itertools import product

    import numpy as np


    # ########## PANORAMA ##########
    from .numerics import _interpol
    from .profiles import _convol, _ewgauss
    from .variables import EW_STEP, FACT_VS, N_EW, N_EW_CUT, R_DELTA, STEP_V


    # #################### OBSERVED LINE ####################
    wx0, fx0, xc, ew, vft = prep['wx0'], prep['fx0'], prep['xc'], prep['ew'], prep['vft']
    vmin = prep.get('vmin', vmin) # vmin of the pre-processing (the keyword only for old files)
    w_intprof = prep.get('w_intprof') # Intrinsic profile of the user (None -> delta-like profile)
    sigma = 1. / prep['snr'] # Noise of the normalised flux


    # #################### GRIDS ####################
    # ########## VSINI & VMAC ##########
    avs = STEP_V * vft # Grid step [km/s]
    vsf, vrtf = FACT_VS * vft, fact_vrt * vft # Maximum vsini and vmac of the grid [km/s]
    vs = vmin + np.arange(int(1 + (vsf - vmin) / avs)) * avs # vsini grid: vmin, vmin + avs, ... <= vsf
    vrt = vmin + np.arange(int(1 + (vrtf - vmin) / avs)) * avs # vmac grid: vmin, vmin + avs, ... <= vrtf


    # ########## EQUIVALENT WIDTH ##########
    n_down, n_up = N_EW_CUT if prep['clip'] == 1 else N_EW # Wider EW range when part of the line was clipped
    ews = ew * (1. + np.arange(-n_down, n_up + 1) * EW_STEP) # EW0*(1 + k*EW_STEP) [mA]
    ew_in, ew_fin = ews[0], ews[-1]

    print(f"{'vsini grid':>13} = {vmin:.1f} - {vsf:.1f} km/s, step avs = {avs:.1f} km/s (vmin, vsf, avs)")
    print(f"{'vmac grid':>13} = {vmin:.1f} - {vrtf:.1f} km/s, step avs = {avs:.1f} km/s (vmin, vrtf, avs)")
    print(f"{'ew grid':>13} = {ew_in:.1f} - {ew_fin:.1f} mA, step = {EW_STEP * ew:.1f} mA (ew = {ew:.1f} mA)\n")


    # #################### INTRINSIC PROFILE ####################
    w_inst, res0 = _model_window(prep['w0'], xc, vsf) # Model wavelengths and output step of the convolutions
    if w_intprof is not None: # User profile, shifted to the best radial velocity
        print('Using an intrinsic profile provided by the user')
        f_inst = _rv_adjust(w_inst, _interpol(prep['f_intprof'], w_intprof, w_inst), wx0, fx0, vft, res0, resol, beta)
    else: # Delta-like profile, built for each EW below
        print('Using a delta profile as intrinsic profile')


    # #################### GOF COMPUTATION ####################
    """
    - Order of the loops, as in IDL: EW (outer) -> vsini -> vmac (inner).
    """
    print("GOF ... Rot+Mac(RT) ... Please, wait for a while ...")
    rows = []
    for kew in ews:
        # ########## DELTA-LIKE PROFILE ##########
        if w_intprof is None:
            w_inst, f_inst = _ewgauss(w_inst, kew * 1.e-3, R_DELTA, xc) # Gaussian with this EW [A]; replaces the wavelength grid, as in IDL
            if em:
                f_inst = 1. - (f_inst - 1.) # Emission line: profile flipped around the continuum
            print(f"{'ew':>13} = {kew:.1f} mA (ew(FT) = {ew:.1f} mA)")


        # ########## VSINI x VMAC ##########
        for vs_i, vrt_j in product(vs, vrt):
            w_mod, f_mod = _convol(w_inst, f_inst, vsini=vs_i, res0=res0, resol=resol, beta=beta, zeta_t=vrt_j) # Synthetic profile
            resid = _interpol(f_mod, w_mod, wx0) - fx0 # Residuals on the observed wavelengths
            rows.append((np.std(resid, ddof=1), # sigma_chi: standard deviation of the residuals
                         np.sum((resid / sigma)**2) / wx0.size, # chi2 per point
                         vs_i, vrt_j, kew))


    # #################### SORTING BY CHI2 ####################
    sigma_chi, chi2, chi2_vs, chi2_vrt, chi2_ew = np.array(rows).T
    order = np.argsort(chi2, kind='stable') # Increasing chi2 (stable for identical values)
    gof = dict(sigma=sigma, sigma_chi=sigma_chi[order], chi2=chi2[order], chi2_vs=chi2_vs[order],
               chi2_vrt=chi2_vrt[order], chi2_ew=chi2_ew[order], avs=avs)
    if w_intprof is not None: # As in IDL, the user profile and the grid limits are also saved
        gof.update(vrtf=vrtf, vsf=vsf, w_intprof=w_intprof, f_intprof=prep['f_intprof'])


    # #################### OUTPUT ####################
    return gof




# =================================================================================
# ============================== ESTIMATES FROM CHI2 ==============================
# =================================================================================

# ==============================
def _env_dist(
              var  : np.ndarray[float],
              chi2 : np.ndarray[float],
              curv : float = 5.
    ) -> Envelope:
    """
    Best value and 1-sigma / 2-sigma limits of a variable from its chi2
    distribution (translation of env_dist).

    ----------
    Parameters
    ----------
        var : np.ndarray[float]
            Value of the variable (vsini, vmac or EW) of each model.

        chi2 : np.ndarray[float]
            chi2 of each model.

        curv : float, optional
            Tension of the spline (high values -> nearly linear fit).
            Defaults to 5.

    -------
    Returns
    -------
        envelope : Envelope
            Spline grid (x), fitted envelope (y) and result dict (resul):
            bfm, min/max/bfm/err for 1x and 2x, var_min/max.

    -----
    Notes
    -----
        - The lower envelope is the minimum chi2 for each distinct value of the
          variable. It is fitted with a spline under tension evaluated on 1000
          points over [min(var), max(var)).

        - Keys of resul:
            >> bfm: value of the best-fitting model (absolute minimum).
            >> bfm_1x: minimum of the spline within the 1-sigma region
               (below min + 1).
            >> min_1x: lower limit of the 1-sigma region.
            >> max_1x: upper limit of the 1-sigma region.
            >> errmin_1x: lower 1-sigma error, abs(min_1x - bfm_1x).
            >> errmax_1x: upper 1-sigma error, abs(max_1x - bfm_1x).
            >> err_1x: 1-sigma error, abs(min_1x - bfm_1x).
            >> bfm_2x: minimum of the spline within the 2-sigma region
               (below min + 4).
            >> min_2x: lower limit of the 2-sigma region.
            >> max_2x: upper limit of the 2-sigma region.
            >> errmin_2x: lower 2-sigma error, abs(min_2x - bfm).
            >> errmax_2x: upper 2-sigma error, abs(max_2x - bfm).
            >> err_2x: 2-sigma error, abs(min_2x - bfm_2x).
            >> var_min: minimum value of the variable.
            >> var_max: maximum value of the variable.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # ########## PANORAMA ##########
    from .numerics import _spline_tension


    # #################### INITIALISATION ####################
    tags = ('bfm', 'min_1x', 'max_1x', 'bfm_1x', 'errmin_1x', 'errmax_1x', 'err_1x',
            'min_2x', 'max_2x', 'bfm_2x', 'errmin_2x', 'errmax_2x', 'err_2x', 'var_min', 'var_max')
    var, chi2 = np.asarray(var, dtype=float), np.asarray(chi2, dtype=float)
    res = dict.fromkeys(tags, 0.) # IDL structure 'resul' (all tags initialised to 0)
    res['var_min'], res['var_max'] = var.min(), var.max()


    # #################### LOWER ENVELOPE ####################
    var_x, inv = np.unique(var, return_inverse=True) # Distinct values of the variable
    var_y = np.full(var_x.size, np.inf)
    np.minimum.at(var_y, inv, chi2) # Minimum chi2 for each distinct value


    # #################### BEST-FITTING MODEL ####################
    res['bfm'] = var[np.argmin(chi2)] # First model with the absolute minimum of chi2
    if var_x.size == 1: # Single distinct value: that value is the result
        res['bfm_1x'] = res['bfm_2x'] = var_x[0]
    if var_x.size < 3: # Not enough values for the spline
        return Envelope(None, None, res)


    # #################### SPLINE FIT ####################
    x_sp = var.min() + np.arange(1000) * (var.max() - var.min()) / 1000. # 1000 points over [min, max)
    y_sp = _spline_tension(var_x, var_y, x_sp, curv) # Spline under tension through the envelope
    best = x_sp[np.argmin(y_sp)] # Minimum of the spline


    # #################### UNCERTAINTIES ####################
    # ########## 1 SIGMA: MIN + 1 ##########
    if (y_sp <= chi2.min() + 1.).any():
        sel = x_sp[y_sp <= y_sp.min() + 1.] # Region where the spline is below its minimum + 1
        res.update(min_1x=sel.min(), max_1x=sel.max(), bfm_1x=best)
        res.update(errmin_1x=abs(res['min_1x'] - best), errmax_1x=abs(res['max_1x'] - best),
                   err_1x=abs(res['min_1x'] - best))


    # ########## 2 SIGMA: MIN + 4 ##########
    sel = x_sp[y_sp <= chi2.min() + 4.] # Region where the spline is below min(chi2) + 4
    if sel.size:
        res.update(min_2x=sel.min(), max_2x=sel.max(), bfm_2x=best)
        res.update(errmin_2x=abs(res['min_2x'] - res['bfm']), errmax_2x=abs(res['max_2x'] - res['bfm']),
                   err_2x=abs(res['min_2x'] - best))


    # #################### OUTPUT ####################
    return Envelope(x_sp, y_sp, res)




# ==============================
def _final_estimates(
                     gof : dict,
                     vft : float
    ) -> Estimates:
    """
    vsini and vmac estimates from the GOF distributions (translation of
    final_estimates).

    ----------
    Parameters
    ----------
        gof : dict
            GOF results (chi2, chi2_vs, chi2_vrt, chi2_ew, avs; see _gof_grid).

        vft : float
            vsini from the first zero of the FT [km/s].

    -------
    Returns
    -------
        estimates : Estimates
            - vs0: float
                vsini from the models with the smallest vmac of the grid (vmac ~ 0).

            - vmft: float
                vmac from the models within half a grid step of vsini(FT).

            - evmft: np.ndarray[float]
                1-sigma errors of vmft.

            - vsgof: float
                vsini from the full grid.

            - evsgof: np.ndarray[float]
                1-sigma errors of vsgof.

            - vmgof: float
                vmac from the full grid.

            - ew: float
                Equivalent width from the full grid.

            NOTE: All of them are the bfm_1x values of _env_dist, and the errors are [errmin_1x, errmax_1x].
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # #################### DISTRIBUTIONS ####################
    chi2, vs, vrt = gof['chi2'], gof['chi2_vs'], gof['chi2_vrt']

    def _err(r): # [lower, upper] 1-sigma errors of a result dict
        return np.array([r['errmin_1x'], r['errmax_1x']])


    # ########## VSINI WITH VMAC = 0 ##########
    m0 = vrt == vrt.min() # Models with the smallest vmac of the grid
    vs0 = _env_dist(vs[m0], chi2[m0]).resul['bfm_1x']


    # ########## VMAC WITH VSINI = VSINI(FT) ##########
    mft = (vs >= vft - 0.5 * gof['avs']) & (vs <= vft + 0.5 * gof['avs']) # Models within half a grid step of vsini(FT)
    rft = _env_dist(vrt[mft], chi2[mft]).resul


    # ########## FULL GRID ##########
    rvs = _env_dist(vs, chi2).resul # vsini
    rvm = _env_dist(vrt, chi2).resul # vmac
    rew = _env_dist(gof['chi2_ew'], chi2).resul # EW


    # #################### OUTPUT ####################
    return Estimates(vs0, rft['bfm_1x'], _err(rft), rvs['bfm_1x'], _err(rvs),
                     rvm['bfm_1x'], _err(rvm), rew['bfm_1x'])




# =================================================================================================
# ============================== SYNTHETIC PROFILES OF THE SOLUTIONS ==============================
# =================================================================================================

# ==============================
def _best_models(
                 prep    : dict,
                 est     : Estimates,
                 vsf     : float,
                 resol   : float,
                 gauss   : bool = False,
                 em      : bool = False,
                 beta    : float | None = None,
                 vs_user : float = 0.,
                 vm_user : float = 0.,
                 sigma_FT0 : float | None = None
    ) -> tuple[dict, dict]:
    """
    Intrinsic profile and synthetic profiles of the solutions found, and
    their Fourier transforms (FIGURE 1 and 2 blocks of iacob_broad_v9.pro).

    ----------
    Parameters
    ----------
        prep : dict
            Pre-processed line (see _vsiniFT).

        est : Estimates
            Final estimates (see _final_estimates).

        vsf : float
            Maximum vsini of the GOF grid [km/s].

        resol : float
            Resolving power of the spectrum.

        gauss : bool, optional
            Gaussian instead of radial-tangential macroturbulence.

        em : bool, optional
            Emission line (the delta-like profile is flipped).

        beta : float | None, optional
            Limb-darkening coefficient (None -> 1.5).

        vs_user : float, optional
            vsini of the user [km/s] (0 -> not computed).

        vm_user : float, optional
            vmac of the user [km/s] (0 -> not computed).

        sigma_FT0 : float | None, optional
            First zero of the rotational FT (None -> 0.660).

    -------
    Returns
    -------
        models : dict
            (wavelength, flux) of each profile:
            - 'inst'  : intrinsic profile,
            - 'ft'    : vsini(FT), vmac = 0,
            - 'gof0'  : vsini(GOF, vmac = 0),
            - 'gof'   : vsini(GOF) + vmac(GOF),
            - 'ftgof' : vsini(FT) + vmac(GOF for vsini(FT)),
            - 'user'  : vs_user + vm_user (only if vs_user is given).

        fts : dict
            Fourier transform (FT named tuple) of the profiles 'ft', 'gof',
            'ftgof', 'user', and 'inst' for an intrinsic profile of the user.

    -----
    Notes
    -----
        - The FT of each synthetic profile is computed as for the observed line:
          the region below 0.999 (plus 10% on each side) is resampled on the
          observed wavelengths and transformed with vfft.

        - As in IDL, vs_user <= 5 km/s means no rotation, and the user model
          uses radial-tangential macroturbulence even with gauss.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # ########## PANORAMA ##########
    from .fourier import vfft
    from .numerics import _interpol
    from .profiles import _convol, _ewgauss
    from .variables import R_DELTA


    # #################### INTRINSIC PROFILE ####################
    xc, w0, vft = prep['xc'], prep['w0'], prep['vft']
    w_inst, res0 = _model_window(w0, xc, vsf) # Model wavelengths and output step of the convolutions
    if 'w_intprof' in prep: # User profile on the model grid, with its minimum moved to the line core
        f_inst = _interpol(prep['f_intprof'], prep['w_intprof'], w_inst)
        w_inst = w_inst - w_inst[np.argmin(f_inst)] + xc
    else: # Delta-like Gaussian with the EW(GOF) [A] (replaces the wavelength grid)
        w_inst, f_inst = _ewgauss(w_inst, est.ew * 1.e-3, R_DELTA, xc)
        if em:
            f_inst = 1. - (f_inst - 1.) # Emission line: flipped around the continuum


    # #################### SYNTHETIC PROFILES ####################
    macro = 'vgauss' if gauss else 'zeta_t' # Macroturbulence profile of the GOF

    def _model(vsini, vmac=None, kind=macro): # Rotation + macroturbulence + instrument
        kw = {kind: vmac} if vmac is not None else {}
        return _convol(w_inst, f_inst, vsini=vsini, res0=res0, resol=resol, beta=beta, **kw)

    models = {'inst': (w_inst, f_inst),
              'ft': _model(vft), # vsini(FT), vmac = 0
              'gof0': _model(est.vs0), # vsini(GOF) with vmac = 0
              'gof': _model(est.vsgof, est.vmgof), # 2D GOF solution
              'ftgof': _model(vft, est.vmft)} # vsini(FT) + vmac(GOF)
    if vs_user: # User values: vs_user <= 5 km/s -> no rotation; zeta_t even with gauss (as in IDL)
        models['user'] = _model(vs_user if vs_user > 5. else 0, vm_user, 'zeta_t')


    # #################### FOURIER TRANSFORMS ####################
    def _ft(w, f): # FT of the region below 0.999 (+10%), resampled on the observed wavelengths
        core = w[f <= 0.999]
        pad = 0.1 * (core.max() - core.min())
        wref = w0[(w0 >= core.min() - pad) & (w0 <= core.max() + pad)]
        return vfft(wref, _interpol(f, w, wref), xc, sigma_FT0)

    fts = {k: _ft(*models[k]) for k in ('ft', 'gof', 'ftgof', 'user') if k in models}
    if 'w_intprof' in prep:
        fts['inst'] = _ft(w_inst, f_inst)


    # #################### OUTPUT ####################
    return models, fts
