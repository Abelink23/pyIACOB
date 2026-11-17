from __future__ import annotations

# ===============================================================================
# ============================== STELLAR PARAMETERS ==============================
# ===============================================================================

# ==============================
def _rlm(
         sol : np.ndarray,
         Mv  : float | None = None,
         Mk  : float | None = None,
         kud : bool = False
    ) -> None:
    """
    Radius, luminosity and spectroscopic mass of each model from the
    absolute magnitude (translation of gbat_rlm.pro). sol is modified in
    place (fields r [Rsun], l [log L/Lsun] and m [Msun]).

    -----
    Notes
    -----
        - R = 10^(-0.2 (Mv - V_syn)) (or with Ks and Mk), with V the synthetic
          magnitude of the model for R = 1 Rsun; with 'kud', R is computed
          with Kudritzki's formula, R = 10^((29.57 - (Mv - vsyn)) / 5).
        - M = g R^2 / G [Msun], log L = 2 log R + 4 log(Teff/Tsun).
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # ########## PANORAMA ##########
    from gbat.variables import G_CGS, MSUN, RSUN, TSUN


    # #################### RADIUS ####################
    glogsun = np.log10(G_CGS * MSUN / RSUN / RSUN)
    if Mv:
        r = 10.**((29.57 - (Mv - sol['vsyn'])) / 5.) if kud else 10.**(-0.2 * (Mv - sol['v']))
    elif Mk:
        r = 10.**(-0.2 * (Mk - sol['ks']))
    else:
        return


    # #################### MASS & LUMINOSITY ####################
    sol['r'] = r
    sol['m'] = 10.**(sol['grav'] - glogsun + 2. * np.log10(r))
    sol['l'] = 2. * np.log10(r) + 4. * np.log10(sol['teff'] * 1e3 / TSUN)




# ==============================
def _massloss(
              sol  : np.ndarray,
              vinf : float
    ) -> None:
    """
    Mass-loss rate, -log Mdot, of each model from log Q, R and vinf
    (translation of gbat_massloss.pro). sol is modified in place.

    -----
    Notes
    -----
        - Q = Mdot / (R vinf)^1.5 -> -log Mdot = |-|log Q| + 1.5 log(R vinf)|.
    """

    # #################### MODULES & PACKAGES ####################
    import numpy as np


    # #################### MASS LOSS ####################
    sol['mdot'] = np.abs(-np.abs(sol['q']) + 1.5 * np.log10(sol['r'] * vinf))




# ==============================
def _gamma_rad(
               teff : float,
               logg : float,
               yhe  : float,
               z    : float
    ) -> float:
    """
    Eddington factor for electron scattering (translation of
    gbat_gamma_rad.pro, M. A. Urbaneja 2006; not used by IACOB-GBAT).
    """

    # #################### CONSTANTS ####################
    c = 2.9979250e10 # [cm/s]
    sigma = 5.67e-5  # Stefan-Boltzmann constant [cgs]
    met_sun = 0.02


    # #################### OUTPUT ####################
    x = (1. - z * met_sun) / (1. + 4. * yhe) # H mass fraction
    return (sigma / c) * 0.20 * (1. + x) * teff**4. / (10.**logg)




# ==============================
def _vescape(
             teff  : float,
             logg  : float,
             r     : float,
             y     : float,
             z     : float = 1.,
             gamma : bool = False
    ) -> float:
    """
    Effective escape velocity [km/s] (translation of gbat_vescape.pro, M. A.
    Urbaneja 2006; not used by IACOB-GBAT).

    -----
    Notes
    -----
        - The solar radius 6.6599e10 cm is the one of the IDL code.
    """

    # #################### MODULES & PACKAGES ####################
    import numpy as np


    # #################### ESCAPE VELOCITY ####################
    rsun = 6.6599e10
    g = _gamma_rad(teff, logg, y, z)
    if gamma:
        print('Gamma : ', g)


    # #################### OUTPUT ####################
    return np.sqrt(2. * 10.**logg * (r * rsun) * (1. - g)) * 1e-5




# =====================================================================================
# ============================== CHI2 DISTRIBUTION OF ONE PARAMETER ==============================
# =====================================================================================

# ==============================
def _labels(
            Mv   : float | None,
            Mk   : float | None,
            vinf : float | None
    ) -> str:
    """
    Labels of the file names with Mv, Mk and vinf (IDL COMMON share3:
    lab_Mv + lab_Mk + lab_vinf), e.g. '_Mv489_vinf2715' for Mv = -4.89 and
    vinf = 2715 km/s.
    """

    # #################### PANORAMA ####################
    from gbat.io import _fi


    # #################### LABELS ####################
    lab_mv = '_Mv' + _fi(-100. * Mv, 4).strip() if Mv else ''
    lab_mk = '_Mk' + _fi(-100. * Mk, 4).strip() if Mk else ''
    lab_vinf = '_vinf' + _fi(vinf, 4).strip() if vinf else ''


    # #################### OUTPUT ####################
    return lab_mv + lab_mk + lab_vinf




# ==============================
def _plot_dist_one(
                   var   : np.ndarray[float],
                   sol   : np.ndarray,
                   title : str,
                   extra : float,
                   base  : str | Path
    ) -> dict:
    """
    Central value and uncertainties of one parameter from the lower envelope
    of its chi2 distribution (computational part of gbat_plot_dist_one.pro).

    ----------
    Parameters
    ----------
        var : np.ndarray[float]
            Value of the parameter for each model (sol sorted by chi2).

        sol : np.ndarray
            Models sorted by chi2.

        title : str
            Title of the parameter (IDL text format), stored in resul.

        extra : float
            Minimum separation between consecutive points of the envelope.

        base : str | Path
            Save file without extension (xdr_results/<obs_run>/<label>).

    -------
    Returns
    -------
        dist : dict
            resul (structure with the results), var_x, var_y (points of the
            lower envelope) and varx_sp, vary_sp (spline fit to them).

    -----
    Notes
    -----
        - Lower envelope: chi2 is scanned in steps of 1 from its minimum; at
          each step, the lowest-chi2 model below the minimum (above the
          maximum) of the envelope by more than 'extra' is added.
        - Fit: spline under tension (CURV = 80, IDL SPLINE) of the envelope on
          NSTEP + 1 points, or a straight line for 2 points.
        - 1 (2) sigma: range of the fit with chi2 <= min(chi2) + 1 (+ 4).
          resul gives min, max, mean (bfm) of the range and the errors with
          respect to the best-fitting model.
        - Saved in '<base>.idl' with resul and var_x.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # ########## PANORAMA ##########
    from broadening.numerics import _spline_tension
    from gbat.io import _fa, _ff, _resul_init, _save_gbat
    from gbat.variables import CURV, NSTEP


    # #################### RESUL ####################
    var = np.asarray(var, dtype=float)
    chi = sol['chi']
    cmin = chi.min()
    resul = _resul_init()
    r = resul[0]
    r['title'] = title
    r['bfm'] = var[np.nonzero(chi <= cmin)[0][0]]
    print(title + ' : ' + _ff(r['bfm'], 5, 2) + '  ')


    # #################### LOWER ENVELOPE ####################
    var_x, var_y = [var[0]], [chi[0]]
    val = cmin
    while val <= chi.max():
        tt = np.nonzero(chi <= val)[0]
        ss1 = np.nonzero(var[tt] < min(var_x) - extra)[0]
        if ss1.size:
            var_x.append(var[ss1[0]])
            var_y.append(chi[ss1[0]])
        ss2 = np.nonzero(var[tt] > max(var_x) + extra)[0]
        if ss2.size:
            var_x.append(var[ss2[0]])
            var_y.append(chi[ss2[0]])
        val += 1.
    order = np.argsort(var_x, kind='stable')
    var_x, var_y = np.asarray(var_x)[order], np.asarray(var_y)[order]


    # #################### SPLINE FIT ####################
    r['var_max'], r['var_min'] = var_x.max(), var_x.min()
    sh = (r['var_max'] - r['var_min']) / NSTEP
    varx_sp = r['var_min'] + np.arange(NSTEP + 1) * sh
    if var_x.size == 1:
        for k in ('min', 'max', 'bfm'):
            r[k + '_1x'] = r[k + '_2x'] = var_x[0]
        varx_sp, vary_sp = var_x, var_y
    elif var_x.size == 2:
        coef = np.polyfit(var_x, var_y, 1)
        vary_sp = np.polyval(coef, varx_sp)
    else:
        vary_sp = _spline_tension(var_x, var_y, varx_sp, sigma=CURV)


    # #################### UNCERTAINTIES (1 AND 2 SIGMA) ####################
    if var_x.size >= 2:
        for lab, dchi in (('1x', 1.), ('2x', 4.)):
            ff = vary_sp <= cmin + dchi
            if ff.any():
                r['min_' + lab] = varx_sp[ff].min()
                r['max_' + lab] = varx_sp[ff].max()
                r['bfm_' + lab] = 0.5 * (varx_sp[ff].min() + varx_sp[ff].max())
                r['errmin_' + lab] = abs(r['min_' + lab] - r['bfm'])
                r['errmax_' + lab] = abs(r['max_' + lab] - r['bfm'])
                r['err_' + lab] = abs(r['min_' + lab] - r['bfm_' + lab])

    fmt3 = lambda a, b, c: _ff(a, 5, 2) + '  ' + _ff(b, 5, 2) + '  ' + _ff(c, 5, 2) + '  '
    print(_fa('min, max, mean [1sigma]:') + fmt3(r['min_1x'], r['max_1x'], r['bfm_1x']))
    print(_fa('min, max, mean [2sigma]:') + fmt3(r['min_2x'], r['max_2x'], r['bfm_2x']))


    # #################### SAVE ####################
    _save_gbat(base, resul=resul, var_x=var_x)


    # #################### OUTPUT ####################
    return dict(resul=resul, var_x=var_x, var_y=var_y, varx_sp=varx_sp, vary_sp=vary_sp)




# ==============================
def _write_resul_one(
                     file : str | Path,
                     fmt  : str
    ) -> dict:
    """
    Final value and uncertainty of one parameter, as written in the figures
    and in the summary table (translation of gbat_write_resul_one.pro).

    ----------
    Parameters
    ----------
        file : str | Path
            Save file of plot_dist_one (.idl, or IDL .xdr).

        fmt : str
            IDL format of the values: '(F4.1)', '(F6.2)', '(I2)' ...

    -------
    Returns
    -------
        out : dict
            - title : title of the parameter (IDL text format),
            - lab : central value used in the figures (min/max of the 1-sigma
                    range for limits),
            - lab_sum : text, e.g. '= 33.4 +/- 0.8', '> 3.70 -/- 0.10' or
                        '< 13.00 +/+ 0.15',
            - lab_bfm : value of the best-fitting model,
            - bfm_text : best-fitting model value formatted with fmt.

    -----
    Notes
    -----
        - For log Q (title containing 'Q') the values are made negative and
          the lower and upper errors swapped.
        - '=' central value (bfm_1x) +/- error (1-sigma half range, or half
          the step of the grid around the best model if larger); '>' or '<'
          when the 1-sigma range reaches the edge of the subgrid; a range
          'min-max' when it reaches both edges; '= value' for a single value.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import re
    import numpy as np


    # ########## PANORAMA ##########
    from gbat.io import RESUL_FIELDS, _ff, _fi, _restore_gbat


    # #################### FORMAT ####################
    code = re.sub(r'[()\s]', '', fmt).upper()
    if code.startswith('I'):
        w = int(code[1:])
        fv = lambda x: _fi(x, w).strip()
    else:
        w, d = (int(v) for v in code[1:].split('.'))
        fv = lambda x: _ff(x, w, d).strip()


    # #################### RESTORE ####################
    data = _restore_gbat(file)
    f32 = lambda v: float(np.float32(v)) # resul is single precision in IDL
    resul = {k: (str(v) if k == 'title' else f32(v)) for k, v in zip(RESUL_FIELDS, data['resul'][0])}
    var_x = np.asarray(data['var_x'], dtype=float)


    # #################### LOG Q ####################
    if 'Q' in resul['title']:
        for k in RESUL_FIELDS[1:7] + RESUL_FIELDS[10:13]:
            resul[k] = f32(-1. * abs(resul[k]))
        var_x = -1. * np.abs(var_x)
        resul['errmin_1x'], resul['errmax_1x'] = resul['errmax_1x'], resul['errmin_1x']
        resul['errmin_2x'], resul['errmax_2x'] = resul['errmax_2x'], resul['errmin_2x']

    r0 = dict(resul)
    for a, b in (('var_min', 'var_max'), ('min_1x', 'max_1x'), ('min_2x', 'max_2x')):
        resul[a], resul[b] = min(r0[a], r0[b]), max(r0[a], r0[b])


    # #################### TEXT ####################
    r = resul
    if r['var_min'] == r['var_max']:
        lab = r['bfm_1x']
        lab_sum = '= ' + fv(r['bfm_1x'])
    elif r['min_1x'] <= r['var_min'] and r['max_1x'] >= r['var_max']:
        lab = r['min_1x']
        lab_sum = '= ' + fv(r['min_1x']) + '-' + fv(r['max_1x'])
    elif r['min_1x'] > r['var_min'] and r['max_1x'] >= r['var_max']:
        lab = r['min_1x']
        lab_sum = '> ' + fv(r['bfm']) + ' -/- ' + fv(r['errmin_1x'])
    elif r['min_1x'] <= r['var_min'] and r['max_1x'] < r['var_max']:
        lab = r['max_1x']
        lab_sum = '< ' + fv(r['bfm']) + ' +/+ ' + fv(r['errmax_1x'])
    else:
        lab = r['bfm_1x']
        val1 = var_x[np.nonzero(var_x <= r['bfm'])[0].max()]
        val2 = var_x[np.nonzero(var_x >= r['bfm'])[0].min()]
        form_err = 0.5 * (val2 - val1)
        if form_err > r['err_1x']:
            r['err_1x'] = f32(form_err)
        lab_sum = '= ' + fv(r['bfm_1x']) + ' +/- ' + fv(r['err_1x'])


    # #################### OUTPUT ####################
    return dict(title=r['title'], lab=lab, lab_sum=lab_sum, lab_bfm=r['bfm'], bfm_text=fv(r['bfm']))




# =======================================================================================
# ============================== CHI2 DISTRIBUTIONS OF ALL PARAMETERS ==============================
# =======================================================================================

# ==============================
def _plot_dist(
               sol      : np.ndarray,
               obs      : str,
               obs_run  : str,
               info     : dict,
               dat_line : np.ndarray,
               Mv       : float | None = None,
               Mk       : float | None = None,
               kud      : bool = False,
               vinf     : float | None = None,
               fw2      : str | None = None,
               fw3      : str | None = None,
               spc      : str = '',
               root     : str | Path = '.'
    ) -> tuple[plt.Figure, dict, np.ndarray]:
    """
    chi2 distributions of all the parameters, final results and HR diagrams
    (translation of gbat_plot_dist.pro).

    -------
    Returns
    -------
        fig : plt.Figure
            Page with the chi2 distributions (second page of the final PDF).

        results : dict
            {parameter: output of _write_resul_one, plus the resul structure}
            for Teff, logg, yhe, micro, Lsp, q, beta, r, l, m and mdot.

        sol : np.ndarray
            Models with R, L, M and Mdot (if Mv/Mk and vinf are given).

    -----
    Notes
    -----
        - Saved files: xdr_dist/<obs_run>/info_dist_<obs><labels>.idl and
          xdr_results/<obs_run>/<obs>_<param><labels>.idl.
        - log L_spec = 4 log Teff - log g - 10.61.
        - As in IDL, sol.q (-log Q) is positive in the distributions and the
          sign of log Q is restored in the text (_write_resul_one).
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    from pathlib import Path
    import matplotlib.pyplot as plt
    import numpy as np


    # ########## PANORAMA ##########
    from gbat.io import _fi, _ff, _found, _idl2mpl, _save_gbat
    from gbat.plots import _PAGE, _axes, _crea_subplot, _plot_dist_panel, _plot_gva, _symbfill
    from gbat.variables import GBAT_DIST_DIR, GBAT_RESULTS_DIR, LSP_CONST


    # #################### LABELS ####################
    root = Path(root)
    lab = _labels(Mv, Mk, vinf)
    sol = sol.copy()


    # #################### R, L, M AND MDOT ####################
    if Mv or Mk:
        _rlm(sol, Mv=Mv, Mk=Mk, kud=kud)
    if vinf and (Mv or Mk):
        _massloss(sol, vinf)
    lsp = 4. * np.log10(sol['teff'] * 1e3) - sol['grav'] - LSP_CONST
    _save_gbat(root / GBAT_DIST_DIR / obs_run / ('info_dist_' + obs + lab), sol=sol, obs=obs, info=info, dat_line=dat_line)


    # #################### FIGURE ####################
    fig = plt.figure(figsize=_PAGE)
    pos = _crea_subplot([7, 5, 5, 5], 0.03, 0.08, xmin=0.03, ymin=0.07)
    lab_lsp = 'log !12' + 'L' + '!X/!12' + 'L' + '!X!9!dn!n!3'
    mask = lambda fw: (np.char.add(np.char.add(sol['model'], '_VT'), sol['micro']) == fw) if fw else None
    params = [
              ('Teff', sol['teff'], 'Teff (kK)', 0, 0.02, '(F4.1)'),
              ('logg', sol['grav'], 'log g (dex)', 1, 0.02, '(F4.2)'),
              ('yhe', 100. * sol['yhe'], 'Y(He) x 10!u2!n', 2, 0.01, '(F4.1)'),
              ('micro', sol['micro'].astype(float), '!7n!3!dt!n (kms!u-1!n)', 3, 0.02, '(F4.1)'),
              ('Lsp', lsp, lab_lsp, 7, 0.02, '(F4.2)'),
              ('q', sol['q'], '-log Q', 8, 0.02, '(F6.2)'),
              ('beta', sol['beta'], '!7b!3', 9, 0.02, '(F4.2)'),
              ('r', sol['r'], 'R (R!d!9n!3!n)', 12, 0.10, '(F4.1)'),
              ('l', sol['l'], 'log (L/L!d!9n!3!n)', 13, 0.02, '(F4.2)'),
              ('m', sol['m'], 'M (M!d!9n!3!n)', 17, 2.50, '(F4.1)'),
              ('mdot', sol['mdot'], '-log Mdot', 18, 0.15, '(F5.2)'),
             ]


    # #################### DISTRIBUTIONS ####################
    results = {}
    for key, var, title, ipos, extra, fmt in params:
        base = root / GBAT_RESULTS_DIR / obs_run / (obs + '_' + key + lab)
        dist = _plot_dist_one(var, sol, title, extra, base)
        if np.mean(np.abs(var)) > 0:
            _plot_dist_panel(_axes(fig, pos[ipos]), np.asarray(var, dtype=float), sol['chi'], dist['varx_sp'],
                             dist['vary_sp'], dist['var_x'], dist['var_y'], title, mask(fw2), mask(fw3))
        results[key] = dict(resul=dist['resul'], fmt=fmt)


    # #################### SUMMARY OF RESULTS ####################
    ax = _axes(fig, (pos[4][0], pos[8][1], pos[6][2], pos[6][3]))
    ax.axis('off')
    ax.set(xlim=[0, 1], ylim=[0, 1])
    sh, i = 0.06, 1
    text = lambda x, y, s, **kw: ax.text(x, y, s, fontsize=8.5, clip_on=False, **kw)

    spec = obs[:obs.find('_v')] if '_v' in obs else obs
    text(0, 1. - i * sh, spec + ' [' + spc + ']')
    i += 2
    if info['zeta_t']:
        text(0, 1. - i * sh, r'vsini (kms$^{-1}$) = ' + _fi(info['vsini'], 3) + r'  ,  $\Theta_{RT}$ (kms$^{-1}$) = '
             + _fi(info['zeta_t'], 3))
    else:
        text(0, 1. - i * sh, r'vsini (kms$^{-1}$) = ' + _fi(info['vsini'], 3) + r'  ,  vg (kms$^{-1}$) = '
             + _fi(info['vgauss'], 3))
    i += 1
    if Mv:
        text(0.45, 1. - i * sh, 'Mv (mag) = ' + _ff(Mv, 5, 2))
    if Mk:
        text(0.45, 1. - i * sh, 'MKs (mag) = ' + _ff(Mk, 5, 2))
    if vinf:
        text(0, 1. - i * sh, r'vinf (kms$^{-1}$) = ' + _fi(vinf, 4))
    i += 2
    text(0, 1. - i * sh, 'Grid = ' + info['grid'])
    i += 2

    for key, var, title, ipos, extra, fmt in params:
        out = _write_resul_one(_found(root / GBAT_RESULTS_DIR / obs_run / (obs + '_' + key + lab)), fmt)
        results[key].update(out)
        y = 1. - i * sh
        text(0.27, y, _idl2mpl(out['title']), ha='right')
        text(0.3, y, out['lab_sum'])
        text(0.75, y, '( ' + out['bfm_text'] + ' )')
        i += 1


    # #################### LABEL OF THE BEST MODEL ####################
    """
    - Printed as in IDL. Note that the IDL expressions
          sol(where(sol.teff-labteff) LE min(sol.teff-labteff)).teff
      select sol(0) or sol(1) (a 0/1 comparison is used as index), i.e. this
      label is not the model closest to the final values.
    """
    def quirk(field, value):
        diff = sol[field] - value
        idx = np.nonzero(diff)[0]
        first = idx[0] if idx.size else -1
        return sol[field][int(first <= diff.min())]

    fwlab1 = str(int(quirk('teff', results['Teff']['lab']) * 10.))
    fwlab2 = str(int(quirk('grav', results['logg']['lab']) * 100.))
    fwlab3 = str(int(quirk('yhe', results['yhe']['lab']) * 100.))
    fwlab4 = str(int(quirk('q', results['q']['lab']) * 100.))
    fwlab5 = int(quirk('beta', results['beta']['lab']) * 10.)
    fwlab5 = str(fwlab5) if fwlab5 >= 10 else '0' + str(fwlab5)
    print('T' + fwlab1 + 'g' + fwlab2 + 'He' + fwlab3 + 'Q' + fwlab4 + 'b' + fwlab5)


    # #################### SHRD AND HRD ####################
    labteff = results['Teff']['lab']
    ax = _axes(fig, (pos[14][0] + 0.03, pos[18][1], pos[15][2] - 0.10, pos[15][3]))
    _plot_gva(ax, diag='shrd', mli=(15., 91.), xli=(4.79, 4.31), yli=(3.01, 4.29))
    _symbfill(ax, np.log10(labteff * 1e3), results['Lsp']['lab'], 8, 1.2, 2, fill=True)
    if np.mean(sol['l']) > 0:
        ax = _axes(fig, (pos[15][0] + 0.13, pos[18][1], pos[16][2], pos[16][3]))
        _plot_gva(ax, diag='hrd', mli=(15., 91.), xli=(4.79, 4.31), yli=(4.15, 6.19))
        _symbfill(ax, np.log10(labteff * 1e3), results['l']['lab'], 8, 1.2, 2, fill=True)


    # #################### OUTPUT ####################
    return fig, results, sol
