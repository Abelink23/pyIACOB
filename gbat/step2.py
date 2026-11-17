from __future__ import annotations

# ======================================================================
# ============================== STEP 2 ==============================
# ======================================================================

# ==============================
def chi_step2(
              obs     : str,
              obs_run : str,
              teff    : float | None = None,
              grav    : float | None = None,
              yhe     : float | None = None,
              q       : float | None = None,
              beta    : float | None = None,
              micro   : str | None = None,
              Mv      : float | None = None,
              Mk      : float | None = None,
              vinf    : float | None = None,
              ps      : bool = False,
              kud     : bool = False,
              fw2     : str | None = None,
              fw3     : str | None = None,
              spc     : str = '',
              inst    : str | None = None,
              root    : str | Path = '.'
    ) -> dict | None:
    """
    Second step of IACOB-GBAT: global chi2 distributions, central values and
    uncertainties of the parameters, and summary figures (translation of
    gbat_chi_step2.pro).

    Once the chi2 of each line is computed for the whole subgrid (step 1,
    chi_global), this step offers the possibility

        [1] to fix some of the stellar parameters to a value of the subgrid
            (teff, grav, yhe, q, beta, micro), and
        [2] to exclude some of the lines (column step2_in of
            gbat_linesHHe_<obs_run>.dat).

    ----------
    Parameters
    ----------
        obs : str
            Name of the results of step 1 (fileout of chi_global), e.g.
            'HD36512_20111109_032351_M_V85000_vr.dat_v13_m33', read from
            xdr_files/<obs_run>/<obs>.idl (or the IDL .xdr).

        obs_run : str
            Observing run.

        teff, grav, yhe, q, beta : float | None, optional
            Fix Teff [K], log g, Y(He), -log Q or beta to a value of the
            subgrid. None (or 0) -> not fixed.

        micro : str | None, optional
            Fix the microturbulence (e.g. '005').

        Mv, Mk : float | None, optional
            Absolute V or Ks magnitude (only one): R, log L and M are computed.

        vinf : float | None, optional
            Terminal velocity [km/s]: Mdot is computed (with Mv or Mk).

        ps : bool, optional
            Open the final PDF with PS_VIEWER.

        kud : bool, optional
            R computed with Kudritzki's formula (from vsyn).

        fw2, fw3 : str | None, optional
            Two models ('<model>_VT<micro>') marked in the chi2 distributions.

        spc : str, optional
            Spectral classification (label).

        inst : str | None
            Instrument (see INSTRUMENTS in gbat.variables).

        root : str | Path, optional
            Working folder. Default '.'.

    -------
    Returns
    -------
        results : dict | None
            - 'sol' : structured array with the models sorted by chi2,
            - 'dat_line' : lines used in the analysis,
            - 'sigma_real' : noise of each line from the best-fitting model,
            - 'results' : dict {parameter: dict(resul, lab, lab_sum, lab_bfm)},
            - 'bfm' : best-fitting model ('<model>VT<micro>'),
            - 'figure' : PDF summary figure.

    -----
    Notes
    -----
        - Files: xdr_dist/<obs_run>/info_dist_<obs>[_Mv..][_vinf..].idl,
          xdr_results/<obs_run>/<obs>_<param>[_Mv..][_vinf..].idl,
          output_figures/<obs_run>/<obs>[_Mv..][_vinf..].pdf (page 1: line profiles,
          page 2: chi2 distributions) and spectra_bfm/<obs_run>/<obs>_fw.asc.
        - As in IDL, the values of casa_red, resol and res0 used for the figures
          and the synthetic spectrum are those of the current configuration
          (chi_input is called again in extract_lines), while 'info' keeps the
          values saved in step 1.
        - If a FASTWIND model is not available (e.g. step 2 run without access
          to the grid), the synthetic profiles and spectrum are skipped with a
          warning.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    from pathlib import Path
    import shutil
    import subprocess
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages


    # ########## PANORAMA ##########
    from broadening.convol import _convol, _crea_sp_fw
    from gbat.dist import _labels, _plot_dist
    from gbat.grid import _load_model
    from gbat.input import _chi_input, _lineshhe, _question, _read_lines_hhe
    from gbat.io import _fi, _found, _restore_gbat, _save_gbat, _wascii
    from gbat.plots import _PAGE, _axes, _crea_subplot, _plot_prof
    from gbat.variables import GBAT_BFM_DIR, GBAT_DIST_DIR, GBAT_FIGURES_DIR, GBAT_FILES_DIR


    # #################### INPUT ####################
    root = Path(root)
    cfg = _chi_input(obs_run, inst, root)
    casa_red0 = cfg.casa_red
    model2plot = 0 # The model whose profiles will be plotted

    print('')
    print(' NOW I PROCEED WITH STEP 2 USING')
    print('')
    print('   - ' + obs)

    if Mv and Mk:
        _question('Both Mv and Mk specified, please use only one')
        return None


    # #################### RESTORE STEP 1 ####################
    d = _restore_gbat(_found(root / GBAT_FILES_DIR / obs_run / obs))
    sol = d['sol'][np.isfinite(d['sol']['chi'])]
    casa_red = d['casa_red']
    vsini, zeta_t, vgauss = d['vsini'], d['zeta_t'], d['vgauss']
    wo0, fo0, pe0 = d['wo0'], d['fo0'], d['pe0']

    tt = casa_red.find('iac_gbatgrid_HHe')
    info = dict(casa=casa_red[:tt], grid=casa_red[tt:tt + 21], res0=d['res0'], resol=d['resol'],
                vsini=vsini, zeta_t=zeta_t, vgauss=vgauss)


    # #################### FIXED PARAMETERS ####################
    sol = _extract_models(sol, teff=teff, grav=grav, yhe=yhe, q=q, beta=beta, micro=micro)


    # #################### LINES OF STEP 2 ####################
    """
    - Lines to be considered/excluded must be indicated in gbat_linesHHe_<obs_run>.dat
    """
    sol, dat_line = _extract_lines(sol, d['dat_line'], root)
    cfg = _chi_input(obs_run, inst, root) # IDL: extract_lines calls chi_input (COMMON variables reset)


    # #################### GLOBAL CHI2 ####################
    sol, sigma_real = _chi_total(sol, dat_line, coef_sigma=1.)
    _save_gbat(root / GBAT_DIST_DIR / obs_run / ('info_dist_' + obs), sol=sol, obs=obs, info=info, dat_line=dat_line)


    # #################### FIGURES ####################
    matplotlib.rcParams.update({'font.size': 8, 'mathtext.default': 'regular'})
    lab = _labels(Mv, Mk, vinf)

    # ########## CHI2 DISTRIBUTIONS ##########
    fig_dist, results, sol = _plot_dist(sol, obs, obs_run, info, dat_line, Mv=Mv, Mk=Mk, kud=kud, vinf=vinf,
                                        fw2=fw2, fw3=fw3, spc=spc, root=root)

    # ########## LINE PROFILES ##########
    """
    - The lines are placed following the column rowx10+column_pos of
      gbat_linesHHe_<obs_run>.dat. The last row holds the information of the
      run.
    """
    dat = _read_lines_hhe(root)
    rc_pos = dat['rc_pos']
    pos_num = np.sort(rc_pos[rc_pos > 0])
    pos_vec = np.zeros(int(rc_pos.max() / 10.) + 1, dtype=int)
    for i in range(pos_vec.size - 1):
        mm = (rc_pos >= (i + 1) * 10) & (rc_pos < (i + 2) * 10)
        pos_vec[i] = max(rc_pos[mm] - (i + 1) * 10) if mm.any() else 0
    pos_vec[-1] = 1
    pos = _crea_subplot(list(pos_vec), 0.04, 0.10, xmin=0.05, xmax=0.95, ymin=0.05, ymax=0.95)

    fig_prof = plt.figure(figsize=_PAGE)
    best = sol[model2plot]
    model_file = _model_file_safe(casa_red0 + best['red'] + best['model'])
    model = _load_model(model_file) if model_file is not None else None
    if model is None:
        print('WARNING ... ' + casa_red0 + best['red'] + best['model'] + '.xdr NOT FOUND! (synthetic profiles not plotted)')
    for i in range(dat_line.size):
        line = dat_line['line'][i]
        mm0 = np.nonzero(dat['fwlab'] == line)[0][0]
        mm1 = int(np.nonzero(pos_num == rc_pos[mm0])[0][0])
        w0 = f0 = None
        if model is not None and (line, best['micro']) in model['lines']:
            x, y = model['lines'][(line, best['micro'])]
            w0, f0 = _convol(x, y, cfg.res0, resol=cfg.resol, vsini=vsini, zeta_t=zeta_t, vgauss=vgauss)
        _plot_prof(_axes(fig_prof, pos[mm1]), _lineshhe(line, root), dat_line['sigma'][i], sigma_real[i],
                   w0, f0, wo0, fo0, pe0)

    ax = _axes(fig_prof, pos[int(pos_vec.sum()) - 1])
    ax.axis('off')
    ax.set(xlim=[0, 1], ylim=[0, 1])
    ax.text(0, 0.90, cfg.casa_red, fontsize=8)
    ax.text(0, 0.65, obs + ' [' + spc + ']' + '   R=' + _fi(cfg.resol, 5), fontsize=8)
    ax.text(0, 0.40, best['model'] + '_VT' + best['micro'], fontsize=8)

    # ########## PDF ##########
    figure = root / GBAT_FIGURES_DIR / obs_run / (obs + lab + '.pdf')
    with PdfPages(figure) as pdf:
        pdf.savefig(fig_prof)
        pdf.savefig(fig_dist)
    plt.close(fig_prof)
    plt.close(fig_dist)
    if ps:
        if shutil.which(cfg.ps_viewer):
            subprocess.Popen([cfg.ps_viewer, str(figure)])
        else:
            print(f"WARNING ... PDF viewer '{cfg.ps_viewer}' not found (see PS_VIEWER in gbat.variables)")


    # #################### SYNTHETIC SPECTRUM OF THE BFM ####################
    model_file = _model_file_safe(cfg.casa_red + best['red'] + best['model'])
    if model_file is None:
        print('WARNING ... ' + cfg.casa_red + best['red'] + best['model'] + '.xdr NOT FOUND! (synthetic spectrum not created)')
    else:
        yfw = _crea_sp_fw(wo0, _load_model(model_file), cfg.res0, cfg.resol, vsini, zeta_t=zeta_t, vgauss=vgauss,
                          micro=best['micro'])
        (root / GBAT_BFM_DIR / obs_run).mkdir(parents=True, exist_ok=True)
        _wascii(root / GBAT_BFM_DIR / obs_run / (obs + '_fw.asc'), wo0, yfw)


    # #################### OUTPUT ####################
    return dict(sol=sol, dat_line=dat_line, sigma_real=sigma_real, results=results,
                bfm=best['model'] + 'VT' + best['micro'], figure=figure)




# ==============================
def _model_file(*args):
    from gbat.grid import _model_file
    return _model_file(*args)




# ==============================
def _model_file_safe(base : str) -> Path | None:
    """
    File of the best-fitting model for the figures and the synthetic spectrum
    of step 2. A connection problem with the grid share gives a warning
    (None) instead of an error, since the results are already saved.
    """

    # #################### FILE ####################
    try:
        return _model_file(base)
    except ConnectionError as e:
        print(f'WARNING ... grid share not available ({e})')
        return None




# ===========================================================================
# ============================== SELECTION ==============================
# ===========================================================================

# ==============================
def _extract_models(
                    sol   : np.ndarray,
                    teff  : float | None = None,
                    grav  : float | None = None,
                    yhe   : float | None = None,
                    q     : float | None = None,
                    beta  : float | None = None,
                    micro : str | None = None
    ) -> np.ndarray:
    """
    Fix some parameters to a value of the subgrid (translation of
    gbat_extract_models.pro).

    -----
    Notes
    -----
        - Tolerances: Teff 0.1 K, log g 0.05, Y(He) 0.01, -log Q 0.01 (q is
          compared with -sol.q, as in IDL, i.e. log Q must be given with its
          sign), beta 0.05. Zero values (or None) mean not fixed.
    """

    # #################### MODULES & PACKAGES ####################
    from gbat.io import _fstr


    # #################### SELECTION ####################
    if teff:
        sol = sol[abs(1e3 * sol['teff'] - teff) < 0.1]
    if grav:
        sol = sol[abs(sol['grav'] - grav) < 0.05]
    if yhe:
        sol = sol[abs(sol['yhe'] - yhe) < 0.01]
    if q:
        sol = sol[abs(-1. * sol['q'] - q) < 0.01]
    if beta:
        sol = sol[abs(sol['beta'] - beta) < 0.05]
    if micro:
        sol = sol[sol['micro'] == micro]


    # #################### MESSAGES ####################
    fixed = any(v for v in (teff, grav, yhe, q, beta, micro))
    if fixed:
        print('')
    for lab, v in (('Teff  fixed to ', teff), ('logg  fixed to ', grav), ('YHe   fixed to ', yhe),
                   ('logQ  fixed to ', q), ('beta  fixed to ', beta)):
        if v:
            print(lab + _fstr(v))
    if micro:
        print('micro fixed to ' + str(micro))
    if fixed:
        print('')
    if sol.size == 0:
        raise ValueError('No models left after fixing the parameters (they must be values of the subgrid).')


    # #################### OUTPUT ####################
    return sol




# ==============================
def _extract_lines(
                   sol      : np.ndarray,
                   dat_line : np.ndarray,
                   root     : str | Path = '.'
    ) -> tuple[np.ndarray, np.ndarray]:
    """
    Exclude the lines with step2_in = 0 from the final analysis (translation
    of gbat_extract_lines.pro).

    -------
    Returns
    -------
        sol : np.ndarray
            New sol with only the chi2 of the lines with step2_in = 1, and chi
            = sum of their chi2.

        dat_line : np.ndarray
            Information of these lines.

    -----
    Notes
    -----
        - A line with step2_in = 1 that was not computed in step 1 (chi_in = 0)
          is NOT a valid option: as in IDL, its entry stays empty.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # ########## PANORAMA ##########
    from gbat.input import _read_lines_hhe
    from gbat.io import _dat_line_init, _fa, _sol_init


    # #################### LINES ####################
    dat = _read_lines_hhe(root)
    lines = list(dat['fwlab'][dat['step2_in'] == 1])
    nlines, nmodels = len(lines), sol.size
    dat_line0, sol0 = dat_line, sol


    # #################### NEW STRUCTURES ####################
    dat_line = _dat_line_init(nlines)
    sol = _sol_init(nlines, nmodels)
    for name in sol.dtype.names:
        if name != 'chi_lin':
            sol[name] = sol0[name]

    print('')
    print(' EXCLUDED LINES IN STEP 2')
    print('')
    for line in dat_line0['line']:
        if line not in lines:
            print(_fa(line, 10) + '  ' + _fa('X', 3))
    print('')

    for i, line in enumerate(lines):
        tt1 = np.nonzero(dat_line0['line'] == line)[0]
        if tt1.size:
            dat_line[i] = dat_line0[tt1[0]]
            sol['chi_lin'][:, i] = sol0['chi_lin'][:, tt1[0]]
        else:
            print(f'WARNING ... {line} has step2_in = 1 but was not computed in step 1 (chi_in = 0): NOT a valid option')
    sol['chi'] = sol['chi_lin'].sum(axis=1)


    # #################### OUTPUT ####################
    return sol, dat_line




# ======================================================================
# ============================== GLOBAL CHI2 ==============================
# ======================================================================

# ==============================
def _chi_total(
               sol        : np.ndarray,
               dat_line   : np.ndarray,
               option     : int = 2,
               coef_sigma : float = 1.
    ) -> tuple[np.ndarray, np.ndarray]:
    """
    Global chi2 of each model and real noise of each line (translation of
    gbat_chi_total.pro).

    ----------
    Parameters
    ----------
        sol : np.ndarray
            Models with the chi2_tilda = sum((F_obs - F_mod)^2) of each line.

        dat_line : np.ndarray
            Information of the lines (sigma from the SNR, np).

        option : int, optional
            Global chi2: 1 (all the points joined together), 2 (sum of the
            reduced chi2 of the lines, default) or 3 (master quantities, for
            the 14 lines of the benchmark; to be improved).

        coef_sigma : float, optional
            Factor applied to the real noise of the lines. Default 1.

    -------
    Returns
    -------
        sol : np.ndarray
            Models sorted by the global chi2, with chi_lin replaced by the
            reduced chi2 of each line.

        sigma_real : np.ndarray
            Real noise of each line.

    -----
    Notes
    -----
        - The real noise of each line is computed iteratively (4 loops):
          first from the best model of each line,
            sigma = sqrt(min(chi2_tilda_line) / NF),
          then from the best global model of the previous loop.
        - The reduced chi2 printed for each model in the final table are
          those of the last loop, computed before the last sorting (as in
          IDL).
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # ########## PANORAMA ##########
    from gbat.io import _fa, _ff, _fi


    # #################### FORMATS ####################
    row_a = lambda lab, w, values, extra='': _fa(lab, w) + values + extra
    fa8 = lambda v: ''.join(_fa(s, 8) + ' ' for s in v)
    fi3 = lambda v, sp: ''.join(_fi(x, 3) + ' ' * 6 for x in v)
    ff74 = lambda v: ''.join(_ff(x, 7, 4) + '  ' for x in v)


    # #################### ZERO APPROACH ####################
    nmodels, nlines = sol.size, dat_line.size
    npts = dat_line['np'].astype(float)
    sol = sol.copy()
    sol['chi'] = sol['chi_lin'].sum(axis=1)
    sol = sol[np.argsort(sol['chi'], kind='stable')]
    name = lambda s: s['model'] + 'VT' + s['micro']


    # #################### BEST MODEL OF EACH LINE ####################
    print(' ')
    print('      LINE :          MODEL            :  sigma0  -- sigma(SNR)')
    print(' --------- :  ------------------------ : --------    ----------')
    for mm in range(nlines):
        cmin = sol['chi_lin'][:, mm].min()
        jj = np.nonzero(sol['chi_lin'][:, mm] <= cmin)[0][0]
        print(_fa(dat_line['line'][mm], 10) + _fa(' : ', 3) + _fa(name(sol[jj]), 25) + _fa(' : ', 3)
              + _ff(np.sqrt(cmin / npts[mm]), 7, 5) + _fa('  --  ', 6) + _ff(dat_line['sigma'][mm], 7, 5))


    # #################### REAL SIGMAS (ITERATIVE) ####################
    total_lab = 'Total-' + str(nlines)
    print(' ')
    print(row_a('Line        :', 16, '   ' + fa8(dat_line['line']), _fa(total_lab, 8)))
    print(row_a('SNR         :', 16, '      ' + fi3(1. / dat_line['sigma'], 6)))
    print(row_a('sigma (SNR) :', 16, '    ' + ff74(dat_line['sigma'])))

    sigma_real = np.zeros(nlines)
    chi2_red_lin = None
    for loop in range(4):
        if loop == 0: # Option 0: chi2_tilda_line (best model for each line) / NF
            sigma_real = coef_sigma * np.sqrt(sol['chi_lin'].min(axis=0) / npts)
        else: # Option 1: chi2_tilda_line (best global model) / NF
            sigma_real = coef_sigma * np.sqrt(sol['chi_lin'][0] / npts)
        print(row_a('sigma_real (' + _fi(loop, 2) + '):', 16, '    ' + ff74(sigma_real), _fa(name(sol[0]), 25)))

        # ########## REDUCED CHI2 PER LINE ##########
        chi2_red_lin = sol['chi_lin'] / sigma_real**2 / npts
        print(row_a('chi2_red (' + _fi(loop, 2) + '):', 16, '    ' + ff74(chi2_red_lin[0])))

        # ########## GLOBAL CHI2 ##########
        if option == 1: # Everything joined together
            sol['chi'] = sol['chi_lin'].sum(axis=1) / npts.sum()
        elif option == 2: # Line per line
            sol['chi'] = chi2_red_lin.sum(axis=1)
        elif option == 3: # Master quantities ** TO BE IMPROVED ** (fixed line indices of the benchmark)
            c = sol['chi_lin']
            mq = [c[:, 0] + c[:, 13], c[:, 1:4].sum(axis=1), c[:, 5] + c[:, 6] + c[:, 11] + c[:, 12],
                  c[:, 6], c[:, 12], c[:, 1]]
            nq = [npts[0] + npts[13], npts[1:4].sum(), npts[5] + npts[6] + npts[11] + npts[12],
                  npts[6], npts[12], npts[1]]
            sol['chi'] = sum(m / n for m, n in zip(mq, nq))
        else:
            raise ValueError(f"Option 'option' must be 1, 2 or 3. Got {option}.")

        # ########## SORTING ##########
        sol = sol[np.argsort(sol['chi'], kind='stable')]


    # #################### FINAL REDUCED CHI2 PER LINE ####################
    sol['chi_lin'] = sol['chi_lin'] / sigma_real**2 / npts
    print(row_a('chi2_red_fin :', 16, '    ' + ff74(sol['chi_lin'][0])))


    # #################### SUMMARY ####################
    print('     ')
    print(row_a('Line        :', 13, '    ' + fa8(dat_line['line']), _fa(total_lab, 8)))
    print(row_a('SNR         :', 13, '      ' + fi3(1. / dat_line['sigma'], 6)))
    print(row_a('sigma (SNR) :', 13, '    ' + ff74(dat_line['sigma'])))
    print(row_a('sigma_real  :', 13, '    ' + ff74(sigma_real)))
    print(' ')
    nf = [_fi(x, 4) + ' ' * 5 for x in npts] # IDL format (A13,4X,15(I4,5X)): 15 values per row
    print(_fa('NF          :', 13) + '    ' + ''.join(nf[:15]))
    for k in range(15, len(nf), 15):
        print(''.join(nf[k:k + 15]))

    ss, cmin = 0, sol['chi'].min()
    while True:
        print(_fa('Model (', 7) + _fi(ss, 3) + _fa(') :', 3) + '    '
              + ''.join(_ff(v, 6, 2) + '   ' for v in chi2_red_lin[ss]) + _ff(sol['chi'][ss], 6, 2) + '  '
              + _fa(name(sol[ss]), 25))
        ss += 1
        if ss >= nmodels or not sol['chi'][ss] <= cmin + 1.5:
            break
    print(' ')


    # #################### OUTPUT ####################
    return sol, sigma_real
