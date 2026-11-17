from __future__ import annotations

# ==========================================================================
# ============================== SUMMARY TABLE ==============================
# ==========================================================================

# ==============================
def summarytable(
                 obsrun : str,
                 rlm    : bool = False,
                 root   : str | Path = '.'
    ) -> Path:
    """
    Table with the results of all the stars of an observing run
    (translation of gbat_summarytable.pro).

    ----------
    Parameters
    ----------
        obsrun : str
            Observing run (sub-folder of xdr_results).

        rlm : bool, optional
            Include R, log L, M and log Mdot. Default False.

        root : str | Path, optional
            Working folder. Default '.'.

    -------
    Returns
    -------
        file : Path
            summary/<obsrun>_summary.txt, with the IDL format
            (A70,4X,n(A20,2X),A30): star, Teff, log g, log Q, Y(He), micro, beta,
            [R, log L, M, Mdot,] and the label of the best-fitting model.

    -----
    Notes
    -----
        - The stars are those with a file xdr_results/<obsrun>/*_Teff*.idl
          (Python) or *_Teff*.xdr (IDL version).
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    from pathlib import Path


    # ########## PANORAMA ##########
    from gbat.dist import _write_resul_one
    from gbat.io import _fa, _fi, _found
    from gbat.variables import GBAT_RESULTS_DIR, GBAT_SUMMARY_DIR


    # #################### FILES ####################
    root = Path(root)
    folder = root / GBAT_RESULTS_DIR / obsrun
    names = sorted({f.name[:-4] for ext in ('.idl', '.xdr') for f in folder.glob('*_Teff*' + ext)})
    print('# STARS IN ' + obsrun + ' : ', len(names))
    (root / GBAT_SUMMARY_DIR).mkdir(parents=True, exist_ok=True)
    out = root / GBAT_SUMMARY_DIR / (obsrun + '_summary.txt')


    # #################### HEADER ####################
    cols = ['Teff(kK)', 'logg(dex)', 'logQ', 'Y(He)', 'micro (km/s)', 'beta'] + (['R', 'logL', 'M', 'Mdot'] if rlm else [])
    row = lambda first, values, last: _fa(first, 70) + '    ' + ''.join(_fa(v, 20) + '  ' for v in values) + _fa(last, 30)
    lines = [row('Star-ID', cols, 'BFM'),
             row('-' * 65, ['------------------'] * len(cols), '-' * 30)]


    # #################### STARS ####################
    for name in names:
        tt2 = name.find('Teff')
        file = name[:tt2 - 1]
        extra_lab = name[tt2 + 4:]
        res = lambda key, fmt: _write_resul_one(_found(folder / (file + '_' + key + extra_lab)), fmt)

        teff = res('Teff', '(F4.1)')
        logg = res('logg', '(F4.2)')
        q = res('q', '(F6.2)')
        yhe = res('yhe', '(F4.1)')
        micro = res('micro', '(I2)')
        beta = res('beta', '(F4.1)')

        lab_bfm = ('T' + _fi(teff['lab_bfm'] * 10., 3).strip()
                   + 'g' + _fi(logg['lab_bfm'] * 100., 3).strip()
                   + 'He' + ('0' + _fi(yhe['lab_bfm'], 2).strip())[-2:]
                   + 'Q' + _fi(abs(q['lab_bfm']) * 100., 4).strip()
                   + 'b' + ('0' + _fi(beta['lab_bfm'] * 10., 2).strip())[-2:]
                   + '_m' + ('0' + _fi(micro['lab_bfm'], 2).strip())[-2:])

        values = [teff, logg, q, yhe, micro, beta]
        if rlm:
            values += [res('r', '(F4.1)'), res('l', '(F4.1)'), res('m', '(F4.1)'), res('mdot', '(F4.1)')]
        lines.append(row(file + extra_lab, [v['lab_sum'] for v in values], lab_bfm))


    # #################### WRITE ####################
    out.write_text('\n'.join(lines) + '\n')
    print('')
    print('I HAVE GENERATED THE FILE ' + str(out) + '  WITH THE SUMMARY OF RESULTS.')
    print('I HOPE YOU WRITE A NICE AND IMPORTANT PAPER WITH IT! ENJOY!')
    print('')


    # #################### OUTPUT ####################
    return out




# ==========================================================================
# ============================== LAUNCHER TABLE ==============================
# ==========================================================================

# ==============================
def _crea_par(
              apar      : str,
              par       : float,
              epar      : float,
              step_par  : float,
              par_range : tuple[float, float]
    ) -> tuple[float, float]:
    """
    New range of a parameter from its result (crea_par of
    gbat_crea_launcher.pro): par -/+ max(epar, step_par), or up to the edge
    of par_range for lower ('>') and upper ('<') limits.
    """

    # #################### RANGE ####################
    delta = max(epar, step_par)
    if apar == '=':
        return par - delta, par + delta
    if apar == '>':
        return par - delta, max(par_range)
    if apar == '<':
        return min(par_range), par + delta


    # #################### OUTPUT ####################
    raise ValueError(f"Unknown result type '{apar}' (expected '=', '>' or '<').")




# ==============================
def crea_launcher(
                  file      : str,
                  lab_extra : str = '',
                  inst      : str = '',
                  root      : str | Path = '.'
    ) -> Path:
    """
    Launcher table for a new run of IACOB-GBAT from the summary of a previous
    one (translation of gbat_crea_launcher.pro).

    ----------
    Parameters
    ----------
        file : str
            Observing run: reads summary/<file>_summary.txt and the step-1
            results in xdr_files/<file>.

        lab_extra : str, optional
            Extra label of the new table. Default ''.

        inst : str, optional
            Instrument used when it cannot be identified from the file name
            ('_F_' feros, '_M_' hermes, '_N_' fies, 'BLOeM' bloems). The
            resolution label (two digits after '_V') is appended.

        root : str | Path, optional
            Working folder. Default '.'.

    -------
    Returns
    -------
        table : Path
            tables_launcher/table_<file>[_<lab_extra>].dat.

    -----
    Notes
    -----
        - New ranges: Teff +/- max(2 err, 2100 K) (+/- 3100 K for limits),
          log g +/- max(2 err, 0.31), log Q +/- max(2 err, 0.41), Y(He) +/-
          max(2 err, 0.026), micro +/- max(1.3 err, 4.1 or 10.1 km/s);
          beta = [0.7, 1.4].
        - As in IDL, the columns fw_ref and SPC are fixed to
          'T370g390He10Q1350b10' and 'O6.5V', and vmac is set to 0 for
          vsini >= 200 km/s.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    from pathlib import Path


    # ########## PANORAMA ##########
    from gbat.io import _fa, _ff, _fi, _found, _idl_round, _readcol, _restore_gbat
    from gbat.variables import GBAT_FILES_DIR, GBAT_LAUNCHER_DIR, GBAT_SUMMARY_DIR


    # #################### SUMMARY ####################
    root = maindir / 'gbat'
    lab_extra = '_' + lab_extra if lab_extra else ''
    (sp, ateff, teff, _, eteff, agrav, grav, _, egrav, alq, lq, _, elq, ayhe, yhe, _, eyhe,
     amic, mic, _, emic) = _readcol(root / GBAT_SUMMARY_DIR / (file + '_summary.txt'),
                                    fmt='(A,A,F,A,F,A,F,A,F,A,F,A,F,A,F,A,F,A,F,A,F)')


    # #################### TABLE ####################
    out = root / GBAT_LAUNCHER_DIR / ('table_' + file + lab_extra + '.dat')
    lines = ['                                              Name             run            inst vsini  vmac   Teffm   '
             'TeffM  loggm  loggM   Qm   QM   bm   bM   Hem   HeM micm micM               fw_ref             SPC    Mv vinf',
             '-------------------------------------------------- --------------- --------------- ----- ----- ------- ------- '
             '------ ------ ---- ---- ---- ---- ----- ----- ---- ---- -------------------- --------------- ----- ----']

    for i, s in enumerate(sp):
        # ########## SPECTRUM ##########
        sp_in = s
        for ext in ('.dat', '.ascii', '.asc'):
            k = s.find(ext)
            if k >= 0:
                sp_in = s[:k + len(ext)]
                break
        k = s.find('_V')
        resol0 = s[k + 2:k + 4] if k >= 0 else ''
        k1, k2 = s.find('_Mv'), s.find('_vinf')
        name = s[:k1] if k1 >= 0 else (s[:k2] if k2 >= 0 else s)

        # ########## STEP 1 ##########
        d = _restore_gbat(_found(root / GBAT_FILES_DIR / file / name))
        vs, vm = float(d['vsini']), float(d['zeta_t'])
        if vs >= 200.:
            vm = 0.

        # ########## MV, VINF ##########
        mv = vinf = 0.
        if k1 >= 0:
            rest = s[k1 + 3:]
            mv = -1e-2 * float(rest.split('_')[0])
        if k2 >= 0:
            rest = s[k2 + 5:]
            vinf = float(rest.split('_')[0])

        # ########## INSTRUMENT ##########
        inst_i = inst
        for key, val in (('_F_', 'feros'), ('_M_', 'hermes'), ('_N_', 'fies'), ('BLOeM', 'bloems')):
            if key in s:
                inst_i = val

        # ########## RANGES ##########
        teff_min, teff_max = _crea_par(ateff[i], _idl_round(teff[i]) * 1e3, 2. * eteff[i] * 1e3, 2100.,
                                       (teff[i] * 1e3 - 3100., teff[i] * 1e3 + 3100.))
        grav_min, grav_max = _crea_par(agrav[i], grav[i], 2. * egrav[i], 0.31, (grav[i] - 0.31, grav[i] + 0.31))
        lq_min, lq_max = _crea_par(alq[i], lq[i], 2. * elq[i], 0.41, (-14.1, -11.7))
        lq_min, lq_max = -abs(lq_min), -abs(lq_max)
        if lq_min <= -13.8:
            lq_min = -14.1
        if lq_max <= -13.0:
            lq_max = -12.9
        yhe_min, yhe_max = _crea_par(ayhe[i], yhe[i], 2. * eyhe[i], 2.6, (5, 31))
        mic_min, mic_max = _crea_par(amic[i], mic[i], 1.3 * emic[i], 4.1 if mic[i] <= 11 else 10.1, (0., 30.))
        for lo, hi in ((25, 30), (20, 25), (15, 20), (11, 15)):
            if lo <= mic_max <= hi:
                mic_max = float(hi)

        # ########## LINE ##########
        lines.append(_fa(sp_in, 50) + ' ' + _fa(file, 15) + ' ' + _fa(inst_i + resol0, 15)
                     + '   ' + _fi(vs, 3) + '   ' + _fi(vm, 3)
                     + '   ' + _fi(teff_min, 5) + '   ' + _fi(teff_max, 5)
                     + '   ' + _ff(grav_min, 4, 2) + '   ' + _ff(grav_max, 4, 2)
                     + ' ' + _ff(-lq_max, 4, 1) + ' ' + _ff(-lq_min, 4, 1)
                     + '  ' + _ff(0.7, 3, 1) + '  ' + _ff(1.4, 3, 1)
                     + '  ' + _ff(_idl_round(yhe_min) * 0.01, 4, 2) + '  ' + _ff(_idl_round(yhe_max) * 0.01, 4, 2)
                     + '   ' + _fi(max(mic_min, 0.), 2) + '   ' + _fi(min(mic_max, 30.), 2)
                     + ' ' + _fa('T370g390He10Q1350b10', 20) + ' ' + _fa('O6.5V', 15)
                     + ' ' + _ff(mv, 5, 2) + ' ' + _fi(vinf, 4))


    # #################### WRITE ####################
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text('\n'.join(lines) + '\n')


    # #################### OUTPUT ####################
    return out




# ===================================================================
# ============================== COPYCUT ==============================
# ===================================================================

# ==============================
def copycut(
            obs_run   : str,
            spectrum  : str,
            spec_comp : str,
            vrad      : bool = False,
            xr        : tuple[float, float] = (5850., 5900.),
            root      : str | Path = '.'
    ) -> Path:
    """
    Copy the masks (line and SNR regions) of a pre-processed spectrum to
    another spectrum of the same star (translation of gbat_copycut.pro).

    ----------
    Parameters
    ----------
        obs_run : str
            Observing run.

        spectrum : str
            Spectrum in the data folder defined in paths.txt.

        spec_comp : str
            Reference spectrum in spectra_cut/<obs_run> (4 columns).

        vrad : bool, optional
            Shift spec_in in wavelength by the difference of two points marked
            with the mouse (first on spec_in, then on the reference). Default
            False.

        xr : tuple[float, float], optional
            Wavelength range of the plot [A]. Default (5850, 5900).

        root : str | Path, optional
            Working folder. Default '.'.

    -------
    Returns
    -------
        file : Path
            spectra_cut/<obs_run>/<spec_in> (4 columns).

    -----
    Notes
    -----
        - Both spectra are cut to their common range and the masks of the
          reference are copied point by point (both spectra should have the
          same sampling), as in IDL.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    from pathlib import Path
    import matplotlib.pyplot as plt
    import numpy as np
    from spec import spec

    # ########## PANORAMA ##########
    from gbat.plots import _xhair
    from gbat.io import _read_spectrum, _write_spectrum_cut
    from gbat.plots import _color
    from gbat.variables import GBAT_CUT_DIR


    # #################### SPECTRA ####################
    root = Path(root)
    wo0, fo0, pe0, sig0 = _read_spectrum(root / GBAT_CUT_DIR / obs_run / spec_comp, 4)
    s = spec(spectrum)
    w1, f1 = s.wave, s.flux


    # #################### PLOT ####################
    fig = plt.figure(num='IACOB-GBAT copycut', figsize=(6., 4.))
    fig.clf()
    ax = fig.add_subplot(111)
    ax.plot(w1, f1, color='black', lw=0.8)
    ax.plot(wo0, fo0, color=_color(3), lw=0.8)
    ax.set(xlim=xr, ylim=[0.4, 1.1], title='First black (PROB), then blue (REF)')
    plt.pause(0.05)


    # #################### RADIAL VELOCITY ####################
    if vrad:
        (x1, _), (x2, _) = _xhair(ax, 2)
        w1 = w1 + (x2 - x1)
        ax.plot(w1, f1, color=_color(4), lw=0.8)


    # #################### COMMON RANGE ####################
    mm = (w1 >= wo0.min()) & (w1 <= wo0.max())
    w1, f1 = w1[mm], f1[mm]
    mm = (wo0 >= w1.min()) & (wo0 <= w1.max())
    wo0, fo0, pe0, sig0 = wo0[mm], fo0[mm], pe0[mm], sig0[mm]
    n = min(wo0.size, w1.size)
    pe1, sig1 = pe0[:n], sig0[:n]
    ax.plot(wo0[pe0 == 1], fo0[pe0 == 1] - 0.1, '.', color=_color(2), ms=2)
    ax.plot(w1[:n][pe1 == 1], f1[:n][pe1 == 1], '.', color=_color(2), ms=2)
    plt.pause(0.5)


    # #################### WRITE ####################
    out = root / GBAT_CUT_DIR / obs_run / spectrum
    _write_spectrum_cut(out, w1[:n], f1[:n], pe1, sig1)


    # #################### OUTPUT ####################
    return out




# =============================================================================
# ============================== CORRELATIONS ==============================
# =============================================================================

# ==============================
def plot_correlations(
                      obs_run : str,
                      file    : str,
                      Mv      : float | None = None,
                      Mk      : float | None = None,
                      root    : str | Path = '.'
    ) -> list[Path]:
    """
    chi2 maps (1, 2 and 3 sigma) in the planes of every pair of parameters
    (translation of gbat_plot_correlations.pro; preliminary in IDL).

    ----------
    Parameters
    ----------
        obs_run : str
            Observing run.

        file : str
            Results of step 2: xdr_dist/<obs_run>/info_dist_<file>.idl (or
            .xdr), e.g. 'HD36512_20111109_032351_M_V85000_vr.dat_v13_m33_Mv356'.

        Mv, Mk : float | None, optional
            If given, a second figure with R, log L and M is created.

        root : str | Path, optional
            Working folder. Default '.'.

    -------
    Returns
    -------
        files : list[Path]
            output_figures/<obs_run>/dist_<file>_1.pdf (and _2.pdf).
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    from pathlib import Path
    import matplotlib.pyplot as plt
    import numpy as np


    # ########## PANORAMA ##########
    from gbat.io import _found, _restore_gbat
    from gbat.plots import _axes, _crea_subplot, _plot_map
    from gbat.variables import GBAT_DIST_DIR, GBAT_FIGURES_DIR, LSP_CONST


    # #################### DATA ####################
    root = Path(root)
    sol = _restore_gbat(_found(root / GBAT_DIST_DIR / obs_run / ('info_dist_' + file)))['sol']
    chi = sol['chi']
    micro = sol['micro'].astype(float)
    lsp = 4. * np.log10(sol['teff'] * 1e3) - sol['grav'] - LSP_CONST


    # #################### FIGURES ####################
    def corner(variables, out):
        n = len(variables)
        pos = _crea_subplot([n - 1] * (n - 1), 0, 0, xmin=0.10, ymin=0.10, xmax=0.95, ymax=0.95)
        fig = plt.figure(figsize=(8.27, 8.27))
        for col in range(n - 1):
            for row in range(col, n - 1):
                (xlab, x), (ylab, y) = variables[col], variables[row + 1]
                ax = _axes(fig, pos[row * (n - 1) + col])
                _plot_map(ax, chi, x, y, xtit=xlab if row == n - 2 else ' ', ytit=ylab if col == 0 else ' ',
                          xticks=row == n - 2, yticks=col == 0)
        fig.savefig(out, format='pdf')
        plt.close(fig)
        return out

    folder = root / GBAT_FIGURES_DIR / obs_run
    folder.mkdir(parents=True, exist_ok=True)
    files = [corner([('Teff [kK]', sol['teff']), ('logg [dex]', sol['grav']), ('logLsp', lsp),
                     (r'Y$_{He}$ x 100', 100. * sol['yhe']), ('micro', micro), ('-logQ', sol['q']),
                     ('beta', sol['beta'])], folder / ('dist_' + file + '_1.pdf'))]
    if Mv or Mk:
        files.append(corner([('Teff [kK]', sol['teff']), ('logg', sol['grav']), ('R', sol['r']),
                             ('logL', sol['l']), ('M', sol['m'])], folder / ('dist_' + file + '_2.pdf')))


    # #################### OUTPUT ####################
    return files
