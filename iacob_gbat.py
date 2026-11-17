from __future__ import annotations

from spec import *

# =========================================================================
# ============================== IACOB-GBAT ==============================
# =========================================================================

# ==============================
def chi_global(
               spectrum   : str | spec,
               obs_run    : str,
               vsini      : float,
               rv0        : float,
               # optional parameters from spec module
               orig       : str = 'IACOB',
               delimiter  : str = ' ',
               outside_dir: str = '',
               # until here from spec module
               zeta_t     : float = 0.,
               vgauss     : float = 0.,
               vr         : float = 0.,
               rv0_cut  : bool = False,
               snr        : float | None = None,
               cut        : bool = False,
               nocut      : bool = False,
               interactive: bool = True,
               fw_ref     : str | None = None,
               Mv         : float | None = None,
               Mk         : float | None = None,
               kud        : bool = False,
               vinf       : float | None = None,
               step2      : bool = False,
               inst       : str | None = None,
               ps         : bool = False,
               spc        : str = '',
               teff       : tuple[float, float] | None = None,
               grav       : tuple[float, float] | None = None,
               yhe        : tuple[float, float] | None = None,
               q          : tuple[float, float] | None = None,
               beta       : tuple[float, float] | None = None,
               micro      : list[str] | None = None,
               root       : str | Path = '.'
    ) -> dict | None:
    """
    Quantitative spectroscopic analysis of an O-type star with a grid of
    FASTWIND models (IACOB-GBAT, step 1 + step 2).

    IACOB-GBAT (Simón-Díaz et al. 2011, JPhCS 328, 012021) is the IACOB
    grid-based automatised tool for the analysis of the H and He lines of
    OB stars. The observed profiles are compared, line by line, with the
    FASTWIND profiles of a grid of models, convolved with vsini, the
    macroturbulence and the instrumental profile. The chi2 distributions of
    the whole (sub)grid give the central values and uncertainties of Teff,
    log g, Y(He), micro, log Q and beta (and R, log L, M and log Mdot if Mv or
    Mk, and vinf, are given).

    This function is the Python implementation of gbat_chi_global.pro (IDL
    v5.0, S. Simón-Díaz, 2025/03). The analysis runs in two steps:

        1. Pre-processing of the lines (only once, saved in
           spectra_cut/<obs_run>/<obs>) and chi2 of each line for every model
           and micro of the subgrid. Saved in xdr_files/<obs_run>/<fileout>.idl.

        2. chi_step2 (see gbat_step2): global chi2 and final results,
           saved in xdr_dist, xdr_results, z_figures and spectra_bfm.

    ----------
    Parameters
    ----------
    spectrum : str | class
        Filename of the input spectrum, or name of the star together with the snr.
        The filename must be unique.

    rv0 : float
        Radial velocity [km/s] used to correct the original spectrum.

    orig, delimiter, outside_dir : see spec module

    obs_run : str
        Observing run: sub-folder of spectra_orig and spectra_cut, and
        label of the line file gbat_linesHHe_<obs_run>.dat.

    vsini : float
        Projected rotational velocity [km/s].

    zeta_t : float, optional
        Radial-tangential macroturbulent velocity [km/s]. Default 0.

    vgauss : float, optional
        Gaussian macroturbulent velocity [km/s] (do not use with zeta_t).
        Default 0.

    vr : float, optional
        Radial velocity [km/s] used to correct the original spectrum.
        Default 0.

    rv0_cut : bool, optional
        Apply the vr correction to an already pre-processed spectrum
        (spectra_cut). Default False.

    snr : float | None, optional
        Global SNR. If None, the SNR is computed line by line from a
        region marked by the user.

    cut : bool, optional
        Pre-process the lines, save them in spectra_cut and exit. Default
        False.

    nocut : bool, optional
        The spectrum is already pre-processed: do not ask for
        re-processing. Default False.

    interactive : bool, optional
        For every line (once per star), show the line and ask whether to
        renormalise a region of it, whether to compute the SNR from a
        region of the continuum, and whether to mark again the limits of
        the line (all marked with the mouse), also with 'nocut' or a
        global 'snr'. The changes are saved in spectra_cut. Needs a
        graphical Matplotlib backend. Default True; False gives the IDL
        behaviour (no questions with 'nocut').

    fw_ref : str | None, optional
        FASTWIND model of the grid used as reference in the
        pre-processing (e.g. 'T370g390He10Q1350b10').

    Mv, Mk : float | None, optional
        Absolute V or Ks magnitude (only one of them), to compute R,
        log L and M in step 2.

    kud : bool, optional
        The synthetic magnitudes of the grid were computed 'ala
        Kudritzki' (R from vsyn). Default False.

    vinf : float | None, optional
        Terminal velocity [km/s], to compute Mdot from Q (needs Mv or Mk).

    step2 : bool, optional
        Go directly to step 2 (step 1 already computed). Default False.

    inst : str | None
        Instrument (see INSTRUMENTS in gbat_variables). Mandatory.

    ps : bool, optional
        Open the final PDF figure with PS_VIEWER. Default False.

    spc : str, optional
        Spectral classification of the star (label). Default ''.

    teff, grav, yhe, q, beta : tuple[float, float] | None, optional
        Ranges of Teff [K], log g [dex], Y(He), -log Q and beta of the
        subgrid. Defaults: [19000, 61000], [2.1, 4.5], [0.04, 0.31],
        [11, 16] and [0.7, 2.1].

    micro : list[str] | None, optional
        Microturbulences of the subgrid, as in the grid ('001', '005' ...).
            Default ['001', '005', '009', '015', '020', '025', '030'].

    root : str | Path, optional
        Working folder, with the folders of IACOB-GBAT (spectra_orig ...)
        and the line files. Default '.'.

    -------
    Returns
    -------
    results : dict | None
        Results of step 2 (see chi_step2), or None with 'cut'.

    -----
    Notes
    -----
    - Output files keep the IDL folders, names and text formats. The IDL
      save files (.xdr) are replaced by pickle files (.idl) with a .json
      copy; IDL .xdr files are read if no .idl file exists, so that the
      results of the IDL version can be used (e.g. step2=True).
    - The subgrid is read from the list of models of the grid
      (see CASA, METAL in gbat_variables).

    ----------
    References
    ----------
    Simón-Díaz, S. et al. 2011, JPhCS 328, 012021.
    Puls, J. et al. 2005, A&A 435, 669.
    Sabín-Sanjulián, C. et al. 2014, A&A 564, A39.
    Holgado, G. et al. 2018, A&A 613, A65.

    --------
    Examples
    --------
    >>> chi_global('HD36512_20111109_032351_M_V85000_vr.dat', 'benchmark', 13., zeta_t=33.,
    ...            inst='hermes85', snr=150, nocut=True, teff=[31900, 36100], grav=[3.89, 4.32],
    ...            q=[12.9, 14.1], yhe=[0.07, 0.19], beta=[0.7, 1.4], Mv=-3.56, kud=True)
    """

    # #################### MODULES & PACKAGES ####################
    from gbat.chi import _chi_one
    from gbat.grid import _load_model, _model_file, _read_models
    from gbat.input import _chi_input, _question, _read_lines_hhe
    from gbat.remote import _prefetch
    from gbat.io import (_dat_line_init, _fa, _ff, _fi, _fstr, _idl_round, _readcol,
                         _save_gbat, _sol_init, _write_spectrum_cut)
    from gbat.step2 import chi_step2
    from gbat.variables import (BETA_RANGE, C_KMS, GBAT_BFM_DIR, GBAT_CUT_DIR,
                                GBAT_DIST_DIR, GBAT_FIGURES_DIR, GBAT_FILES_DIR,
                                GBAT_RESULTS_DIR, GBAT_TRASH_DIR, GRAV_RANGE,
                                MICRO_DEFAULT, Q_RANGE, TEFF_RANGE, YHE_RANGE)


    # #################### INPUTS VALIDATION ####################
    if not isinstance(spectrum, (str, spec)):
        raise TypeError(f"Input spectrum 'spectrum' must be a string or a spec class. Got {type(spectrum).__name__}.")

    for name, value in (('vsini', vsini), ('zeta_t', zeta_t), ('vgauss', vgauss), ('vr', vr)):
        if not isinstance(value, (int, float)) or not np.isfinite(value):
            raise TypeError(f"'{name}' must be a numeric value (int or float). Got {type(value).__name__}.")
    for name, value in (('snr', snr), ('Mv', Mv), ('Mk', Mk), ('vinf', vinf)):
        if value is not None and (not isinstance(value, (int, float)) or not np.isfinite(value)):
            raise TypeError(f"'{name}' must be a numeric value or None. Got {type(value).__name__}.")
    for name, value in (('teff', teff), ('grav', grav), ('yhe', yhe), ('q', q), ('beta', beta)):
        if value is not None and (not isinstance(value, (list, tuple, np.ndarray)) or len(value) != 2):
            raise TypeError(f"Range '{name}' must be a list, tuple or array of two elements. Got {value!r}.")

    # use root as the gbat/ subfolder of the working folder, to avoid problems with relative paths
    root = maindir / 'gbat'

    cfg = _chi_input(obs_run, inst, root) # ** IMPORTANT ** Do not forget to review gbat_variables!


    # ################### WELCOME ####################
    msg.info('Welcome to the GBAT package\nThe grid based automatized tool for the spectroscopic analysis of O-type stars\n')


    # #################### FOLDERS ####################
    for folder in (GBAT_CUT_DIR, GBAT_BFM_DIR, GBAT_FIGURES_DIR, GBAT_DIST_DIR, GBAT_FILES_DIR, GBAT_RESULTS_DIR):
        (root / folder / obs_run).mkdir(parents=True, exist_ok=True)
    (root / GBAT_TRASH_DIR).mkdir(parents=True, exist_ok=True)


    # #################### PARAMETERS ####################
    if zeta_t and vgauss:
        _question('Both zeta_t and vgauss specified, please use only one')
        return None
    teff = teff if teff is not None and np.any(teff) else TEFF_RANGE
    grav = grav if grav is not None and np.any(grav) else GRAV_RANGE
    q = q if q is not None and np.any(q) else Q_RANGE
    yhe = yhe if yhe is not None and np.any(yhe) else YHE_RANGE
    beta = beta if beta is not None and np.any(beta) else BETA_RANGE
    micro = [str(m) for m in micro] if micro is not None and len(micro) else list(MICRO_DEFAULT)
    sigma = 1. / snr if snr else 0.

    if zeta_t:
        fileout = spectrum + '_v' + str(int(vsini)) + '_m' + str(int(zeta_t))
    elif vgauss:
        fileout = spectrum + '_v' + str(int(vsini)) + '_mg' + str(int(vgauss))
    else:
        fileout = spectrum + '_v' + str(int(vsini))


    # #################### STEP 2 DIRECTLY ####################
    if step2:
        print('')
        print(' I GO DIRECTLY TO STEP 2 ... ')
        return chi_step2(fileout, obs_run, Mv=Mv, Mk=Mk, ps=ps, kud=kud, vinf=vinf, spc=spc, inst=inst, root=root)

    print('')
    print(' THE SPECTRUM ' + spectrum + ' WILL BE ANALYZED')


    # #################### LINES ####################
    lines = cfg.lines
    n_lines = len(lines)
    dat_line = _dat_line_init(n_lines, sigma)
    dat_line['line'] = lines
    dat = _read_lines_hhe(root)


    # #################### OBSERVED SPECTRUM ####################
    """
    - Read from spectra_cut (if pre-processed in the current run) or load it from spec module.
    - If vr is given but rv0_cut=False, no vrad correction is applied to a
      pre-processed spectrum; with rv0_cut=True it is.
    """
    cut_dir = root / GBAT_CUT_DIR / obs_run  # Directory with the processed observed spectra
    # check that the input spectrum replacing whatever original extension with '.dat' exists in the cut_dir
    cut_file = cut_dir / Path(spectrum).with_suffix('.dat')
    if cut_file.exists():
        print('')
        print('   - Spectrum read from ' + str(cut_dir) + '/')
        try:
            data = np.loadtxt(cut_file, usecols=range(4), ndmin=2, comments=None)
            wo0, fo0, pe0, sig0 = [data[:, i].copy() for i in range(4)]
        except ValueError:
            wo0, fo0, pe0, sig0 = _readcol(cut_file, fmt='(' + ','.join(['F'] * 4) + ')')

        if rv0_cut:
            wo0 = wo0 * (1. - vr / C_KMS)
    else:
        print('')
        s = spec(spectrum, rv0=rv0, orig=orig, delimiter=delimiter, outside_dir=outside_dir)
        wo0, fo0 = s.wave, s.flux
        pe0, sig0 = np.zeros(wo0.size), np.zeros(wo0.size)
        wo0 = wo0 * (1. - vr / C_KMS)

    # #################### SUBGRID ####################
    print('')
    print(' GRID OF FASTWIND MODELS USED FOR THE ANALYSIS')
    print('')
    print('   - ' + cfg.file_models)
    grid = _read_models(cfg.file_models, teff, grav, yhe, q, beta)
    n_models, n_micro = grid['model'].size, len(micro)
    n_tot = n_models * n_micro
    if n_tot == 0:
        raise ValueError('No models of the grid within the given ranges of parameters.')

    print('')
    print(' SOME NUMBERS OF INTEREST GIVEN THE CONSIDERED INPUT PARAMETERS ')
    print('')
    print('   # models : ' + str(n_models) + '   # micro  : ' + str(n_micro) + '   # lines  : ' + str(n_lines))
    print('')

    sol = _sol_init(n_lines, n_tot)
    tt, ss = 0, 1.
    if fw_ref:
        _question('Model ' + fw_ref + ' used as reference for pre-processing of observed spectrum')


    # #################### DOWNLOAD OF THE SUBGRID (GRID IN A SHARE) ####################
    _prefetch([cfg.casa_red + r + m for r, m in zip(grid['red'], grid['model'])])


    # #################### CHI2 OF EVERY MODEL ####################
    t00 = time.time()
    for i in range(n_models):
        file = _model_file(cfg.casa_red + grid['red'][i] + grid['model'][i])
        if file is None:
            print(cfg.casa_red + grid['red'][i] + grid['model'][i] + '.xdr NOT FOUND!')
            continue
        model = _load_model(file)

        for j in range(n_micro):
            # ########## MODEL ##########
            s = sol[tt]
            s['chi'] = 0.
            s['red'] = grid['red'][i]
            s['model'] = grid['model'][i]
            s['teff'] = model['teff'] * 1e-3
            s['grav'] = model['logg']
            s['ehe'] = model['yhe']
            s['yhe'] = model['yhe']
            s['q'] = 0.1 * _idl_round(-1e1 * model['q'])
            s['beta'] = model['beta']
            s['micro'] = micro[j]
            s['vsyn'] = model['vsyn']
            s['v'] = model['v']
            s['ks'] = model['ks']

            # ########## LINES ##########
            for mm in range(n_lines):
                chi0, err_mes = _chi_one(dat_line[mm], tt, wo0, fo0, pe0, sig0, model, micro[j], vsini, cfg, dat,
                                         zeta_t=zeta_t, vgauss=vgauss, fw_ref=fw_ref, snr=snr, nocut=nocut,
                                         interactive=interactive)
                if err_mes == 1:
                    print('WARNING ... ' + s['model'] + ' ' + dat_line['line'][mm] + ' ' + s['micro'])
                s['chi_lin'][mm] = chi0     # Computed chi2_tilda for line mm (model tt)
                s['chi'] = s['chi'] + chi0  # Computed chi2_tilda for model tt (1st pass)

            # ########## SAVE THE PRE-PROCESSED SPECTRUM ##########
            if tt == 0:
                print('')
                print(' OBSERVED SPECTRUM PROCESSED AND SAVED IN ')
                print('   - ' + str(cut_dir / spectrum))
                print('')
                _write_spectrum_cut(cut_dir / spectrum, wo0, fo0, pe0, sig0)
                if cut:
                    return None
                print(' IACOB-GBAT running ... please wait')
                print(' ')
                t00 = time.time()
                print('  N_TOT   #MODELS    #TIME(min)  #ESTIMAT. TIME (min)')
                print(' ------   -------    ----------  --------------------')
            elif tt == _idl_round(ss * n_tot / 20.):
                dt = time.time() - t00
                print(' ' + _fi(n_tot, 6) + '    ' + _fi(tt, 6) + '       ' + _ff(dt / 60., 7, 3) + '      '
                      + _ff(dt * 20. / ss / 60., 7, 3) + ' ' * 15 + _fa('(', 1) + _fa(_fstr(ss * 5.), 2) + _fa('%)', 2))
                ss += 1.

            tt += 1
    print(' ')


    # #################### MODELS NOT FOUND ####################
    """
    - In IDL, the entries of sol of the models not found are kept empty
      (chi = 0), which breaks step 2. They are removed here.
    """
    if tt < n_tot:
        print(f' WARNING ... {n_tot - tt} model/micro combinations not found in the grid (removed from sol)')
        sol = sol[:tt]
    if tt == 0:
        raise FileNotFoundError(f"None of the models of the subgrid was found in '{cfg.casa_red}'.")


    # #################### SAVE STEP 1 ####################
    base = root / GBAT_FILES_DIR / obs_run / fileout
    print(' SAVING RESULTS FROM STEP 1 IN')
    print('')
    print('   - ' + str(base) + '.idl')
    _save_gbat(base, sol=sol, casa_red=cfg.casa_red, vsini=float(vsini), zeta_t=float(zeta_t), vgauss=float(vgauss),
               resol=cfg.resol, res0=cfg.res0, wo0=wo0, fo0=fo0, pe0=pe0, dat_line=dat_line)


    # #################### STEP 2 ####################
    print('')
    print(' STEP 1 COMPLETED ... ')


    # #################### OUTPUT ####################
    return chi_step2(fileout, obs_run, Mv=Mv, Mk=Mk, ps=ps, kud=kud, vinf=vinf, spc=spc, inst=inst, root=root)




# ==============================
def launcher(
             file    : str,
             nocut   : bool = False,
             interactive : bool = True,
             cut     : bool = False,
             step2   : bool = False,
             star    : str | None = None,
             vr      : float = 0.,
             rv0_cut : bool = False,
             it0     : int | None = None,
             it1     : int | None = None,
             ps      : bool = False,
             root    : str | Path = '.',
             orig    : str = 'IACOB', # for spec module
             delimiter : str = ' ', # for spec module
             outside_dir : str = '' # for spec module
    ) -> list[dict | None]:
    """
    Run IACOB-GBAT for a list of stars (translation of gbat_launcher.pro).

    ----------
    Parameters
    ----------
        file : str
            Table in z_tables_launcher, with one star per line and the columns
                Name run inst vsini vmac Teffm TeffM loggm loggM Qm QM bm bM
                Hem HeM micm micM fw_ref SPC Mv vinf
            (e.g. 'table_benchmark.dat'). Lines that do not follow this
            format (e.g. the header) are skipped.

        nocut, interactive, cut, step2, vr, rv0_cut, ps
            See chi_global.

        star : str | None, optional
            Analyse only the first star whose file name contains '<star>_'.

        it0, it1 : int | None, optional
            First and last line (0-based) of the table to analyse.

        root : str | Path, optional
            Working folder. Default '.'.

    -------
    Returns
    -------
        results : list[dict | None]
            Results of chi_global for each star.

    -----
    Notes
    -----
        - The microturbulences of the subgrid are those within [micm, micM]
          of ['001', '005', '009', '015', '020', '025', '030'] if
          mean([micm, micM]) >= 10, or of ['001', '003', '005', '007', '009',
          '011', '015', '020', '025', '030'] otherwise.
        - The analysis uses snr=150 and kud=1, as in IDL.

    --------
    Examples
    --------
        >>> launcher('table_benchmark.dat', nocut=True, ps=True)
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    from pathlib import Path
    import numpy as np


    # ########## PANORAMA ##########
    from gbat.io import _readcol
    from gbat.variables import GBAT_LAUNCHER_DIR, LAUNCHER_SNR, MICRO_HIGH, MICRO_LOW


    # #################### TABLE ####################
    root = maindir / 'gbat'
    table = root / GBAT_LAUNCHER_DIR / file
    if not table.exists():
        raise FileNotFoundError(f"Launcher table '{table}' not found.")
    (sp_file, run, inst, vs, vm, teffmin, teffmax, gravmin, gravmax, qmin, qmax, bmin, bmax, hemin, hemax,
     micromin, micromax, fwref0, spc0, mv0, vinf0, rv0) = _readcol(table, fmt='(A,A,A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,A,A,F,F,F)')
    if sp_file.size == 0:
        raise ValueError(f"No valid lines found in '{table}'.")


    # #################### STARS ####################
    it0 = it0 or 0
    it1 = sp_file.size - 1 if it1 is None else it1
    if it1 >= sp_file.size:
        raise IndexError(f"it1 ({it1}) is out of bounds for the table '{table}' with {sp_file.size} entries.")

    if star:
        tt = [k for k, s in enumerate(sp_file) if star + '_' in s]
        if not tt:
            raise ValueError(f"Star '{star}' not found in '{table}'.")
        it0 = it1 = tt[0]


    # #################### ANALYSIS ####################
    results = []
    for i in range(it0, it1 + 1):
        micro0 = MICRO_HIGH if np.mean([micromin[i], micromax[i]]) >= 10 else MICRO_LOW
        micro = [m for m in micro0 if micromin[i] <= float(m) <= micromax[i]]
        results.append(chi_global(
                                  spectrum=str(sp_file[i]),
                                  obs_run=str(run[i]),
                                  vsini=float(vs[i]),
                                  rv0=float(rv0[i]), # for spec module
                                  orig=orig, # for spec module
                                  delimiter=delimiter, # for spec module
                                  outside_dir=outside_dir, # for spec module
                                  zeta_t=float(vm[i]),
                                  nocut=nocut,
                                  interactive=interactive,
                                  cut=cut, step2=step2, inst=str(inst[i]), teff=[teffmin[i], teffmax[i]],
                                  grav=[gravmin[i], gravmax[i]], q=[qmin[i], qmax[i]], yhe=[hemin[i], hemax[i]],
                                  beta=[bmin[i], bmax[i]], micro=micro, fw_ref=str(fwref0[i]), snr=LAUNCHER_SNR,
                                  spc=str(spc0[i]), kud=True, Mv=float(mv0[i]) or None, vinf=float(vinf0[i]) or None,
                                  vr=vr, rv0_cut=rv0_cut, ps=ps, root=root))


    # #################### OUTPUT ####################
    return results
