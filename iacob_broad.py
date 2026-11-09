from __future__ import annotations

from spec import *

def iacob_broad(
                # spectrum is a string or a class from pyIACOB spec module
                spectrum   : str | spec,
                snr        : float | None = None,
                rv0        : float = 0.0,
                offset     : float = 0.0,
                cut_edges  : bool = False,
                orig       : str = 'IACOB',
                delimiter  : str = ' ',
                outside_dir: str = '',
                # until here, the items from pyIACOB spec module
                line       : str = '',
                output_dir : str = maindir / 'broadening',
                cut        : bool = False,
                noproc     : bool = False,
                resul      : bool = False,
                xl         : np.ndarray[float] | None = None,
                xyc        : np.ndarray[float] | None = None,
                xrec       : np.ndarray[float] | None = None,
                xsnr       : np.ndarray[float] | None = None,
                vs_user    : float = 0.,
                vm_user    : float = 0.,
                vmin       : float | None = None,
                silen      : bool = False,
                rerec      : bool = False,
                em         : bool = False,
                beta       : float = 1.5,
                sigma_FT0  : float | None = None,
                w_intprof  : np.ndarray | None = None,
                f_intprof  : np.ndarray | None = None,
                qib_step   : bool = True,
                ) -> dict | None:

    """
    Projected rotational velocity (vsini) and macroturbulent broadening (vmac) of a
    spectral line with the Fourier Transform + Goodness-of-Fit method.

    IACOB-BROAD (Simón-Díaz & Herrero 2014, A&A 562, A135) is a
    semi-automatic tool for the line-broadening characterisation of OB-type
    stars. It disentangles the rotational and macroturbulent contributions to
    the profile of a spectral line by combining two techniques:

        - Fourier Transform (FT): vsini is obtained from the first zero of the FT
        of the line profile (sigma_FT0 = 0.660 for a linear limb darkening with
        epsilon = 0.6).

        - Goodness-of-Fit (GOF): the observed profile is compared, through chi2,
        with a grid of synthetic profiles built by convolving an intrinsic
        profile with rotation, macroturbulence (radial-tangential or Gaussian)
        and the instrumental profile. The best values and their 1-sigma
        uncertainties come from the lower envelope of the chi2 distributions.

    This function is the Python implementation of the original IDL tool
    (v9.0, S. Simón-Díaz, 2024/11/20) developed by D. Galán together with the
    other sub-modules included in the broadening/ folder.

    The analysis of one line runs in three stages, each one saved to disk so that
    it can be skipped in later runs:

        1. Pre-processing: spectral window, continuum level, SNR, interpolation to
        a fixed wavelength step and vsini(FT), marked on the FT by the user.
        Saved in '<output_dir>/output_idl/<star>_<line>.idl'; the analysis stops here with 'cut'.

        2. GOF grid in vsini x vmac x EW. Saved in '<output_dir>/output_idl/<star>_<line>_resul.idl'.
        Each .idl file is accompanied by a human-readable copy '<output_dir>/output_json/*.json'.

        3. Final estimates, summary figure '<output_dir>/output_pdf/<star>_<line>.pdf'
        and table of results in the terminal.

    Parameters
    ----------
    spectrum : str | class
        Filename of the input spectrum, or name of the star together with the snr.
        The filename must be unique.

    snr, rv0, offset, cut_edges, orig, delimiter, outside_dir : see spec module

    line : str, optional
        Name of the line in pyiacob_broad_lines.dat (e.g. 'SiIII4552'). It sets
        which line is analysed, the initial +/- 20 A window and the output file names.
        Default '' (the region is marked interactively).

    output_dir : str | Path, optional
        Folder where the output files are saved, created if it does not
        exist. Inside it, three subfolders are always created:
        'output_idl' (.idl save files), 'output_json' (their readable .json
        copies) and 'output_pdf' (summary figure).
        Default is the broadening subfolder inside the main directory.

    cut : bool, optional
        Stop after the pre-processing and vsini computation (Fourier Transform) of the line.
        Default False.

    noproc : bool, optional
        Skip the pre-processing and use the saved file '<output_dir>/output_idl/<star>_<line>.idl'
        (or, if not found, '<star>_<line>.xdr' in the same folder, e.g. from the
        IDL version), only if it exists. Default False.

    resul : bool, optional
        Skip the pre-processing and the GOF and use the saved files in
        '<output_dir>/output_idl' (.idl, or else .xdr), only if they exist;
        implies 'noproc'. Default False.

    xl : np.ndarray[float] | None = None, optional
        Limits of the initial spectral window [A] around the line. If None, marked with two clicks.
        Format: [x1, x2]. Requires 'line'. Default None.

    xyc : np.ndarray[float] | None = None, optional
        Wavelength of the line core [A] and flux of the adjacent continuum.
        If None, marked with one click. Format: [xc, yc]. Requires 'line'. Default None.

    xrec : np.ndarray[float] | None = None, optional
        Two continuum windows [A] for the local renormalisation (only used
        with 'rerec'). If None, marked with four clicks. Format: [x1, x2, x3, x4].
        Requires 'line'. Default None.

    xsnr : np.ndarray[float] | None = None, optional
        Continuum window [A] used to compute the SNR. If None, marked with two clicks. Format:
        [x1, x2]. Requires 'line'. Default None.

    vs_user : float, optional
        Projected rotational velocity [km/s] provided by the user, drawn in the final
        figure for comparison (vs_user <= 5 km/s means no rotation). Default 0 (not drawn).

    vm_user : float, optional
        Macro-turbulence [km/s] provided by the user, drawn in the final figure
        for comparison (vs_user <= 5 km/s means no macroturbulence). Default 0 (not
        drawn).

    vmin : float | None, optional
        Minimum projected rotational velocity of the GOF grid [km/s], > 0. Default: one pixel
        of the line in velocity units.

    silen : bool, optional
        Do not print the table of results. Default False.

    rerec : bool, optional
        Local renormalisation of the continuum with a straight line fitted to
        the 'xrec' windows. Default False.

    rv0 : float, optional
        Radial velocity correction applied to the wavelengths [km/s].
        Default 0.

    em : bool, optional
        Emission line (the delta-like intrinsic profile is flipped).
        Default False.

    beta : float, optional
        Limb-darkening coefficient of the rotational profile,
        beta = epsilon/(1 - epsilon), in [0, 2].
        Default is 1.5 (epsilon = 0.6, linear limb darkening).

    sigma_FT0 : float, optional
        Position of the first zero of the rotational FT, in units of
        lambda*vsini/c, > 0. Default None (0.660).

    w_intprof : np.ndarray[float] | None, optional
        Wavelength of the intrinsic profile provided by the user: wavelength [A], strictly
        increasing, and normalised flux (same size). Both or none must be
        given. Default None (delta-like Gaussian profile).

    f_intprof : np.ndarray[float] | None, optional
        Flux of the intrinsic profile provided by the user: wavelength [A], strictly
        increasing, and normalised flux (same size). Both or none must be
        given. Default None (delta-like Gaussian profile).

    qib_step : bool, optional
        Ask the user for a quality index (QIB) to the analysis (0 to 5).
        Default is True.

    Returns
    -------
    results : dict or None
    A dictionary containing the following key-value pairs:
        - 'filename' : str
            Identification of the spectrum.

        - 'line' : str
            Identification of the spectral line.

        - 'ew' : float
            Equivalent width from the GOF [mA].

        - 'snr' : float
            Signal-to-noise ratio of the adjacent continuum.

        - 'vs0' : float
            vsini from the GOF assuming vmac = 0 [km/s].

        - 'vft' : float
            vsini from the first zero of the FT [km/s].

        - 'vmft' : float
            vmac from the GOF for vsini = vsini(FT) [km/s

        - 'evmft' : list[float, float]
            Lower and upper 1-sigma errors of 'vmft' [km/s].

        - 'vsgof' : float
            vsini from the GOF [km/s].

        - 'evsgof' : list[float, float]
            Lower and upper 1-sigma errors of 'vsgof' [km/s].

        - 'vmgof' : float
            vmac from the GOF [km/s].

        - 'evmgof' : list[float, float]
            Lower and upper 1-sigma errors of 'vmgof' [km/s].

    References
    ----------
        Simón-Díaz, S. & Herrero, A. 2014, A&A 562, A135.

    Examples
    --------
        >>> from iacob_broad import *
        >>> results = iacob_broad('HD37042', snr='best', line='SiIII')              # interactive
        >>> results = iacob_broad('HD37042', snr='best', line='SiIII', cut=True)    # pre-processing only
        >>> results = iacob_broad('HD37042', snr='best', line='SiIII', noproc=True) # skip pre-processing if already done and GOF is computed
    """

    # #################### IMPORT SUB-PACKAGES ####################
    from broadening.gof import _best_models, _env_dist, _final_estimates, _gof_grid
    from broadening.interactive import _read_linelist
    from broadening.results import _plot_summary, _print_results, _restore_idl, _restore_json, _save
    from broadening.variables import FACT_VRT, FACT_VRT_HEI, FACT_VS
    from broadening.vsiniFT import _vsiniFT

    # #################### INPUTS VALIDATION ####################
    num_types = (int, float, np.integer, np.floating)

    # Validate input spectrum 'spectrum' that can be a string or a spec class
    if not isinstance(spectrum, (str, spec)):
        raise TypeError(f"Input spectrum 'spectrum' must be a string or a spec class. Got {type(spectrum).__name__}.")

    # Validate spectral line 'line'
    if not isinstance(line, str):
        raise TypeError(f"Spectral line 'line' must be a string. Got {type(line).__name__}.")

    # Validate output folder 'output_dir'
    if not isinstance(output_dir, (str, Path)):
        raise TypeError(f"Output folder 'output_dir' must be a string or Path. Got {type(output_dir).__name__}.")

    output_dir = Path(output_dir)
    if output_dir.exists() and not output_dir.is_dir():
        raise NotADirectoryError(f"Output folder 'output_dir' is not a directory: '{output_dir}'.")
    output_dir.mkdir(parents=True, exist_ok=True)
    for sub_dir in ('output_json', 'output_pdf', 'output_idl'):
        (output_dir / sub_dir).mkdir(exist_ok=True)

    # Validate boolean options
    if not isinstance(cut, bool):
        raise TypeError(f"Option for stop after the pre-processing process 'cut' must be a boolean. Got {type(cut).__name__}.")
    if not isinstance(noproc, bool):
        raise TypeError(f"Option for skip the pre-processing process 'noproc' must be a boolean. Got {type(noproc).__name__}.")
    if not isinstance(resul, bool):
        raise TypeError(f"Option for skip the pre-processing and the GOF process 'resul' must be a boolean. Got {type(resul).__name__}.")

    # Validate limits of the line 'xl'
    if xl is not None:
        if not (isinstance(xl, (list, tuple, np.ndarray)) and len(xl) == 2):
            raise TypeError(f"Limits of the line 'xl' must be a list, tuple or array of two elements. Got {type(xl).__name__} with {len(xl)} elements.")
        xl = np.asarray(xl, dtype=float)

    # Validate wavelength of the line core and flux of the adjacent continuum 'xyc'
    if xyc is not None:
        if not (isinstance(xyc, (list, tuple, np.ndarray)) and len(xyc) == 2):
            raise TypeError(f"Wavelength of the line core and flux of the adjacent continuum 'xyc' must be a list, tuple or array of two elements. Got {type(xyc).__name__} with {len(xyc)} elements.")
        xyc = np.asarray(xyc, dtype=float)

    # Validate continuum windows used for the local renormalisation 'xrec'
    if xrec is not None:
        if not (isinstance(xrec, (list, tuple, np.ndarray)) and len(xrec) == 4):
            raise TypeError(f"Two continuum windows for the local renormalisation 'xrec' must be a list, tuple or array of four elements. Got {type(xrec).__name__} with {len(xrec)} elements.")
        xrec = np.asarray(xrec, dtype=float)

    # Validate continuum window used to compute the SNR 'xsnr'
    if xsnr is not None:
        if not (isinstance(xsnr, (list, tuple, np.ndarray)) and len(xsnr) == 2):
            raise TypeError(f"Continuum window used to compute the SNR 'xsnr' must be a list, tuple or array of two elements. Got {type(xsnr).__name__} with {len(xsnr)} elements.")
        xsnr = np.asarray(xsnr, dtype=float)

    # Validate user-provided vsini 'vs_user'
    if isinstance(vs_user, (bool, np.bool_)) or not isinstance(vs_user, num_types):
        raise TypeError(f"User-provided vsini 'vs_user' must be a numeric value (int or float). Got {type(vs_user).__name__}.")
    vs_user = float(vs_user)

    # Validate user-provided vmac 'vm_user'
    if isinstance(vm_user, (bool, np.bool_)) or not isinstance(vm_user, num_types):
        raise TypeError(f"User-provided vmac 'vm_user' must be a numeric value (int or float). Got {type(vm_user).__name__}.")
    vm_user = float(vm_user)

    # Validate minimum vsini and vmac of the GOF grid vmin
    if vmin is not None:
        if isinstance(vmin, (bool, np.bool_)) or not isinstance(vmin, num_types):
            raise TypeError(f"Minimum vsini and vmac of the GOF grid 'vmin' must be a numeric value (int or float). Got {type(vmin).__name__}.")
        vmin = float(vmin)

    # Validate verbosity & continuum options
    if not isinstance(silen, bool):
        raise TypeError(f"Option for silence the table of results 'silen' must be a boolean. Got {type(silen).__name__}.")
    if not isinstance(rerec, bool):
        raise TypeError(f"Option for local renormalisation of the continuum 'rerec' must be a boolean. Got {type(rerec).__name__}.")

    # Validate radial velocity correction 'rv0'
    if isinstance(rv0, (bool, np.bool_)) or not isinstance(rv0, num_types):
        raise TypeError(f"Radial velocity correction 'rv0' must be a numeric value (int or float). Got {type(rv0).__name__}.")
    rv0 = float(rv0)

    # Validate emission line option
    if not isinstance(em, bool):
        raise TypeError(f"Option for emission line 'em' must be a boolean. Got {type(em).__name__}.")

    # Validate limb-darkening coefficient 'beta'
    if isinstance(beta, (bool, np.bool_)) or not isinstance(beta, num_types):
        raise TypeError(f"Limb-darkening coefficient 'beta' must be a numeric value (int or float). Got {type(beta).__name__}.")
    if not (0.0 <= beta <= 2.0):
        raise ValueError(f"Limb-darkening coefficient 'beta' must be in the range [0, 2]. Got {beta}.")
    beta = float(beta)

    # Validate sigma_FT0
    if sigma_FT0 is not None:
        if isinstance(sigma_FT0, (bool, np.bool_)) or not isinstance(sigma_FT0, num_types):
            raise TypeError(f"Position of the first zero of the rotational FT 'sigma_FT0' must be a numeric value (int or float). Got {type(sigma_FT0).__name__}.")
        sigma_FT0 = float(sigma_FT0)

    # Validate intrinsic profile pair
    if (w_intprof is None) != (f_intprof is None):
        raise ValueError(f"Both 'w_intprof' and 'f_intprof' must be provided, or none of them. Got w_intprof={w_intprof} and f_intprof={f_intprof}.")
    if w_intprof is not None and f_intprof is not None:
        if not (isinstance(w_intprof, (list, tuple, np.ndarray)) and isinstance(f_intprof, (list, tuple, np.ndarray))):
            raise TypeError(f"Intrinsic profile 'w_intprof' and 'f_intprof' must be lists, tuples or arrays. Got {type(w_intprof).__name__} and {type(f_intprof).__name__}.")
        w_intprof = np.asarray(w_intprof, dtype=float)
        f_intprof = np.asarray(f_intprof, dtype=float)
        if w_intprof.shape != f_intprof.shape:
            raise ValueError(f"Intrinsic profile 'w_intprof' and 'f_intprof' must have the same shape. Got {w_intprof.shape} and {f_intprof.shape}.")
        if w_intprof.ndim != 1:
            raise ValueError(f"Intrinsic profile 'w_intprof' and 'f_intprof' must be 1D arrays. Got {w_intprof.ndim}D.")
        if not np.all(np.diff(w_intprof) > 0):
            raise ValueError(f"Intrinsic profile wavelength 'w_intprof' must be strictly increasing. Got {w_intprof}.")

    # Validate option for QIB prompt step 'qib_step'
    if not isinstance(qib_step, bool):
        raise TypeError(f"Option for QIB prompt step 'qib_step' must be a boolean. Got {type(qib_step).__name__}.")


    # #################### SETTINGS ####################
    # ########## SPECTRAL LINES ##########
    """
    - The line list is read from the file 'iacob_broad/lines.dat'.
      It contains the line names and their central wavelengths [A].
    """
    lines = _read_linelist() # dict with line names and central wavelengths [A]
    if line and line not in lines:
        raise ValueError(f"Spectral line 'line' = '{line}' is not in lines.dat. Available lines: {', '.join(lines)}.")
    line_lab = f'_{line}' if line else '' # Label for the line in the output files

    # ########## FIT PROCESS ##########
    xr0 = lines[line] + np.array([-20., 20.]) if line else None # Initial window of +/- 20 A around the line (None -> marked interactively)
    fact_vrt = FACT_VRT_HEI if 'HEI' in line.upper() else FACT_VRT # Maximum vmac of the GOF grid relative to vsini(FT) (larger for HeI lines)

    # #################### CONSTRUCT THE SPEC CLASS FOR THE STAR ####################
    if isinstance(spectrum, str):
        if snr is None:
            raise ValueError("Input spectrum 'spectrum' is a string, but 'snr' is None. Please provide the SNR of the spectrum.")
        if not isinstance(snr, (int, float)):
            raise TypeError(f"Input SNR 'snr' must be a numeric value (int or float). Got {type(snr).__name__}.")
        snr = float(snr)
        s = spec(spectrum, snr=snr, rv0=rv0, offset=offset, cut_edges=cut_edges,
                 orig=orig, delimiter=delimiter, outside_dir=outside_dir)
    else:
        s = spectrum

    filename = Path(s.filename).stem
    res = s.resolution
    while res is None or (isinstance(res, float) and np.isnan(res)):
        try:
            val = float(input("Resolving power of the spectrum (lambda/dlambda) > 0: "))
            if val <= 0:
                msg.error("Error: Resolving power must be > 0.")
                continue
            res = val
        except ValueError:
            msg.error("Error: Please enter a valid numeric value.")


    # ########## FILES ##########
    idl_dir = output_dir / 'output_idl'
    json_dir = output_dir / 'output_json'
    pdf_dir = output_dir / 'output_pdf'

    # Pre-processing of the line (step 1)
    prep_idl = idl_dir / f'{filename}{line_lab}.idl'
    prep_xdr = idl_dir / f'{filename}{line_lab}.xdr'
    prep_json = json_dir / f'{filename}{line_lab}.json'

    # Goodness-of-fit (step 2)
    gof_file = idl_dir / f'{filename}{line_lab}_resul.idl'
    gof_xdr = idl_dir / f'{filename}{line_lab}_resul.xdr'
    gof_json = json_dir / f'{filename}{line_lab}_resul.json'

    # Summary figure (step 3)
    figure = str(pdf_dir / f'{filename}{line_lab}.pdf')

    # Saved file to read: .idl, or else .xdr (IDL version / old runs) in 'output_idl'
    found = lambda f: f if f.exists() or not f.with_suffix('.xdr').exists() else f.with_suffix('.xdr')

    noproc = noproc or resul


    # #################### PREPROCESSING & VSINI(FT) ####################
    """
    - In this section, the line is pre-processed to prepare it for the computation. That is,
        >> the spectral window is defined,
        >> the continuum level is determined,
        >> the SNR is computed, and
        >> the line is interpolated to a fixed wavelength step.

    - Then, the user marks the vsini(FT) on the Fourier transform of the line.

    - The results are saved in 'output_idl/<star>_<line>.idl' (and 'output_json/<star>_<line>.json').

    - If the 'cut' option is True, the function stops here and returns None.

    - Skipped with 'noproc' (or 'resul') only when '<star>_<line>.idl' (or '.xdr') already exists in 'output_idl'.
    - New results are always written as .idl.
    """
    exists_prep = prep_idl.exists() or prep_xdr.exists() or prep_json.exists()

    if exists_prep and not noproc:
        while True:
            choice = input(f"\nPre-processing files for '{filename}', line '{line}' already exist. Do you want to run the analysis again and overwrite them? [y/n]: ").strip().lower()
            if choice in ('y', 'yes'):
                msg.info("Re-running pre-processing and overwriting existing files...")
                noproc = False
                break
            elif choice in ('n', 'no'):
                msg.info("Skipping pre-processing (using existing files)...")
                noproc = True
                break
            else:
                msg.error("Please enter 'y' or 'n'.")

    if not noproc or not exists_prep:
        if not exists_prep:
            msg.info('\nPre-processing file not found ... lets create it ...')
        prep = _vsiniFT(
                        list(lines.values()),
                        w = s.wave,
                        f = s.flux,
                        res = res,
                        xr0 = xr0,
                        xl = xl,
                        xyc = xyc,
                        xsnr = xsnr,
                        xrec = xrec,
                        rerec = rerec,
                        vmin = vmin,
                        beta = beta,
                        sigma_FT0 = sigma_FT0,
                        w_intprof = w_intprof,
                        f_intprof = f_intprof
                        )
        _save(prep_idl, output_dir=json_dir, **prep)
        if cut:
            # Stop after the pre-processing (no results yet)
            return None
    else:
        # Pre-processing was skipped: if the binary .idl/.xdr is missing, restore from .json
        if not prep_idl.exists() and not prep_xdr.exists() and prep_json.exists():
            msg.info(f"\nPre-processing binary file not found. Restoring from '{prep_json.name}' and regenerating '{prep_idl.name}'...")
            prep_data = _restore_json(prep_json)
            _save(prep_idl, readable=False, **prep_data)


    # #################### GOODNESS OF FIT ####################
    """
    - The observed line is compared, through chi2, with a grid of synthetic profiles:
        intrinsic profile x rotation x radial-tangential macroturbulence x instrument,
        that is, intrinsic profile x vsini x vmac x EW.

    - The results are saved in 'output_idl/<star>_<line>_resul.idl' (and 'output_json/<star>_<line>_resul.json').

    - Skipped with 'resul' only when '<star>_<line>_resul.idl' (or '.xdr') already exists in 'output_idl'.
    """
    exists_gof = gof_file.exists() or gof_xdr.exists() or gof_json.exists()

    if not resul or not exists_gof:
        if not exists_gof:
            msg.info(f'\n{gof_file.name} file not found ... lets create it ...')
        gof = _gof_grid(
                        _restore_idl(found(prep_idl)),
                        res,
                        fact_vrt,
                        vmin = vmin,
                        em = em,
                        beta = beta
                        )
        _save(gof_file, output_dir=json_dir, **gof)
    else:
        # resul=True: if the binary GOF file is missing, restore from .json
        if not gof_file.exists() and not gof_xdr.exists() and gof_json.exists():
            msg.info(f"\nGOF binary file not found. Restoring from '{gof_json.name}' and regenerating '{gof_file.name}'...")
            gof_data = _restore_json(gof_json)
            _save(gof_file, readable=False, **gof_data)

    # #################### QUALITY INDEX (QIB) STEP ####################
    if qib_step:
        qib = -1.0
        while qib is -1.0:
            try:
                val = input("Give a quality index to the analysis (0 to 5): ")
                qib = float(val)
                if not (0.0 <= qib <= 5.0):
                    msg.error("Error: QIB must be a value between 0 and 5.")
                    qib = -1.0
            except ValueError:
                msg.error("Error: Please enter a valid numeric value.")
    else:
        qib = -1.0


    # #################### FINAL ESTIMATES ####################
    # ########## SAVED RESULTS ##########
    prep, gof = _restore_idl(found(prep_idl)), _restore_idl(found(gof_file))
    vft = prep['vft'] # vsini from the first zero of the FT [km/s]
    vsf = gof.get('vsf', FACT_VS * vft) # Maximum vsini of the grid [km/s]
    vrtf = gof.get('vrtf', fact_vrt * vft) # Maximum vmac of the grid [km/s]


    # ########## ESTIMATES ##########
    est = _final_estimates(gof, vft) # vsini and vmac (with errors) from the chi2 distributions
    env_vs = _env_dist(gof['chi2_vs'], gof['chi2']) # Envelope of chi2 vs vsini (for the figure)
    env_vm = _env_dist(gof['chi2_vrt'], gof['chi2']) # Envelope of chi2 vs vmac (for the figure)


    # #################### SUMMARY FIGURE ####################
    """
    - Synthetic profiles of the solutions and their Fourier transforms,
      compared with the observed line and the chi2 distributions.
    """
    models, fts = _best_models(
                               prep,
                               est,
                               vsf,
                               res,
                               em = em,
                               beta = beta,
                               vs_user = vs_user,
                               vm_user = vm_user,
                               sigma_FT0 = sigma_FT0
                               )

    _plot_summary(
                  figure,
                  filename,
                  line,
                  prep,
                  gof,
                  est,
                  models,
                  fts,
                  env_vs,
                  env_vm,
                  vsf,
                  vs_user = vs_user,
                  vm_user = vm_user
                 )


    # #################### RESULTS ####################
    results = dict(
                   filename = filename,
                   line = line,
                   ew = est.ew,
                   snr = prep['snr'],
                   vs0 = est.vs0,
                   vft = vft,
                   vmft = est.vmft,
                   evmft = est.evmft,
                   vsgof = est.vsgof,
                   evsgof = est.evsgof,
                   vmgof = est.vmgof,
                   evmgof = est.evmgof,
                  )
    if not silen:
        _print_results(results, 'w_intprof' in prep)


    # ########## FINAL SAVE ##########
    results = gof.get('results', results)
    _save(
          gof_file,
          output_dir = json_dir,
          sigma_chi = gof['sigma_chi'],
          chi2 = gof['chi2'],
          chi2_vs = gof['chi2_vs'],
          chi2_vrt = gof['chi2_vrt'],
          chi2_ew = gof['chi2_ew'],
          Avs = gof['avs'],
          vrtf = vrtf,
          vsf = vsf,
          results = results,
          qib = qib
         )

    return results




def ib_results(solution_dir=None, table=None, check_best=False, format='fits'):

    """
    Function to generate a table with the results from IACOB-broad.
    If a table is given with a list of stars, it will iterate over the table to
    search only for the corresponding results.

    Parameters
    ----------
    solution_dir : str
        Directory of IACOB-Broad directory where to search for the solution files.
        It must be the root where the output_idl/output_json/output_pdf folders are.
        Default is the broadening subfolder inside the main directory.

    table : str
        Name of the input table containing the list of stars to search, with the
        column 'ID', 'filename', 'line' and optionally 'QIB'.
        Defeult is None, meaning that it will gather all the results in the folder.
        Note: the input table must be in fits or csv format.

    check_best : boolean, optional
        True if each spectra from the .xdr file is checked against the best
        spectrum in the database. Default is False.

    format : str, optional
        Enter the output format for the table: 'fits' (default), 'ascii' or 'csv'.

    Note: the option for using this over ascii input spectra is not yet implemented.

    Returns
    -------
    Nothing but the output table with the IACOB broad results is generated.
    """

    import json

    solution_dir = Path(solution_dir or maindir / 'broadening').expanduser() / 'output_json'

    files = []
    QIB = []
    if table is not None:
        table = findtable(table, delimiter=' ', guess=True)

        for row in table:
            matches = [f for f in solution_dir.glob('*_resul.json') if f.name.split('_')[0] == row['ID']]

            if not matches:
                msg.warn(f"\nNo .json file found for {row['ID']}. Continuing...")
                continue
            elif len(matches) > 1:
                msg.warn(f"\More than one .json (line) file found for {row['ID']}:")
                msg.warn([f.stem.replace('_resul', '') for f in matches])
                msg.info('Using the one provided in the input table...')

            line_matches = [f for f in matches if row['line'] in f.name]
            if not line_matches:
                msg.warn(f"\nNo .json file matching line '{row['line']}' found for {row['ID']}. Continuing...")
                continue

            files.append(line_matches[0])
            QIB.append(row['QIB'] if 'QIB' in row.colnames else -1.0)
    else:
        files = sorted(solution_dir.glob('*_resul.json'))
        QIB = [-1.0] * len(files)

    bar = pb.ProgressBar(maxval=len(files),
          widgets=[pb.Bar('=', '[', ']'), ' ', pb.Percentage()])
    bar.start()

    data_rows = []
    for i, (file, qib) in enumerate(zip(files, QIB)):
        with open(file, 'r') as f:
            raw_data = json.load(f)

        res = raw_data['results']
        res_id = Path(res['filename']).stem.split('_')[0]
        data_rows.append(
            {
                'ID': res_id,
                'filename': res['filename'],
                'vsini_FT': res['vft'],
                'vmac_FT': res['vmft'],
                'vmac_FT_eDW': res['evmft'][0],
                'vmac_FT_eUP': res['evmft'][1],
                'vsini_GF': res['vsgof'],
                'vsini_GF_eDW': res['evsgof'][0],
                'vsini_GF_eUP': res['evsgof'][1],
                'vmac_GF': res['vmgof'],
                'vmac_GF_eDW': res['evmgof'][0],
                'vmac_GF_eUP': res['evmgof'][1],
                'line_IB': res['line'],
                'EW_IB': res['ew'],
                'SNR_IB': res['snr'],
                'QIB': qib
            }
        )

        if check_best and res['filename'] != spec(res_id, snr='bestHF').filename:
            msg.warn(f"\n{res['filename']} does not match with best spectrum available.")

        bar.update(i)

    bar.finish()

    data_rows = Table(data_rows)

    # Post calculation modifications:
        # No values above 130 km/s were found in Simon-Diaz (2017)
        # we prevent the parameter to take values too large.
    data_rows['vmac_GF'][data_rows['vmac_GF'] > 130] = 130
    data_rows['SNR_IB'] = [int(round(row)) for row in data_rows['SNR_IB']]

    # Saving the results:
    full_path = maindir + 'tables/IB_results_new.' + format
    if format == 'ascii':
        format += '.fixed_width_two_line'

    data_rows.write(full_path, format=format, overwrite=True)

    return 'DONE'
