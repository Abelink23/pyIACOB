from __future__ import annotations


# ========================================================================
# ===================== SAVE FILES OR RESTORE SOLUTIONS ==================
# ========================================================================

# ==============================
def _to_json(value):
    """
    Convert a saved variable into a JSON-compatible object.

    ----------
    Parameters
    ----------
        value : any
            Variable to convert (NumPy array or scalar, dict, list, str...).

    -------
    Returns
    -------
        json_value : any
            JSON-compatible object.

    -----
    Notes
    -----
        - NumPy arrays become lists and NumPy scalars become Python numbers.
        - Dicts, lists and tuples are converted recursively.

        - Records of IDL structures (e.g. d restored from an IDL file) become
          dicts of tags.

        - Non-finite values (e.g. the Inf of the FT at zero frequency) become
          null, because JSON does not accept them.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # #################### CONVERSION ####################
    # ########## CONTAINERS ##########
    if isinstance(value, dict):
        return {k: _to_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_to_json(v) for v in value]


    # ########## ARRAYS ##########
    if isinstance(value, np.ndarray):
        if value.dtype.names: # IDL structure -> dict of tags
            return {n.lower(): _to_json(value[n][0] if value.size == 1 else value[n]) for n in value.dtype.names}
        return [_to_json(v) for v in value.tolist()] if value.ndim else _to_json(value.item())


    # ########## SCALARS ##########
    if isinstance(value, np.generic): # NumPy scalar -> Python number
        value = value.item()
    if isinstance(value, bytes): # IDL strings are read as bytes
        return value.decode()
    if isinstance(value, float) and not np.isfinite(value): # Inf/NaN -> null
        return None


    # #################### OUTPUT ####################
    return value




# ==============================
def _save(
          output_file : str | Path,
          readable    : bool = True,
          output_dir  : str | Path | None = None,
          **variables
    ) -> None:
    """
    Save variables to disk (replacement of IDL SAVE).

    ----------
    Parameters
    ----------
        output_file : str | Path
            Output file (e.g. HD37042_SiIII.idl). It is a pickle file, read
            back by _restore.

        readable : bool, optional
            Also write a human-readable copy with the same name and extension
            .json (e.g. HD37042_SiIII.json), which can be opened with any
            text editor or read with json.load / pandas. Defaults to True.

        output_dir : str | Path | None, optional
            Folder of the .json copy. If None, the folder of 'output_file'.
            Defaults to None.

        **variables
            Variables to save. Names are stored in lower case, as IDL does.

    -------
    Returns
    -------
        None

    -----
    Notes
    -----
        - The .idl (pickle) file keeps the exact NumPy types and is the one used
          by IACOB-BROAD. The .json copy is only for inspection and exchange.

        - IDL XDR save files cannot be written from Python, so the .idl files
          written here cannot be read by the IDL version.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import json
    import pickle
    from pathlib import Path


    # #################### VARIABLES ####################
    data = {k.lower(): v for k, v in variables.items()} # Lower-case names, as IDL


    # #################### BINARY FILE ####################
    with open(output_file, 'wb') as fh: # Used by _restore
        pickle.dump(data, fh, protocol=4)


    # #################### READABLE COPY ####################
    if readable:
        json_file = Path(output_dir or Path(output_file).parent) / (Path(output_file).stem + '.json')
        with open(json_file, 'w') as fh:
            json.dump(_to_json(data), fh, indent=1)




# ==============================
def _restore_idl(
             file : str | Path
    ) -> dict:
    """
    Read the variables of a save file (replacement of IDL RESTORE).

    ----------
    Parameters
    ----------
        file : str | Path
            File written by _save (pickle) or by the original IDL
            IACOB-BROAD (XDR save file).

    -------
    Returns
    -------
        data : dict
            Saved variables, with lower-case names.

    -----
    Notes
    -----
        - IDL save files are recognised by their signature 'SR' and read with
          scipy.io.readsav, so lines pre-processed with the IDL version can be
          used directly.

        - Single-precision values of IDL files are converted to double
          precision, the working precision of the package.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import pickle

    import numpy as np


    # #################### FILE FORMAT ####################
    with open(file, 'rb') as fh:
        is_idl = fh.read(2) == b'SR' # IDL save files start with 'SR'


    # #################### PICKLE FILE ####################
    if not is_idl:
        with open(file, 'rb') as fh:
            return pickle.load(fh)


    # #################### IDL FILE ####################
    from scipy.io import readsav
    data = {}
    for key, val in readsav(str(file)).items():
        if isinstance(val, np.ndarray) and val.dtype.kind == 'f': # Float arrays -> native double precision
            val = val.astype(float)
        elif isinstance(val, np.floating): # Float scalars -> Python float
            val = float(val)
        data[key.lower()] = val

    # if vmin is not provided, set it to 0.0
    if 'vmin' not in data:
        print(f"Warning: 'vmin' not found in {file}. Setting 'vmin' to 0.0.")
        data['vmin'] = 0.0

    # #################### OUTPUT ####################
    return data


def _restore_json(
                json_file : str | Path
        ) -> dict:
    """
    Restore variables saved in a JSON file by _save(), converting numeric
    lists back to NumPy arrays so they can be consumed by GOF and plotting routines.
    """
    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import json

    import numpy as np


    with open(json_file, 'r') as fh:
        raw = json.load(fh)

    def _convert(val):
        if isinstance(val, dict):
            return {k.lower(): _convert(v) for k, v in val.items()}
        if isinstance(val, list):
            if len(val) == 0:
                return np.array([])
            # Numerical 1D arrays or lists containing numbers / nulls
            if all(isinstance(x, (int, float, type(None))) for x in val):
                clean = [np.nan if x is None else x for x in val]
                if all(isinstance(x, int) for x in val):
                    return np.array(clean, dtype=int)
                return np.array(clean, dtype=float)
            # Nested structures
            converted = [_convert(x) for x in val]
            try:
                arr = np.array(converted)
                if arr.dtype.kind in ('i', 'f', 'u'):
                    return arr
            except Exception:
                pass
            return converted
        return val

    return {k.lower(): _convert(v) for k, v in raw.items()}




# =============================================================================
# ============================== TERMINAL OUTPUT ==============================
# =============================================================================

# ==============================
def _print_results(
                   results      : dict,
                   user_profile : bool
    ) -> None:
    """
    Print the final table of results of IACOB-BROAD.

    ----------
    Parameters
    ----------
        results : dict
            Summary of results (see pyiacob_broad).

        user_profile : bool
            True if an intrinsic profile of the user was used, False for
            the delta-like profile.

    -------
    Returns
    -------
        None

    -----
    Notes
    -----
        - EW and SNR are truncated to integers, as the IDL I format.

        - Columns: vs0 = vsini(GOF) with vmac = 0, vft = vsini(FT),
          vmft = vmac(GOF) for vsini = vft, vsgof / vmgof = vsini / vmac (GOF).

        - Values with errors are written as "value - lower + upper".
    """

    # #################### FORMAT ####################
    def _pm(v, e): # "value - lower + upper"
        return f"{v:.1f} - {e[0]:.1f} + {e[1]:.1f}"


    # #################### TABLE ####################
    """
    - Column names are the keys of 'results' (units in brackets), so the
      values are easy to find in the returned dictionary.
    - Every column is as wide as its longest entry, and right-aligned.
    """
    columns = [
               ('filename',    str(results['filename'])),
               ('line',        str(results['line'])),
               ('ew [mA]',     f"{int(results['ew'])}"),
               ('snr',         f"{int(results['snr'])}"),
               ('vs0 [km/s]',  f"{results['vs0']:.1f}"),
               ('vft [km/s]',  f"{results['vft']:.1f}"),
               ('vmft [km/s]', _pm(results['vmft'], results['evmft'])),
               ('vsgof [km/s]', _pm(results['vsgof'], results['evsgof'])),
               ('vmgof [km/s]', _pm(results['vmgof'], results['evmgof'])),
               ('QIB',         f"{results.get('qib', 0.0):.1f}"),
               ('profile',     'USER PROFILE' if user_profile else 'DELTA FUNCTION'),
              ]
    widths = [max(len(name), len(value)) for name, value in columns]
    print('\n')
    print('  '.join(f"{name:>{w}}" for (name, _), w in zip(columns, widths)))
    print('  '.join('-' * w for w in widths))
    print('  '.join(f"{value:>{w}}" for (_, value), w in zip(columns, widths)))




# ============================================================================
# ============================== SUMMARY FIGURE ==============================
# ============================================================================

# ==============================
def _panel_positions(
                     m    : tuple[int, ...] = (2, 2, 2),
                     xgap : float = 0.10,
                     ygap : float = 0.07,
                     xmin : float = 0.1,
                     xmax : float = 0.95,
                     ymin : float = 0.4,
                     ymax : float = 0.96
    ) -> list[list[float]]:
    """
    Positions of panels arranged in rows (translation of crea_subplot.pro).

    ----------
    Parameters
    ----------
        m : tuple[int, ...], optional
            Number of panels in each row, from top to bottom.
            Defaults to (2, 2, 2).

        xgap : float, optional
            Horizontal gap between panels (figure fraction).

        ygap : float, optional
            Vertical gap between panels (figure fraction).

        xmin : float, optional
            Lower horizontal limit of the plotting area (figure fraction).

        xmax : float, optional
            Upper horizontal limit of the plotting area (figure fraction).

        ymin : float, optional
            Lower vertical limit of the plotting area (figure fraction).

        ymax : float, optional
            Upper vertical limit of the plotting area (figure fraction).

    -------
    Returns
    -------
        positions : list[list[float]]
            [left, bottom, width, height] of each panel, row by row from the
            top-left one.
    """

    # #################### PANEL POSITIONS ####################
    height = (ymax - ymin - ygap * (len(m) - 1)) / len(m) # Height of each row
    positions = []
    for i, ncol in enumerate(m):
        width = (xmax - xmin - xgap * (ncol - 1)) / ncol # Width of each panel of row i
        bottom = ymax - (i + 1) * height - i * ygap # Bottom of row i
        positions += [[xmin + j * (width + xgap), bottom, width, height] for j in range(ncol)]


    # #################### OUTPUT ####################
    return positions




# ==============================
def _plot_map(
              ax   : plt.Axes,
              chi2 : np.ndarray[float],
              vs   : np.ndarray[float],
              vrt  : np.ndarray[float],
              vsf  : float
    ) -> None:
    """
    Filled contours of chi2 - min(chi2) in the (vsini, vmac) plane
    (translation of plot_map).

    ----------
    Parameters
    ----------
        ax : matplotlib.axes.Axes
            Axis of the panel.

        chi2 : np.ndarray[float]
            chi2 of each model.

        vs : np.ndarray[float]
            vsini of each model [km/s].

        vrt : np.ndarray[float]
            vmac of each model [km/s].

        vsf : float
            Maximum vsini of the grid [km/s].

    -------
    Returns
    -------
        None

    -----
    Notes
    -----
        - Levels at min + [0, 1, 4, 9] (1, 2 and 3 sigma), filled with the
          ILUPAL greys 9, 10, 11 and white above the last level, as the IDL
          CONTOUR with /IRREGULAR and CELL_FILL.
        - The vmac range is set by the models within 3 sigma.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # ########## PANORAMA ##########
    from .interactive import _COLORS


    # #################### RANGES ####################
    chi0 = chi2 - chi2.min() # chi2 relative to the minimum
    near = vrt[chi2 <= chi2.min() + 9] # vmac of the models within 3 sigma
    yr = [-0.05 * (near.max() + near.min()), 1.05 * (near.max() + near.min())]
    xr = [-0.05 * vsf, 1.05 * vsf]


    # #################### CONTOURS ####################
    levels = chi0.min() + np.array([0, 1, 4, 9]) # 0, 1, 2 and 3 sigma
    ax.tricontourf(vs, vrt, chi0, levels=[*levels, np.inf], colors=[_COLORS[9], _COLORS[10], _COLORS[11], _COLORS[1]])
    cs = ax.tricontour(vs, vrt, chi0, levels=levels, colors='black', linewidths=1)
    ax.clabel(cs, levels=levels[1:], fmt=dict(zip(levels[1:], [r'1$\sigma$', r'2$\sigma$', r'3$\sigma$'])),
              fontsize=8)


    # #################### GRID LIMITS & LABELS ####################
    ax.plot([vs.min()] * 2, yr, 'k--') # Minimum vsini of the grid
    ax.plot(xr, [vrt.min()] * 2, 'k--') # Minimum vmac of the grid
    ax.set(xlim=xr, ylim=yr, title='Goodness-of-fit', xlabel=r'v sini [km s$^{-1}$]',
           ylabel=r'v$_{\rm mac}$ (R-T) [km s$^{-1}$]')




# ==============================
def _plot_summary(
                  output_file : str,
                  filename    : str,
                  line        : str,
                  prep        : dict,
                  gof         : dict,
                  est         : Estimates,
                  models      : dict,
                  fts         : dict,
                  env_vs      : Envelope,
                  env_vm      : Envelope,
                  vsf         : float,
                  gauss       : bool = False,
                  vs_user     : float = 0.,
                  vm_user     : float = 0.
    ) -> None:
    """
    Summary PDF figure of IACOB-BROAD.

    Portrait page with 3 rows x 2 panels:

        - top-left: observed profile and synthetic profiles of the solutions,
        - top-right: Fourier transforms of the observed and synthetic profiles,
        - middle-left: text summary of the results,
        - middle-right: chi2 vs vsini,
        - bottom-left: chi2 vs vmac,
        - bottom-right: chi2 map in the (vsini, vmac) plane.

    ----------
    Parameters
    ----------
        output_file : str
            Output PDF file.

        filename : str
            Filename of the spectrum.

        line : str
            Line label of the summary.

        prep : dict
            Pre-processed line (see _vsiniFT).

        gof : dict
            GOF results (see _gof_grid).

        est : Estimates
            Final estimates (see _final_estimates).

        models : dict
            Synthetic profiles (see _best_models).

        fts : dict
            Fourier transforms of the synthetic profiles (see _best_models).

        env_vs : Envelope
            Spline envelope of chi2 vs vsini.

        env_vm : Envelope
            Spline envelope of chi2 vs vmac.

        vsf : float
            Maximum vsini of the grid [km/s].

        gauss : bool, optional
            Gaussian macroturbulence (only changes the labels).

        vs_user : float, optional
            vsini of the user [km/s] (0 -> not drawn).

        vm_user : float, optional
            vmac of the user [km/s] (0 -> not drawn).

    -------
    Returns
    -------
        None

    -----
    Notes
    -----
        - Colours follow the ILUPAL table of IDL: red for vsini(FT), blue for
          the 2D GOF solution, green for vsini(GOF, vmac = 0), purple for
          vsini(FT) + vmac(GOF) and pink for the user values.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import matplotlib.pyplot as plt
    import numpy as np


    # ########## PANORAMA ##########
    from .interactive import _COLORS, _DASH_DOT3


    # #################### FIGURE ####################
    plt.rcParams.update({'font.size': 9})
    fig = plt.figure(figsize=(20. / 2.54, 27. / 2.54)) # Portrait page of 20 x 27 cm
    pos = _panel_positions()
    wx0, fx0, snr, vft = prep['wx0'], prep['fx0'], prep['snr'], prep['vft']
    chi2, vs, vrt = gof['chi2'], gof['chi2_vs'], gof['chi2_vrt']
    cmin = chi2.min()
    user = vs_user > 5. # User values drawn in the GOF panels only with rotation
    vmac_lab = r'v$_{\rm mac,G}$)' if gauss else r'v$_{\rm mac,RT}$' # IDL labels (the stray ')' is in the IDL code)


    # #################### LINE PROFILE ####################
    ax = fig.add_axes(pos[0])
    span = wx0.max() - wx0.min()
    ax.set(xlim=[wx0.min() - 0.3 * span, wx0.max() + 0.3 * span], ylim=[0.98 * fx0.min(), 1.02 * fx0.max()],
           xlabel=r'$\lambda$ [A]', ylabel='Norm. flux')
    ax.xaxis.set_major_locator(plt.MaxNLocator(3))


    # ########## OBSERVATIONS ##########
    ax.plot([prep['w0'].min(), prep['w0'].max()], [1. - 3. / snr] * 2, color=_COLORS[8], ls='--') # 3-sigma noise level
    ax.plot(*models['inst'], color=_COLORS[10], lw=3) # Intrinsic profile
    ax.plot(prep['w0'], prep['f0'], color=_COLORS[9], lw=3) # Full observed spectrum
    ax.errorbar(wx0, fx0, yerr=np.full(wx0.size, 1. / snr), fmt='o', ms=1.5, color='black', elinewidth=0.4) # Line points with 1/SNR error bars


    # ########## SYNTHETIC PROFILES ##########
    styles = {'ft': dict(color=_COLORS[2], ls='--'), 'gof0': dict(color=_COLORS[4]),
              'gof': dict(color=_COLORS[3]), 'ftgof': dict(color=_COLORS[7], ls='--'), 'user': dict(color=_COLORS[6])}
    for key, style in styles.items():
        if key in models:
            ax.plot(*models[key], **style)


    # #################### TEXT SUMMARY ####################
    ax = fig.add_axes(pos[2])
    ax.axis('off')
    ax.set(xlim=[0, 1], ylim=[0, 1])
    text = [(1.10, f'Spectrum: {filename}', 0),
            (0.95, f'Line: {line}   EW={est.ew:.1f} mA   SNR={int(snr)}', 0),
            (0.70, f'FT:  vsini = {vft:.1f} km/s  ({vmac_lab}=0)', 2),
            (0.50, f'GOF: vsini = {est.vsgof:.1f} km/s , {vmac_lab} = {est.vmgof:.1f} km/s', 3),
            (0.35, f'GOF: vsini ({vmac_lab}=0) = {est.vs0:.1f} km/s', 4),
            (0.20, f'GOF: {vmac_lab} (vsini(FT) = {vft:.1f}) = {est.vmft:.1f} km/s', 7),
            (-0.20, 'GOF: INTR. PROFILE -> ' + ('USER' if 'w_intprof' in prep else r'$\delta$ FUNCTION'), 0),
            (-0.30, rf'$\Delta$vsini= {gof["avs"]:.2f} km/s ,  $\Delta${vmac_lab} = {gof["avs"]:.2f} km/s ', 0)]
    if vs_user:
        text.append((0.00, f'USER: vsini = {vs_user:.1f} km/s , {vmac_lab} = {vm_user:.1f} km/s', 6))
    for y, label, c in text:
        ax.text(-0.1, y, label, color=_COLORS[c], fontsize=9, clip_on=False)


    # #################### FOURIER TRANSFORMS ####################
    ax = fig.add_axes(pos[1])
    ax.set(xlim=[-0.05 * vsf, 1.05 * vsf], ylim=[-4., 0.], xlabel=r'v sini [km s$^{-1}$]',
           ylabel='ABS (FT AMPLITUDE)', title='Fourier transform')
    ax.plot([prep['max_sigma_obs']] * 2, [-8, -0.7], color='black', ls=_DASH_DOT3) # Nyquist limit
    ax.text(prep['max_sigma_obs'], -0.6, 'Nyq.', fontsize=6, rotation=90)
    if 'inst' in fts: # FT of the intrinsic profile of the user
        ax.plot(fts['inst'].vsini, fts['inst'].log_amp, color=_COLORS[10], lw=2.5)
    ax.plot(prep['xfou_obs'], prep['yfou_obs'], color='black') # FT of the observed line
    ax.axvline(vft, color=_COLORS[2], ls='--') # vsini(FT)
    ft_styles = {'ft': dict(color=_COLORS[2], ls='--'), 'gof': dict(color=_COLORS[3]),
                 'ftgof': dict(color=_COLORS[7], ls='-.'), 'user': dict(color=_COLORS[6])}
    for key, style in ft_styles.items():
        if key in fts:
            ax.plot(fts[key].vsini, fts[key].log_amp, **style)


    # #################### CHI2 MAP ####################
    ax = fig.add_axes(pos[5])
    _plot_map(ax, chi2, vs, vrt, vsf)
    ax.plot(est.vsgof, est.vmgof, 's', mfc='none', mec=_COLORS[3]) # 2D GOF solution
    ax.plot(vft, est.vmft, 's', mfc='none', mec=_COLORS[2]) # vsini(FT) + vmac(GOF)
    if user:
        ax.plot(vs_user, vm_user, 's', mfc='none', mec=_COLORS[6])


    # #################### CHI2 VS VSINI ####################
    ax = fig.add_axes(pos[3])
    ax.set(xlim=[-0.05 * vsf, 1.05 * vsf], ylim=cmin + np.array([-1., 4.5]), xlabel=r'vsini [km s$^{-1}$]',
           ylabel=r'$\chi^2$')
    ax.plot(vs, chi2, 'o', ms=4, color=_COLORS[10]) # All the models
    ax.plot([vs.min()] * 2, cmin + np.array([4.5, -1]), 'k--') # Minimum vsini of the grid
    ax.plot([est.vsgof] * 2, cmin + np.array([-1, -0.02]), color=_COLORS[3], lw=2) # vsini(GOF)
    ax.plot([vft] * 2, cmin + np.array([-1, -0.02]), color=_COLORS[2], ls='-.', lw=2) # vsini(FT)
    if user:
        ax.plot([vs_user] * 2, cmin + np.array([-1, -0.02]), color=_COLORS[6], ls='--', lw=2)
    ax.plot([0, vsf], [cmin + 1] * 2, 'k--') # 1 sigma
    ax.plot([0, vsf], [cmin + 4] * 2, 'k-.') # 2 sigma
    ax.plot(env_vs.x, env_vs.y, color='black', lw=1.5) # Spline envelope
    ax.text(0.05 * vsf, cmin - 0.4, rf'$\sigma$ (BFM) = {gof["sigma_chi"].min():.3f}', fontsize=7) # Rms of the best model
    ax.text(0.05 * vsf, cmin - 0.8, rf'$\sigma$ (SNR) = {1. / snr:.3f}', fontsize=7) # Noise of the continuum


    # #################### CHI2 VS VMAC ####################
    ax = fig.add_axes(pos[4])
    near = vrt[chi2 <= cmin + 16] # vmac of the models within 4 sigma
    vmax = near.max() + near.min()
    ax.set(xlim=cmin + np.array([4.5, -1]), ylim=[-0.05 * vmax, 1.05 * vmax], xlabel=r'$\chi^2$',
           ylabel=r'v$_{\rm mac}$ (' + ('G' if gauss else 'R-T') + r') [km s$^{-1}$]')
    ax.plot(chi2, vrt, 'o', ms=4, color=_COLORS[10]) # All the models (chi2 on the x-axis, decreasing)
    ax.plot(cmin + np.array([4.5, -1]), [vrt.min()] * 2, 'k--') # Minimum vmac of the grid
    ax.plot(cmin + np.array([-1, -0.02]), [est.vmgof] * 2, color=_COLORS[3], lw=2) # vmac(GOF)
    ax.plot(cmin + np.array([-1, -0.02]), [est.vmft] * 2, color=_COLORS[7], ls='-.', lw=2) # vmac(GOF, vsini(FT))
    if user:
        ax.plot(cmin + np.array([-1, -0.02]), [vm_user] * 2, color=_COLORS[6], ls='--', lw=2)
    ax.plot([cmin + 1] * 2, [vrt.min(), vrt.max()], 'k--') # 1 sigma
    ax.plot([cmin + 4] * 2, [vrt.min(), vrt.max()], 'k-.') # 2 sigma
    ax.plot(env_vm.y, env_vm.x, color='black', lw=1.5) # Spline envelope


    # #################### SAVE FIGURE ####################
    fig.savefig(output_file, format='pdf')
    plt.close(fig)
