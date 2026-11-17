from __future__ import annotations

# ============================================================================
# ============================== LIST OF MODELS ==============================
# ============================================================================

# ==============================
def _read_models(
                 file_models : str | Path,
                 teff        : tuple[float, float],
                 grav        : tuple[float, float],
                 yhe         : tuple[float, float],
                 q           : tuple[float, float],
                 beta        : tuple[float, float]
    ) -> dict[str, np.ndarray]:
    """
    Read the list of models of the grid and select those within the given
    ranges of parameters (part of gbat_chi_global.pro).

    ----------
    Parameters
    ----------
        file_models : str | Path
            File 'iac_gbatgrid_HHe_<metal>_allHe.lis', one model per line, e.g.
            'iac_gbatgrid_HHe_Z100_He08_xdr/T320g390He08Q1300b10'.

        teff, grav, yhe, q, beta : tuple[float, float]
            Ranges of Teff [K], log g, Y(He), -log Q and beta.

    -------
    Returns
    -------
        models : dict[str, np.ndarray]
            red (sub-folder, with the final '/'), model (name), and teff [kK],
            grav, he, q, beta read from the name of the models.

    -----
    Notes
    -----
        - As in IDL, the parameters are read at fixed positions of the name
          (T320g390He08Q1300b10 -> Teff = 32.0 kK, log g = 3.90, Y(He) = 0.08,
          -log Q = 13.0, beta = 1.0), the position of the sub-folder being
          that of the first model.
    """

    # #################### MODULES & PACKAGES ####################
    import numpy as np


    # #################### PANORAMA ####################
    from gbat.io import _readcol
    from gbat.remote import _resolve


    # #################### READING ####################
    local = _resolve(file_models) # Downloaded first if the grid is in a share
    if local is None:
        raise FileNotFoundError(f"List of models of the grid not found: '{file_models}'. "
                                "Check CASA and METAL in gbat.variables.")
    (model,) = _readcol(local, fmt='A')
    if model.size == 0:
        raise ValueError(f"No models found in '{file_models}'.")
    hh = str(model[0]).find('/') # Position of the sub-folder separator (IDL: strpos of the first model)

    xred = np.array([m[:hh + 1] for m in model])
    xmodel = np.array([m[hh + 1:hh + 1 + 20] for m in model])
    num = lambda a, b, f: np.array([float(m[hh + 1 + a:hh + 1 + a + b]) * f for m in model])
    xteff = num(1, 3, 0.1)
    xgrav = num(5, 3, 0.01)
    xhe = num(10, 2, 0.01)
    xq = num(13, 3, 0.1)
    xbeta = num(18, 2, 0.1)


    # #################### SELECTION ####################
    nn = ((xteff * 1e3 >= teff[0]) & (xteff * 1e3 <= teff[1]) &
          (xgrav >= grav[0]) & (xgrav <= grav[1]) &
          (xhe >= yhe[0]) & (xhe <= yhe[1]) &
          (xq >= q[0]) & (xq <= q[1]) &
          (xbeta >= beta[0]) & (xbeta <= beta[1]))


    # #################### OUTPUT ####################
    return dict(red=xred[nn], model=xmodel[nn], teff=xteff[nn], grav=xgrav[nn], he=xhe[nn], q=xq[nn], beta=xbeta[nn])




# ===================================================================
# ============================== MODELS ==============================
# ===================================================================

# ==============================
def _model_file(base : str | Path) -> Path | None:
    """
    File of a FASTWIND model: '<base>.xdr' (IDL save file of the grid), or
    else '<base>.idl' (pickle). None if neither exists. If the grid is in a
    share (base is a link), the file is downloaded to the cache first.
    """

    # #################### PANORAMA ####################
    from gbat.remote import _resolve


    # #################### FILE ####################
    for ext in ('.xdr', '.idl'):
        file = _resolve(str(base) + ext)
        if file is not None:
            return file


    # #################### OUTPUT ####################
    return None




# ==============================
def _get(obj, name : str):
    """
    Field of a structure restored from an IDL file (record array, scalar
    value taken for 1-element structures) or of a dict.
    """

    # #################### MODULES & PACKAGES ####################
    import numpy as np


    # #################### FIELD ####################
    if isinstance(obj, dict):
        value = obj[name.lower()] if name.lower() in obj else obj[name.upper()]
    else:
        value = obj[name.lower()]
        if isinstance(value, np.ndarray) and value.ndim and value.size == 1 and value.dtype == object:
            value = value[0]
        elif isinstance(value, np.ndarray) and value.shape == (1,):
            value = value[0]
    if isinstance(value, bytes):
        value = value.decode()


    # #################### OUTPUT ####################
    return value




# ==============================
def _find_model(
                cfg  : Config,
                name : str
    ) -> Path | None:
    """
    File of a model of the grid given only its name (e.g. the reference
    model fw_ref, 'T370g390He10Q1350b10'), looking for its sub-folder in the
    list of models (IDL: findfile(casa_red+'*/'+fw_ref+'.xdr')).
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PANORAMA ##########
    from gbat.io import _readcol
    from gbat.remote import _resolve


    # #################### LIST OF MODELS ####################
    local = _resolve(cfg.file_models)
    if local is None:
        return None
    (model,) = _readcol(local, fmt='A')


    # #################### OUTPUT ####################
    for m in model:
        if m.rsplit('/', 1)[-1] == name:
            return _model_file(cfg.casa_red + m)
    return None




# ==============================
_MODEL_CACHE: dict = {}




# ==============================
def _load_model(file : str | Path) -> dict:
    """
    Read a FASTWIND model of the grid (IDL: restore of the model .xdr file).

    ----------
    Parameters
    ----------
        file : str | Path
            IDL save file (.xdr) or pickle file (.idl) of the model, with the
            variables
                >> m    : structure with teff [K], logg, yhe, q (log Q), beta
                >> mag  : structure with the synthetic magnitudes v, ks
                >> vsyn : synthetic V magnitude 'ala Kudritzki'
                >> l    : structure array with the lines: linea (FASTWIND
                          label), micro (e.g. '005'), x (wavelength [A]) and
                          y (normalised flux).
            In a pickle file, m and mag may be dicts and l a list of dicts.

    -------
    Returns
    -------
        model : dict
            teff, logg, yhe, q, beta, vsyn, v, ks and lines, a dict
            {(linea, micro): (x, y)} with the profiles of the lines.

    -----
    Notes
    -----
        - The last models read are kept in memory, since step 2 reads the same
          model once per line.
        - When a line appears more than once with the same micro, the first
          one is used.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    from pathlib import Path
    import numpy as np


    # ########## PANORAMA ##########
    from idl_bridge import _restore_idl


    # #################### CACHE ####################
    key = str(Path(file).resolve())
    if key in _MODEL_CACHE:
        return _MODEL_CACHE[key]


    # #################### READING ####################
    data = _restore_idl(file)
    m, mag, l = data['m'], data['mag'], data['l']

    lines = {}
    if isinstance(l, (list, tuple)):
        records = [(str(_get(r, 'linea')), str(_get(r, 'micro')), _get(r, 'x'), _get(r, 'y')) for r in l]
    else:
        dec = lambda s: s.decode() if isinstance(s, bytes) else str(s)
        records = [(dec(a), dec(b), x, y) for a, b, x, y in zip(l['linea'], l['micro'], l['x'], l['y'])]
    for linea, micro, x, y in records:
        lines.setdefault((linea.strip(), micro.strip()), (np.asarray(x, dtype=float), np.asarray(y, dtype=float)))

    model = dict(teff=float(_get(m, 'teff')), logg=float(_get(m, 'logg')), yhe=float(_get(m, 'yhe')),
                 q=float(_get(m, 'q')), beta=float(_get(m, 'beta')), vsyn=float(np.ravel(data['vsyn'])[0]),
                 v=float(_get(mag, 'v')), ks=float(_get(mag, 'ks')), lines=lines)


    # #################### CACHE ####################
    if len(_MODEL_CACHE) > 8:
        _MODEL_CACHE.clear()
    _MODEL_CACHE[key] = model


    # #################### OUTPUT ####################
    return model
