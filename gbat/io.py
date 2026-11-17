from __future__ import annotations

"""
iacob_gbat_io — Input/output of IACOB-GBAT.

    - IDL formatted output (A, I, F codes and free format), so that the
      terminal output and the text files keep the IDL formats.
    - Reading of tables (replacement of readcol.pro and the auxiliary
      gettok, remchar, strnumber and numlines of the IDL astro library).
    - Reading/writing of the observed and synthetic spectra.
    - IDL structures (sol, dat_line, resul, info) as NumPy structured arrays.
    - Save/restore of results (pickle .idl + .json, IDL .xdr read as fallback).
    - Conversion of IDL text-formatting commands ('!7a!3', '!dt!n' ...) into
      Matplotlib mathtext.
"""

# ====================================================================================
# ============================== IDL FORMATTED OUTPUT ==============================
# ====================================================================================

# ==============================
def _idl_round(x : float | np.ndarray) -> float | np.ndarray:
    """
    Round half away from zero (IDL ROUND), unlike the round-half-to-even
    rule of numpy.round.

    ----------
    Parameters
    ----------
        x : float | np.ndarray
            Value(s) to round.

    -------
    Returns
    -------
        rounded : float | np.ndarray
            Rounded value(s), as floats.
    """

    # #################### MODULES & PACKAGES ####################
    import numpy as np


    # #################### ROUND ####################
    rounded = np.sign(x) * np.floor(np.abs(x) + 0.5)


    # #################### OUTPUT ####################
    return rounded




# ==============================
def _fa(
        s : str,
        w : int | None = None
    ) -> str:
    """
    IDL 'Aw' output: string right-justified in w characters, or truncated to
    its first w characters if longer. Without w, the string is unchanged.
    """

    # #################### FORMAT ####################
    s = str(s)
    if w is None:
        return s


    # #################### OUTPUT ####################
    return s[:w] if len(s) > w else s.rjust(w)




# ==============================
def _fi(
        x : float,
        w : int
    ) -> str:
    """
    IDL 'Iw' output: integer right-justified in w characters, '*' * w if it
    does not fit.

    -----
    Notes
    -----
        - Floating-point values are rounded to the nearest integer (half away
          from zero), as the IDL formatted output does (e.g. 1/0.00666667 =
          149.99993 is written as 150).
        - Non-finite values are written as '*' * w.
    """

    # #################### MODULES & PACKAGES ####################
    import numpy as np


    # #################### FORMAT ####################
    if not np.isfinite(x):
        return '*' * w
    s = f'{int(_idl_round(float(x))):{w}d}'


    # #################### OUTPUT ####################
    return '*' * w if len(s) > w else s




# ==============================
def _ff(
        x : float,
        w : int,
        d : int
    ) -> str:
    """
    IDL 'Fw.d' output: fixed-point number in w characters with d decimals,
    '*' * w if it does not fit.
    """

    # #################### MODULES & PACKAGES ####################
    import numpy as np


    # #################### FORMAT ####################
    x = float(x)
    s = f'{x:{w}.{d}f}' if np.isfinite(x) else f'{x:{w}}'


    # #################### OUTPUT ####################
    return '*' * w if len(s) > w else s




# ==============================
def _fg(
        x      : float,
        double : bool = False
    ) -> str:
    """
    IDL free-format output of a number (PRINT/PRINTF without FORMAT, or
    STRING without FORMAT).

    ----------
    Parameters
    ----------
        x : float
            Value to write.

        double : bool, optional
            True for IDL double precision (G16.8-like), False for single
            precision (G13.6-like). Defaults to False.

    -------
    Returns
    -------
        s : str
            Formatted value, e.g. '      3763.78', '      0.00000' or
            '  1.00000e-05' (single precision).

    -----
    Notes
    -----
        - The C format '%#w.pg' keeps the trailing zeros, as IDL does.
    """

    # #################### FORMAT ####################
    s = ('%#16.8g' if double else '%#13.6g') % float(x)


    # #################### OUTPUT ####################
    return s




# ==============================
def _fstr(
          x      : float,
          double : bool = False
    ) -> str:
    """
    IDL strtrim(string(x),2): free-format number without blanks.
    """

    # #################### OUTPUT ####################
    return _fg(x, double).strip()




# ===========================================================================
# ============================== TABLE READING ==============================
# ===========================================================================

# ==============================
def _strnumber(s : str) -> tuple[bool, float]:
    """
    Check whether a string is a valid number (IDL strnumber.pro).

    -------
    Returns
    -------
        valid : bool
            True if the string can be read as a number.

        value : float
            The number (nan if not valid).

    -----
    Notes
    -----
        - IDL double-precision exponents (1.0d0, 2D-3) are accepted.
    """

    # #################### CONVERSION ####################
    try:
        return True, float(s.strip().replace('d', 'e').replace('D', 'E'))
    except ValueError:
        return False, float('nan')




# ==============================
def _readcol(
             file : str | Path,
             fmt  : str | None = None,
             ncol : int | None = None
    ) -> list[np.ndarray]:
    """
    Read the columns of a free-format text file (replacement of readcol.pro).

    ----------
    Parameters
    ----------
        file : str | Path
            Text file.

        fmt : str | None, optional
            IDL format of the columns, e.g. '(A,F,F)': A (string), F/D
            (float), I/L/B (integer) or X (column skipped). If None, ncol
            columns of floats are read.

        ncol : int | None, optional
            Number of columns when fmt is None.

    -------
    Returns
    -------
        columns : list[np.ndarray]
            One array per column of the format (X columns excluded).

    -----
    Notes
    -----
        - As in readcol.pro, the blanks of each line are compressed (tabs
          included) and the line is split at blanks and commas. Lines with
          fewer fields than the format, or with a non-numeric value in a
          numeric column, are skipped. Extra fields are ignored.
        - There is no comment character (as the IDL calls of IACOB-GBAT):
          comment lines are skipped because they are not numeric, or because
          they do not have enough fields.
    """

    # #################### MODULES & PACKAGES ####################
    import re
    import numpy as np


    # #################### FORMAT ####################
    if fmt is None:
        codes = ['F'] * (ncol or 1)
    else:
        codes = [c[:1] if c else 'F' for c in re.sub(r'[\s()]', '', fmt.upper()).split(',')]
    keep = [i for i, c in enumerate(codes) if c != 'X'] # Columns returned
    nfmt = len(codes)


    # #################### READING ####################
    rows = []
    with open(file, 'r', errors='replace') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if len(line) < len(keep): # At least 1 character per output column
                continue
            fields = re.split(r'[ ,]', ' '.join(line.strip().split())) # strcompress + strsplit(' ,', /preserve_null)
            if len(fields) < nfmt:
                continue
            row, valid = [], True
            for i in keep:
                if codes[i] == 'A':
                    row.append(fields[i])
                else:
                    ok, val = _strnumber(fields[i])
                    if not ok:
                        valid = False
                        break
                    row.append(int(val) if codes[i] in 'ILB' else val)
            if valid:
                rows.append(row)


    # #################### COLUMNS ####################
    columns = []
    for k, i in enumerate(keep):
        values = [r[k] for r in rows]
        if codes[i] == 'A':
            columns.append(np.array(values, dtype=str))
        elif codes[i] in 'ILB':
            columns.append(np.array(values, dtype=int))
        else:
            columns.append(np.array(values, dtype=float))


    # #################### OUTPUT ####################
    return columns




# =====================================================================
# ============================== SPECTRA ==============================
# =====================================================================

# ==============================
def _write_spectrum_cut(
                        file : str | Path,
                        wo0  : np.ndarray[float],
                        fo0  : np.ndarray[float],
                        pe0  : np.ndarray[float],
                        sig0 : np.ndarray[float]
    ) -> None:
    """
    Write the pre-processed spectrum in spectra_cut (4 columns: wavelength,
    normalised flux, line mask and SNR mask), with the IDL free format of
    single-precision numbers (13 characters per column).
    """

    # #################### MODULES & PACKAGES ####################
    import numpy as np
    from pathlib import Path

    # Replace any original extension with '.dat'
    file = Path(file).with_suffix('.dat')

    # #################### WRITE ####################
    data = np.column_stack([wo0, fo0, pe0, sig0])
    np.savetxt(file, data, fmt='%#13.6g', delimiter='')




# ==============================
def _wascii(
            file : str | Path,
            w    : np.ndarray[float],
            f    : np.ndarray[float]
    ) -> None:
    """
    Write a 2-column spectrum with the IDL format '(D,D)', i.e. D25.16 for
    each column (translation of gbat_wascii.pro).
    """

    # #################### MODULES & PACKAGES ####################
    import numpy as np


    # #################### WRITE ####################
    np.savetxt(file, np.column_stack([w, f]), fmt='%25.16f', delimiter='')




# ==============================
def _read_spectrum(
                   file  : str | Path,
                   ncols : int
    ) -> list[np.ndarray[float]]:
    """
    Read the first ncols numeric columns of a spectrum (readcol with format
    F).
    """

    # #################### MODULES & PACKAGES ####################
    import numpy as np


    # #################### READING ####################
    """
    - Fast path: numpy.loadtxt on purely numeric files; otherwise the
      line-by-line rules of readcol.
    """
    try:
        data = np.loadtxt(file, usecols=range(4), ndmin=2, comments=None)
        columns = [data[:, i].copy() for i in range(4)]
    except ValueError:
        columns = _readcol(file, fmt='(' + ','.join(['F'] * 4) + ')')


    # #################### OUTPUT ####################
    return columns




# ========================================================================
# ============================== STRUCTURES ==============================
# ========================================================================

# ==============================
def _sol_init(
              nlines  : int,
              nmodels : int
    ) -> np.ndarray:
    """
    Structure 'sol' with the results of each model of the grid (translation
    of gbat_sol.pro), as a NumPy structured array.

    -----
    Notes
    -----
        - Fields: red (sub-folder of the model), model, teff [kK], grav,
          ehe, yhe, Q (-log Q, positive), beta, micro (string, e.g. '005'),
          V, Vsyn, Ks (synthetic magnitudes), R, L, M, mdot, chi_lin (chi2 of
          each line) and chi (global chi2). Names are lower case.
    """

    # #################### MODULES & PACKAGES ####################
    import numpy as np


    # #################### STRUCTURE ####################
    dtype = [('red', 'U128'), ('model', 'U64'),
             ('teff', 'f8'), ('grav', 'f8'), ('ehe', 'f8'), ('yhe', 'f8'),
             ('q', 'f8'), ('beta', 'f8'), ('micro', 'U8'),
             ('v', 'f8'), ('vsyn', 'f8'), ('ks', 'f8'),
             ('r', 'f8'), ('l', 'f8'), ('m', 'f8'), ('mdot', 'f8'),
             ('chi_lin', 'f8', (nlines,)), ('chi', 'f8')]


    # #################### OUTPUT ####################
    return np.zeros(nmodels, dtype=dtype)




# ==============================
def _dat_line_init(
                   nlines : int,
                   sigma  : float = 0.
    ) -> np.ndarray:
    """
    Structure 'dat_line' with the information of each line: line (FASTWIND
    label), sigma (noise of the line) and np (number of points of the line).
    """

    # #################### MODULES & PACKAGES ####################
    import numpy as np


    # #################### STRUCTURE ####################
    dat_line = np.zeros(nlines, dtype=[('line', 'U32'), ('sigma', 'f8'), ('np', 'f8')])
    dat_line['sigma'] = sigma


    # #################### OUTPUT ####################
    return dat_line




# ==============================
RESUL_FIELDS = ('title', 'var_min', 'var_max', 'bfm',
                'min_1x', 'max_1x', 'bfm_1x', 'errmin_1x', 'errmax_1x', 'err_1x',
                'min_2x', 'max_2x', 'bfm_2x', 'errmin_2x', 'errmax_2x', 'err_2x')




# ==============================
def _resul_init() -> np.ndarray:
    """
    Structure 'resul' with the results for one parameter (see
    gbat_plot_dist_one.pro), as a 1-element NumPy structured array. The
    field order is the IDL one (write_resul_one uses tag numbers).

    -----
    Notes
    -----
        - The fields are single precision, as in IDL. This matters: in
          write_resul_one the (float32) best value is compared with the
          (double) grid values var_x, and the rounding of resul.bfm can
          select the previous grid point (e.g. -log Mdot = 6.42 +/- 0.15 in
          the IDL benchmark, instead of +/- 0.07 with double precision).
    """

    # #################### MODULES & PACKAGES ####################
    import numpy as np


    # #################### OUTPUT ####################
    return np.zeros(1, dtype=[('title', 'U64')] + [(f, 'f4') for f in RESUL_FIELDS[1:]])




# ==============================
def _from_idl(value):
    """
    Convert a variable restored from an IDL save file (scipy.io.readsav) into
    the Python types of IACOB-GBAT.

    -----
    Notes
    -----
        - Strings (bytes) -> str.
        - Structures (record arrays) -> structured arrays with lower-case
          names, string fields as unicode and array fields (e.g. chi_lin)
          stacked into fixed-size sub-arrays.
        - Floats -> double precision.
    """

    # #################### MODULES & PACKAGES ####################
    import numpy as np


    # #################### SCALARS ####################
    if isinstance(value, bytes):
        return value.decode()
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.integer):
        return int(value)
    if not isinstance(value, np.ndarray):
        return value


    # #################### STRUCTURES ####################
    if value.dtype.names:
        n = value.size
        columns, dtype = {}, []
        for name in value.dtype.names:
            col = value[name]
            key = name.lower()
            if col.dtype == object:
                first = col[0] if n else None
                if isinstance(first, bytes):
                    data = np.array([c.decode() for c in col])
                    dtype.append((key, f'U{max(1, max(len(s) for s in data))}'))
                elif isinstance(first, np.ndarray):
                    data = np.stack([np.asarray(c, dtype=float) for c in col])
                    dtype.append((key, 'f8', data.shape[1:]))
                else:
                    data = np.asarray(col.tolist())
                    dtype.append((key, data.dtype))
            else:
                data = col.astype(float) if col.dtype.kind == 'f' else col
                dtype.append((key, data.dtype, data.shape[1:]) if data.ndim > 1 else (key, data.dtype))
            columns[key] = data
        out = np.zeros(n, dtype=dtype)
        for key, data in columns.items():
            out[key] = data
        return out


    # #################### ARRAYS ####################
    if value.dtype.kind == 'f':
        return value.astype(float)
    if value.dtype == object and value.size and isinstance(value.flat[0], bytes):
        return np.array([v.decode() for v in value.flat]).reshape(value.shape)


    # #################### OUTPUT ####################
    return value




# ==============================
def _struct_to_dict(value) -> dict:
    """
    Convert a 1-element structure (structured array or dict) into a dict of
    Python scalars (e.g. 'info').
    """

    # #################### MODULES & PACKAGES ####################
    import numpy as np


    # #################### CONVERSION ####################
    if isinstance(value, dict):
        return value
    out = {}
    for name in value.dtype.names:
        v = value[name][0] if value.ndim else value[name]
        out[name.lower()] = v.item() if isinstance(v, np.generic) else v


    # #################### OUTPUT ####################
    return out




# ======================================================================
# ============================== SAVE FILES ==============================
# ======================================================================

# ==============================
def _found(base : str | Path) -> Path:
    """
    Saved file to read: '<base>.idl' (Python), or else '<base>.xdr' (IDL
    version or old runs). The base name is given without extension (the
    names of IACOB-GBAT contain dots, e.g. 'HD36512_..._vr.dat_v13_m33').

    -------
    Returns
    -------
        file : Path
            The '.idl' file if it exists or if the '.xdr' does not exist,
            otherwise the '.xdr' file.
    """

    # #################### MODULES & PACKAGES ####################
    from pathlib import Path


    # #################### PANORAMA ####################
    from gbat.variables import IDL_EXT, SAVE_EXT


    # #################### FILE ####################
    base = str(base)
    idl, xdr = Path(base + SAVE_EXT), Path(base + IDL_EXT)


    # #################### OUTPUT ####################
    return idl if idl.exists() or not xdr.exists() else xdr




# ==============================
def _save_gbat(
               base : str | Path,
               **variables
    ) -> Path:
    """
    Save variables (replacement of IDL SAVE) in '<base>.idl' (pickle), with a
    readable copy '<base>.json' in the same folder if JSON_COPIES is True.

    -------
    Returns
    -------
        file : Path
            The saved '.idl' file.
    """

    # #################### MODULES & PACKAGES ####################
    from pathlib import Path


    # #################### PANORAMA ####################
    from idl_bridge import _save
    from gbat.variables import JSON_COPIES, SAVE_EXT


    # #################### SAVE ####################
    file = Path(str(base) + SAVE_EXT)
    file.parent.mkdir(parents=True, exist_ok=True)
    _save(file, readable=JSON_COPIES, **variables)


    # #################### OUTPUT ####################
    return file




# ==============================
def _restore_gbat(file : str | Path) -> dict:
    """
    Read a save file of IACOB-GBAT (replacement of IDL RESTORE): pickle
    '.idl' written by _save_gbat, or IDL '.xdr' file of the IDL version.

    -------
    Returns
    -------
        data : dict
            Saved variables (lower-case names), with the IDL structures
            converted by _from_idl and 'info' as a dict.
    """

    # #################### MODULES & PACKAGES ####################
    from pathlib import Path


    # #################### PANORAMA ####################
    from idl_bridge import _restore_idl


    # #################### RESTORE ####################
    if not Path(file).exists():
        raise FileNotFoundError(f"Save file not found: '{file}'.")
    data = {k: _from_idl(v) for k, v in _restore_idl(file).items()}
    if 'info' in data:
        data['info'] = _struct_to_dict(data['info'])


    # #################### OUTPUT ####################
    return data




# ======================================================================================
# ============================== IDL TEXT -> MATHTEXT ==============================
# ======================================================================================

# ==============================
_GREEK_LOWER = dict(zip('abcdefghijklmnopqrstuvwx',
                        ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta', 'theta',
                         'iota', 'kappa', 'lambda', 'mu', 'nu', 'xi', 'o', 'pi', 'rho', 'sigma',
                         'tau', 'upsilon', 'phi', 'chi', 'psi', 'omega']))
_GREEK_UPPER = dict(zip('ABCDEFGHIJKLMNOPQRSTUVWX',
                        ['A', 'B', r'\Gamma', r'\Delta', 'E', 'Z', 'H', r'\Theta', 'I', 'K',
                         r'\Lambda', 'M', 'N', r'\Xi', 'O', r'\Pi', 'P', r'\Sigma', 'T',
                         r'\Upsilon', r'\Phi', 'X', r'\Psi', r'\Omega']))




# ==============================
def _idl2mpl(text : str) -> str:
    """
    Convert IDL (Hershey/PostScript) text-formatting commands into Matplotlib
    mathtext.

    ----------
    Parameters
    ----------
        text : str
            IDL text, e.g. 'H!7a!3', '!7n!3!dt!n (kms!u-1!n)' or
            'log !12L!X/!12L!X!9!dn!n!3'.

    -------
    Returns
    -------
        mpl : str
            Mathtext string, e.g. 'H$\\alpha$', '$\\xi_{t}$ (kms$^{-1}$)' or
            'log $\\mathcal{L}$/$\\mathcal{L}_{\\odot}$'.

    -----
    Notes
    -----
        - Fonts: !7 (Greek), !9 (symbols: n -> Sun), !12 (script: L), and !3,
          !X (back to the normal font).
        - Positions: !d/!i (subscript), !u/!e (superscript), !n (normal).
        - '!!' is a literal '!'.
    """

    # #################### PARSER ####################
    segments = [] # (kind, content): 'text' or 'math'
    group = None  # Current sub/superscript: [prefix, list of math tokens]
    font = 'normal'

    def math_token(ch):
        if font == 'greek':
            if ch in _GREEK_LOWER:
                name = _GREEK_LOWER[ch]
                return 'o' if name == 'o' else '\\' + name
            if ch in _GREEK_UPPER:
                return _GREEK_UPPER[ch]
        if font == 'symbol' and ch == 'n':
            return r'\odot'
        if font == 'script' and ch.isalpha():
            return r'\mathcal{' + ch + '}'
        return None

    def add(ch):
        tok = math_token(ch)
        if group is not None:
            if tok is None:
                tok = r'\ ' if ch == ' ' else (r'\mathrm{' + ch + '}' if ch.isalnum() else ch)
            group[1].append(tok)
        elif tok is not None:
            segments.append(('math', tok))
        else:
            segments.append(('text', '\\$' if ch == '$' else ch))

    def close():
        nonlocal group
        if group is not None and group[1]:
            segments.append(('math', group[0] + '{' + ''.join(group[1]) + '}'))
        group = None

    i = 0
    while i < len(text):
        ch = text[i]
        if ch == '!' and i + 1 < len(text):
            cmd = text[i + 1]
            if cmd == '!':
                add('!')
                i += 2
                continue
            if text[i + 1:i + 3] == '12':
                font = 'script'
                i += 3
                continue
            code = cmd.upper()
            if code == '7':
                font = 'greek'
            elif code == '9':
                font = 'symbol'
            elif code in ('3', 'X'):
                font = 'normal'
            elif code in ('D', 'I'):
                close()
                group = ['_', []]
            elif code in ('U', 'E'):
                close()
                group = ['^', []]
            elif code == 'N':
                close()
            i += 2
            continue
        add(ch)
        i += 1
    close()


    # #################### MERGE ####################
    out, math = [], []
    for kind, content in segments:
        if kind == 'math':
            math.append(content)
        else:
            if math:
                out.append('$' + ''.join(math) + '$')
                math = []
            out.append(content)
    if math:
        out.append('$' + ''.join(math) + '$')


    # #################### OUTPUT ####################
    return ''.join(out)
