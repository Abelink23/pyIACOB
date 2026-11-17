from __future__ import annotations

# ========================================================================
# ===================== SAVE FILES OR RESTORE SOLUTIONS ==================
# ========================================================================

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

    return value


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



def _restore_idl(
                file : str | Path) -> dict:
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
    from pathlib import Path


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
        data['vmin'] = 0.0

    return data


def _restore_json(
                json_file : str
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

