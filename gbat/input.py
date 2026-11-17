from __future__ import annotations
from dataclasses import dataclass, field

# =========================================================================
# ============================== CONFIGURATION ==============================
# =========================================================================

# ==============================
@dataclass
class Config:
    """
    Main input information of the computation (IDL COMMON blocks share1 and
    share2, filled by chi_input).

    -----
    Notes
    -----
        - resol and res0 can be modified line by line by _extra_info_obs, and
          the modified values are kept for the next lines (as the IDL COMMON
          variables).
    """
    root: Path                      # Working folder (IDL: current directory)
    casa_red: str                   # Folder of the FASTWIND grid
    metal: str                      # Metallicity of the grid
    metal_extra: str                # Extra label of the grid folder
    resol: float                    # Resolving power of the observations
    res0: float                     # Step [A] of the convolved FASTWIND profiles
    extra_info: int                 # Use EXTRA_INFO_OBS (0/1)
    file_models: str                # List of models of the grid
    lines: list[str] = field(default_factory=list) # FASTWIND labels of the lines with chi_in > 0
    ps_viewer: str = 'evince'       # PDF viewer




# ==============================
def _chi_input(
               obs_run : str,
               inst    : str | None,
               root    : str | Path = '.'
    ) -> Config:
    """
    Main input information of the computation (translation of
    gbat_chi_input.pro).

    ----------
    Parameters
    ----------
        obs_run : str
            Observing run (sub-folder of spectra_orig, spectra_cut ...). The
            lines of the analysis are read from 'gbat_linesHHe_<obs_run>.dat'.

        inst : str | None
            Instrument, one of the keys of INSTRUMENTS in
            gbat.variables.

        root : str | Path, optional
            Working folder. Defaults to '.'.

    -------
    Returns
    -------
        cfg : Config
            Configuration of the run.

    -----
    Notes
    -----
        - The user-editable information (PDF viewer, location of the grids,
          metallicity, instruments and lines) is in gbat.variables.
        - As in IDL, the line file is copied to 'gbat_linesHHe.dat', which
          is the file read by _read_lines_hhe.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import shutil
    from pathlib import Path


    # ########## PANORAMA ##########
    from gbat.variables import (CASA, INSTRUMENTS, LINES_FILE, LINES_FILE_ACTIVE,
                                      METAL, METAL_EXTRA, PS_VIEWER)


    # #################### INSTRUMENT ####################
    if not inst:
        raise ValueError("Please select one inst option ('inst'). "
                         f"Available options: {', '.join(INSTRUMENTS)}.")
    if inst not in INSTRUMENTS:
        raise ValueError(f"Instrument 'inst' = '{inst}' is not defined in INSTRUMENTS (gbat.variables). "
                         f"Available options: {', '.join(INSTRUMENTS)}.")
    resol, res0, extra_info = INSTRUMENTS[inst]


    # #################### GRID ####################
    casa_red = CASA + '/iac_gbatgrid_HHe_' + METAL + METAL_EXTRA + '/'
    file_models = casa_red + 'iac_gbatgrid_HHe_' + METAL + METAL_EXTRA + '_allHe.lis'


    # #################### LINES ####################
    root = Path(root)
    file_lines = root / LINES_FILE.format(obs_run=obs_run)
    if not file_lines.exists():
        raise FileNotFoundError(f"Line file '{file_lines}' not found. Create it by modifying "
                                "gbat_linesHHe_benchmark.dat.")
    active = root / LINES_FILE_ACTIVE
    if file_lines.resolve() != active.resolve():
        shutil.copy(file_lines, active) # IDL: spawn, 'cp '+file_lines+' gbat_linesHHe.dat'
    dat = _read_lines_hhe(root)
    lines = [lab for lab, chi_in in zip(dat['fwlab'], dat['chi_in']) if chi_in > 0]


    # #################### OUTPUT ####################
    return Config(root=root, casa_red=casa_red, metal=METAL, metal_extra=METAL_EXTRA,
                  resol=float(resol), res0=float(res0), extra_info=int(extra_info),
                  file_models=file_models, lines=lines, ps_viewer=PS_VIEWER)




# ==============================
def _read_lines_hhe(root : str | Path = '.') -> dict[str, np.ndarray]:
    """
    Read the lines considered in the analysis from 'gbat_linesHHe.dat'
    (translation of gbat_read_linesHHe.pro).

    -------
    Returns
    -------
        dat : dict[str, np.ndarray]
            lab (label for the plots, with IDL text commands), lam (central
            wavelength [A]), fwlab (name of the line in FASTWIND), chi_in (1:
            the line is used in the chi2 computation), step2_in (1: used in
            the final results) and rc_pos (row x 10 + column of the line in
            the profiles figure).
    """

    # #################### MODULES & PACKAGES ####################
    from pathlib import Path


    # #################### PANORAMA ####################
    from gbat.io import _readcol
    from gbat.variables import LINES_FILE_ACTIVE


    # #################### READING ####################
    lab, lam, fwlab, chi_in, step2_in, rc_pos = _readcol(Path(root) / LINES_FILE_ACTIVE, fmt='(A,F,A,I,I,I)')


    # #################### OUTPUT ####################
    return dict(lab=lab, lam=lam, fwlab=fwlab, chi_in=chi_in, step2_in=step2_in, rc_pos=rc_pos)




# ==============================
def _extra_info_obs(
                    cfg : Config,
                    x0  : float
    ) -> None:
    """
    Wavelength-dependent resolving power and step of the convolved profiles
    (translation of gbat_extra_info_obs.pro). Only applied if
    cfg.extra_info is 1; cfg.resol and cfg.res0 are modified in place.
    """

    # #################### PANORAMA ####################
    from gbat.variables import EXTRA_INFO_OBS


    # #################### RANGES ####################
    if cfg.extra_info:
        for lmin, lmax, resol, res0 in EXTRA_INFO_OBS:
            if lmin <= x0 <= lmax:
                cfg.resol, cfg.res0 = float(resol), float(res0)




# ==============================
def _lineshhe(
              line : str,
              root : str | Path = '.'
    ) -> dict:
    """
    Label and wavelength limits of a line for the profiles figure
    (translation of gbat_lineshhe.pro).

    -------
    Returns
    -------
        labels : dict
            li (label), la and lb (limits [A], lam -/+ X1_PROF).
    """

    # #################### MODULES & PACKAGES ####################
    import numpy as np


    # #################### PANORAMA ####################
    from gbat.variables import X1_PROF


    # #################### LABELS ####################
    dat = _read_lines_hhe(root)
    i = np.nonzero(dat['fwlab'] == line)[0][0]


    # #################### OUTPUT ####################
    return dict(li=str(dat['lab'][i]), la=float(dat['lam'][i]) - X1_PROF, lb=float(dat['lam'][i]) + X1_PROF)




# =====================================================================
# ============================== QUESTIONS ==============================
# =====================================================================

# ==============================
def _question(lab : str) -> None:
    """
    Print a message in a box (translation of gbat_question.pro).
    """

    # #################### PRINT ####################
    print(' ')
    print('----------------------------------------------------------------------')
    print(lab)
    print('----------------------------------------------------------------------')
    print(' ')




# ==============================
def _ask(
         lab   : str,
         valid : tuple[int, ...] = (0, 1)
    ) -> int:
    """
    Ask a question (gbat_question + IDL READ) until a valid integer answer is
    given.
    """

    # #################### QUESTION ####################
    _question(lab)
    while True:
        answer = input(': ').strip()
        try:
            value = int(answer)
        except ValueError:
            value = None
        if value in valid:
            return value
        print(f"Invalid answer '{answer}'. Type one of: {', '.join(str(v) for v in valid)}.")
