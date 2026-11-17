from __future__ import annotations

"""
iacob_gbat_remote — Access to a FASTWIND grid stored in a Nextcloud public share
(e.g. https://cloud.iac.es/index.php/s/<token>).

When CASA (gbat.variables) is the link of a public share, every file of
the grid (list of models, model files) is downloaded the first time it is
needed and kept in a local cache (GRID_CACHE/<token>/...). Only the models of
the subgrid of each analysis are downloaded.

Two layouts of the share are accepted:
    >> the share holds the folder of the grid:   <share>/iac_gbatgrid_HHe_Z100/...
    >> the share IS the folder of the grid:      <share>/iac_gbatgrid_HHe_Z100_allHe.lis ...
The layout is detected with the first file found and then reused.

Files are read, in this order, with the public WebDAV of Nextcloud >= 29
(<host>/public.php/dav/files/<token>/<path>), the download link of the share
(<host>/index.php/s/<token>/download?path=<folder>&files=<file>), or the old
public WebDAV (<host>/public.php/webdav/<path>, user = token, password =
GRID_PASSWORD). The first one that works is used for the next files.
"""

# ======================================================================
# ============================== URLS ==============================
# ======================================================================

# ==============================
_MISSING: set = set()     # Remote files known not to exist (no repeated requests)
_LAYOUT: dict = {}        # token -> number of leading path components to drop (0 or 1)




# ==============================
def _is_url(path) -> bool:
    """
    True if path is an http(s) link.
    """

    # #################### OUTPUT ####################
    return str(path).startswith(('http://', 'https://'))




# ==============================
def _split_share(url : str) -> tuple[str, str, str]:
    """
    Split the link of a file inside a public share.

    ----------
    Parameters
    ----------
        url : str
            e.g. 'https://cloud.iac.es/index.php/s/TOKEN//iac_gbatgrid_HHe_Z100/x.lis'.

    -------
    Returns
    -------
        host : str
            'https://cloud.iac.es'.

        token : str
            Token of the share.

        rel : str
            Path of the file inside the share, without repeated or leading
            slashes ('iac_gbatgrid_HHe_Z100/x.lis').
    """

    # #################### MODULES & PACKAGES ####################
    import re


    # #################### PARSING ####################
    m = re.match(r'^(https?://[^/]+)(?:/index\.php)?/s/([^/?#]+)/*(.*)$', str(url))
    if not m:
        raise ValueError(f"'{url}' is not the link of a Nextcloud public share "
                         "(expected https://<host>/index.php/s/<token>).")
    host, token, rel = m.groups()
    rel = re.sub(r'/+', '/', rel).strip('/')


    # #################### OUTPUT ####################
    return host, token, rel




# =============================================================================
# ============================== DOWNLOAD ==============================
# =============================================================================

# ==============================
def _requests(
              host  : str,
              token : str,
              rel   : str
    ) -> list[tuple[str, dict]]:
    """
    URLs (and headers) tried to download a file of the share, in order.

    -----
    Notes
    -----
        1. Public WebDAV of Nextcloud >= 29: <host>/public.php/dav/files/<token>/<path>
           (no authentication for public shares).
        2. Download link of the share: <host>/index.php/s/<token>/download?path=<folder>&files=<file>
           (the one used by the web page; no authentication).
        3. Old public WebDAV: <host>/public.php/webdav/<path>, user = token,
           password = GRID_PASSWORD (older Nextcloud/ownCloud, and protected shares).
        - The Authorization header is only sent to the old WebDAV endpoint:
          sending it to the others makes some servers answer 401 even for
          public shares.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import base64
    import urllib.parse


    # ########## PANORAMA ##########
    from gbat.variables import GRID_PASSWORD


    # #################### URLS ####################
    folder, _, name = rel.rpartition('/')
    quoted = urllib.parse.quote(rel)
    auth = base64.b64encode(f'{token}:{GRID_PASSWORD}'.encode()).decode()
    common = {'User-Agent': 'IACOB-GBAT', 'X-Requested-With': 'XMLHttpRequest'}


    # #################### OUTPUT ####################
    return [
            (f'{host}/public.php/dav/files/{token}/{quoted}', dict(common)),
            (f'{host}/index.php/s/{token}/download?' + urllib.parse.urlencode({'path': '/' + folder, 'files': name}),
             {'User-Agent': 'IACOB-GBAT'}),
            (f'{host}/public.php/webdav/{quoted}', dict(common, Authorization=f'Basic {auth}')),
           ]




# ==============================
_WORKING: dict = {} # token -> index of the URL scheme that worked (tried first afterwards)




# ==============================
def _download(
              host  : str,
              token : str,
              rel   : str,
              out   : Path
    ) -> bool:
    """
    Download one file of the share into out.

    -------
    Returns
    -------
        ok : bool
            True if downloaded, False if the file does not exist in the share.

    -----
    Notes
    -----
        - The three URL schemes of _requests are tried in order (the one that
          worked before first). A 401/403 of one of them is not final: the
          next one is tried.
        - Nextcloud answers some errors with an HTML page: HTML answers are
          treated as failures.
        - Raises ConnectionError if the file could not be read for any reason
          other than not existing (network, password, downloads disabled ...).
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import os
    import urllib.error
    import urllib.request


    # ########## PANORAMA ##########
    from gbat.variables import GRID_TIMEOUT


    # #################### ORDER ####################
    requests = _requests(host, token, rel)
    order = list(range(len(requests)))
    if token in _WORKING:
        order.remove(_WORKING[token])
        order.insert(0, _WORKING[token])


    # #################### DOWNLOAD ####################
    errors = []
    for k in order:
        url, headers = requests[k]
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=GRID_TIMEOUT) as resp:
                if 'text/html' in resp.headers.get('Content-Type', ''):
                    errors.append((url, 'HTML'))
                    continue
                out.parent.mkdir(parents=True, exist_ok=True)
                tmp = out.with_name(out.name + '.part')
                size = int(resp.headers.get('Content-Length') or 0)
                done, step = 0, 0.1
                with open(tmp, 'wb') as fh:
                    while chunk := resp.read(1 << 20):
                        fh.write(chunk)
                        done += len(chunk)
                        if size > 20e6 and done / size >= step: # Progress of big files (sets of models)
                            print(f'       {out.name}: {done / 1e6:7.1f} / {size / 1e6:.1f} MB')
                            step += 0.1
                os.replace(tmp, out)
                _WORKING[token] = k
                return True
        except urllib.error.HTTPError as e:
            errors.append((url, e.code))
        except urllib.error.URLError as e:
            raise ConnectionError(f"Could not connect to '{host}': {e.reason}") from e


    # #################### NOT FOUND / ERRORS ####################
    codes = [c for _, c in errors]
    if any(c == 404 for c in codes) and all(c in (404, 'HTML', 401, 403, 405) for c in codes):
        return False # At least one scheme works for this share and says the file does not exist
    details = '\n  '.join(f'{u} -> {c}' for u, c in errors)
    hint = ''
    if all(c in (401, 403) for c in codes):
        hint = ("\nThe share refused all the requests. Check that the link is correct, that downloads are "
                "allowed in the share ('Hide download' disabled), and set GRID_PASSWORD only if the share "
                "has a password.")
    raise ConnectionError(f'Could not download {rel} from the grid share:\n  {details}{hint}')




# ==============================
def _resolve(path : str | Path) -> Path | None:
    """
    Local file of a file of the grid.

    ----------
    Parameters
    ----------
        path : str | Path
            Local path, or link of a file inside a public share.

    -------
    Returns
    -------
        local : Path | None
            The local file (downloaded to the cache if path is a link), or
            None if it does not exist.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    from pathlib import Path


    # ########## PANORAMA ##########
    from gbat.variables import GRID_CACHE


    # #################### LOCAL FILE ####################
    if not _is_url(path):
        p = Path(path)
        return p if p.is_file() else None


    # #################### REMOTE FILE ####################
    host, token, rel = _split_share(path)
    cache = Path(GRID_CACHE).expanduser() / token
    parts = rel.split('/')


    # ########## MODEL FILE: SET OF MODELS (.tgz) FIRST ##########
    if rel.endswith(('.xdr', '.idl')) and len(parts) > 1:
        folder = str(path).rsplit('/', 1)[0]
        local_set = _set_folder(folder)
        if local_set is None and _fetch_set(folder):
            local_set = _set_folder(folder)
        if local_set is not None: # The set is unpacked: the model is there or it does not exist
            local = local_set / parts[-1]
            return local if local.is_file() else None


    # ########## SINGLE FILE ##########
    drops = [_LAYOUT[token]] if token in _LAYOUT else ([0, 1] if len(parts) > 1 else [0])
    for drop in drops:
        cand = '/'.join(parts[drop:])
        local = cache / cand
        if local.is_file():
            _LAYOUT.setdefault(token, drop)
            return local
        if (token, cand) in _MISSING:
            continue
        if _download(host, token, cand, local):
            _LAYOUT.setdefault(token, drop)
            return local
        _MISSING.add((token, cand))


    # #################### OUTPUT ####################
    return None




# ==============================
_SETS_TRIED: set = set() # Sets of models whose .tgz was already looked for in this session




# ==============================
def _set_folder(folder : str) -> Path | None:
    """
    Local folder of an unpacked set of models (link of the folder of the set
    in the share), or None if the set has not been downloaded.

    -----
    Notes
    -----
        - A set is complete when its folder holds the marker file
          '.iacob_gbat_complete', written after unpacking its .tgz.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    from pathlib import Path


    # ########## PANORAMA ##########
    from gbat.variables import GRID_CACHE


    # #################### FOLDER ####################
    host, token, rel = _split_share(folder)
    parts = rel.split('/')
    drops = [_LAYOUT[token]] if token in _LAYOUT else [0, 1]
    for drop in drops:
        local = Path(GRID_CACHE).expanduser() / token / '/'.join(parts[drop:])
        if (local / '.iacob_gbat_complete').is_file():
            return local


    # #################### OUTPUT ####################
    return None




# ==============================
def _fetch_set(folder : str) -> bool:
    """
    Download and unpack the .tgz file of a set of models.

    ----------
    Parameters
    ----------
        folder : str
            Link of the folder of the set in the share, e.g.
            '<share>/iac_gbatgrid_HHe_Z100/iac_gbatgrid_HHe_Z100_He08_xdr'.
            The set is looked for as '<folder><ext>', ext in GRID_TGZ_EXT
            (e.g. iac_gbatgrid_HHe_Z100_He08_xdr.tgz).

    -------
    Returns
    -------
        ok : bool
            True if the set was downloaded and unpacked, False if there is no
            .tgz for this set in the share.

    -----
    Notes
    -----
        - The model files (.xdr, .idl) of the archive are written in the
          folder of the set, whatever the folders inside the archive (only
          their names are used, so the archive cannot write outside the
          cache). Files inside another '*_xdr' folder of the archive go to
          that folder.
        - The .tgz is deleted after unpacking unless GRID_KEEP_TGZ is True.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import shutil
    import tarfile
    from pathlib import Path


    # ########## PANORAMA ##########
    from gbat.variables import GRID_KEEP_TGZ, GRID_TGZ_EXT


    # #################### ALREADY TRIED ####################
    if folder in _SETS_TRIED:
        return False
    _SETS_TRIED.add(folder)


    # #################### DOWNLOAD ####################
    name = folder.rstrip('/').rsplit('/', 1)[-1]
    for ext in GRID_TGZ_EXT:
        host, token, rel = _split_share(folder + ext)
        print(f'   - Looking for the set of models {name}{ext} in the grid share ...')
        archive = _resolve(folder + ext)
        if archive is not None:
            break
    else:
        return False


    # #################### UNPACK ####################
    dest = archive.parent / name
    dest.mkdir(parents=True, exist_ok=True)
    print(f'   - Unpacking {archive.name} in {dest} ...')
    n = 0
    with tarfile.open(archive, 'r:*') as tar:
        for member in tar:
            if not member.isfile() or not member.name.endswith(('.xdr', '.idl', '.lis')):
                continue
            pieces = Path(member.name).parts
            parent = pieces[-2] if len(pieces) > 1 else ''
            target = archive.parent / parent if parent.endswith('_xdr') and parent != name else dest
            target.mkdir(parents=True, exist_ok=True)
            with tar.extractfile(member) as src, open(target / pieces[-1], 'wb') as out:
                shutil.copyfileobj(src, out)
            n += 1
    (dest / '.iacob_gbat_complete').write_text(f'{n} files from {archive.name}\n')
    print(f'   - {n} model files unpacked')
    if not GRID_KEEP_TGZ:
        archive.unlink()


    # #################### OUTPUT ####################
    return True




# ==============================
def _prefetch(
              bases   : list[str],
              exts    : tuple[str, ...] = ('.xdr', '.idl')
    ) -> None:
    """
    Download in parallel (GRID_THREADS) the model files of a subgrid, so that
    the chi2 loop does not wait for each download.

    ----------
    Parameters
    ----------
        bases : list[str]
            Links of the models without extension.

        exts : tuple[str, ...], optional
            Extensions tried for each model.
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    from concurrent.futures import ThreadPoolExecutor
    from pathlib import Path


    # ########## PANORAMA ##########
    from gbat.variables import GRID_CACHE, GRID_THREADS


    # #################### MODELS TO DOWNLOAD ####################
    remote = [b for b in bases if _is_url(b)]
    if not remote:
        return


    # ########## SETS OF MODELS (.tgz) ##########
    folders = list(dict.fromkeys(b.rsplit('/', 1)[0] for b in remote))
    missing = [f for f in folders if _set_folder(f) is None]
    if missing:
        print(f'   - The subgrid needs {len(folders)} set(s) of models; {len(missing)} not downloaded yet')
    for folder in missing:
        _fetch_set(folder)
    remote = [b for b in remote if _set_folder(b.rsplit('/', 1)[0]) is None] # Models without a .tgz set
    if not remote:
        return


    # ########## SINGLE MODEL FILES ##########
    _resolve(remote[0] + exts[0]) or _resolve(remote[0] + exts[-1]) # Detects the layout of the share
    host, token, _ = _split_share(remote[0])
    drop = _LAYOUT.get(token, 0)
    cache = Path(GRID_CACHE).expanduser() / token
    todo = [b for b in remote
            if not any((cache / '/'.join(_split_share(b + e)[2].split('/')[drop:])).is_file() for e in exts)]
    if not todo:
        return
    print(f'   - Downloading {len(todo)} models of the subgrid from the grid share (cache: {cache}) ...')


    # #################### DOWNLOAD ####################
    def one(base):
        for e in exts:
            if _resolve(base + e) is not None:
                return True
        return False

    with ThreadPoolExecutor(max_workers=GRID_THREADS) as pool:
        found = sum(pool.map(one, todo))
    print(f'   - {found} of {len(todo)} models downloaded')




# ===========================================================================
# ============================== LISTING ==============================
# ===========================================================================

# ==============================
def list_grid(
              url  : str | None = None,
              path : str = ''
    ) -> list[str]:
    """
    List a folder of the grid share (WebDAV PROPFIND), to check its layout.

    ----------
    Parameters
    ----------
        url : str | None, optional
            Link of the share. Default: CASA (gbat.variables).

        path : str, optional
            Folder inside the share. Default '' (root).

    -------
    Returns
    -------
        names : list[str]
            Files and folders (folders end with '/').

    --------
    Examples
    --------
        >>> from gbat.remote import list_grid
        >>> list_grid()                          # root of the share
        >>> list_grid(path='iac_gbatgrid_HHe_Z100')
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import base64
    import urllib.error
    import urllib.parse
    import urllib.request
    import xml.etree.ElementTree as ET


    # ########## PANORAMA ##########
    from gbat.variables import CASA, GRID_PASSWORD, GRID_TIMEOUT


    # #################### REQUEST ####################
    host, token, rel = _split_share((url or CASA).rstrip('/') + '/' + path)
    auth = base64.b64encode(f'{token}:{GRID_PASSWORD}'.encode()).decode()
    quoted = urllib.parse.quote(rel)
    xhr = {'X-Requested-With': 'XMLHttpRequest'}
    basic = {'Authorization': f'Basic {auth}'}
    attempts = [(f'{host}/public.php/dav/files/{token}/{quoted}', xhr),
                (f'{host}/public.php/dav/files/{token}/{quoted}', dict(xhr, **basic)),
                (f'{host}/public.php/dav/files/{token}/{quoted}', {}),
                (f'{host}/public.php/webdav/{quoted}', dict(xhr, **basic)),
                (f'{host}/public.php/webdav/{quoted}', basic)]
    tree, errors = None, []
    for u, h in attempts:
        req = urllib.request.Request(u, method='PROPFIND', headers=dict(h, Depth='1', **{'User-Agent': 'IACOB-GBAT'}))
        try:
            with urllib.request.urlopen(req, timeout=GRID_TIMEOUT) as resp:
                tree = ET.fromstring(resp.read())
                break
        except urllib.error.HTTPError as e:
            errors.append(f'{u} [{", ".join(h) or "no extra headers"}] -> {e.code}')
    if tree is None:
        raise ConnectionError('Could not list the grid share (its WebDAV interface is not available; the analysis '
                              'can still download the files with the download link of the share):\n  ' + '\n  '.join(errors))


    # #################### NAMES ####################
    names = []
    for href in tree.iter('{DAV:}href'):
        text = urllib.parse.unquote(href.text)
        if text.rstrip('/').endswith(('/webdav', '/' + token) if not rel else '/' + rel.rstrip('/')):
            continue # The listed folder itself
        names.append(text.rstrip('/').rsplit('/', 1)[-1] + ('/' if text.endswith('/') else ''))
    for n in names:
        print(n)


    # #################### OUTPUT ####################
    return names




# ==============================
def check_grid(
               metal : str | None = None,
               url   : str | None = None
    ) -> str | None:
    """
    Diagnose the access to the grid share: try to download the list of models
    of the grid with every URL scheme (see _requests) and print the answer of
    the server to each one.

    ----------
    Parameters
    ----------
        metal : str | None, optional
            Metallicity label of the grid (default: METAL + METAL_EXTRA).

        url : str | None, optional
            Link of the share (default: CASA).

    -------
    Returns
    -------
        working : str | None
            The first URL that returned the file, or None.

    -----
    Notes
    -----
        - Both layouts of the share are tested: the list of models inside the
          folder iac_gbatgrid_HHe_<metal>/, or at the root of the share.
        - Nothing is written to the cache.

    --------
    Examples
    --------
        >>> from gbat.remote import check_grid
        >>> check_grid()
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import urllib.error
    import urllib.request


    # ########## PANORAMA ##########
    from gbat.variables import CASA, GRID_TIMEOUT, METAL, METAL_EXTRA


    # #################### FILES TO TEST ####################
    grid = 'iac_gbatgrid_HHe_' + (metal or METAL + METAL_EXTRA)
    lis = grid + '_allHe.lis'
    host, token, _ = _split_share(url or CASA)
    print(f'Share: {host}  token: {token}')


    # #################### REQUESTS ####################
    for rel in (grid + '/' + lis, lis):
        print(f'\n--- {rel}')
        for u, h in _requests(host, token, rel):
            try:
                with urllib.request.urlopen(urllib.request.Request(u, headers=h), timeout=GRID_TIMEOUT) as resp:
                    ctype = resp.headers.get('Content-Type', '')
                    head = resp.read(300)
                    ok = 'text/html' not in ctype
                    print(f'  {"OK  " if ok else "HTML"} {resp.status} {ctype:30s} {u}')
                    print('       first bytes: ' + head[:120].decode(errors='replace').replace('\n', ' | '))
                    if ok:
                        print(f'\n==> Access works with: {u}')
                        return u
            except urllib.error.HTTPError as e:
                print(f'  FAIL {e.code} {u}')
            except urllib.error.URLError as e:
                print(f'  FAIL {e.reason} {u}')


    # #################### OUTPUT ####################
    print('\n==> None of the URLs returned the list of models. Please send this output.')
    return None




# ==============================
def download_grid(
                  sets  : str | int | list | None = None,
                  teff  : tuple[float, float] | None = None,
                  grav  : tuple[float, float] | None = None,
                  yhe   : tuple[float, float] | None = None,
                  q     : tuple[float, float] | None = None,
                  beta  : tuple[float, float] | None = None
    ) -> list[Path]:
    """
    Choose and download sets of models (.tgz files) of the grid share.

    The sets are the sub-folders of the list of models (.lis), e.g.
    iac_gbatgrid_HHe_Z100_He08_xdr, each one stored in the share as a .tgz
    with its .xdr files. They are downloaded and unpacked in GRID_CACHE, where
    IACOB-GBAT reads them. Sets already downloaded are not downloaded again.

    ----------
    Parameters
    ----------
        sets : str | int | list | None, optional
            Sets to download: 'all', the number of a set in the table (or a
            list of them), or set names. If None and no range is given, the
            table of sets is shown and the choice is asked in the terminal.

        teff, grav, yhe, q, beta : tuple[float, float] | None, optional
            Download the sets with at least one model in these ranges
            (Teff [K], log g, Y(He), -log Q, beta), as in chi_global.

    -------
    Returns
    -------
        folders : list[Path]
            Local folders of the selected sets.

    -----
    Notes
    -----
        - The grid is the one of CASA, METAL and METAL_EXTRA
          (gbat.variables).
        - chi_global downloads automatically the sets needed by each subgrid,
          so this function is only needed to download sets in advance (e.g.
          before working offline).

    --------
    Examples
    --------
        >>> from gbat.remote import download_grid
        >>> download_grid()                    # table of sets + choice in the terminal
        >>> download_grid([1, 2])              # sets 1 and 2 of the table
        >>> download_grid('all')
        >>> download_grid(yhe=(0.07, 0.12))    # sets with Y(He) in [0.07, 0.12]
    """

    # #################### MODULES & PACKAGES ####################
    # ########## PYTHON ##########
    import numpy as np


    # ########## PANORAMA ##########
    from gbat.grid import _read_models
    from gbat.variables import CASA, METAL, METAL_EXTRA


    # #################### GRID ####################
    if not _is_url(CASA):
        print(f'CASA is a local folder ({CASA}): nothing to download.')
        return []
    casa_red = CASA + '/iac_gbatgrid_HHe_' + METAL + METAL_EXTRA + '/'
    file_models = casa_red + '/iac_gbatgrid_HHe_' + METAL + METAL_EXTRA + '_allHe.lis'
    full = (-1e30, 1e30)
    allm = _read_models(file_models, full, full, full, full, full)
    names = list(dict.fromkeys(allm['red']))


    # #################### TABLE OF SETS ####################
    rng = lambda a, f: f'{a.min() * f:g}-{a.max() * f:g}' if a.min() != a.max() else f'{a.min() * f:g}'
    print('')
    print(f'  #  {"set of models":40s} {"models":>6s}  {"Teff(kK)":>9s}  {"logg":>9s}  {"Y(He)":>9s}  '
          f'{"-logQ":>9s}  {"beta":>7s}  downloaded')
    for k, red in enumerate(names, 1):
        mm = allm['red'] == red
        done = _set_folder(casa_red + red.rstrip('/')) is not None
        print(f'{k:3d}  {red.rstrip("/"):40s} {mm.sum():6d}  {rng(allm["teff"][mm], 1):>9s}  {rng(allm["grav"][mm], 1):>9s}  '
              f'{rng(allm["he"][mm], 1):>9s}  {rng(allm["q"][mm], 1):>9s}  {rng(allm["beta"][mm], 1):>7s}  '
              f'{"yes" if done else "no"}')
    print('')


    # #################### SELECTION ####################
    if any(r is not None for r in (teff, grav, yhe, q, beta)):
        sub = _read_models(file_models, teff or full, grav or full, yhe or full, q or full, beta or full)
        chosen = list(dict.fromkeys(sub['red']))
    else:
        if sets is None:
            sets = input("Sets to download (numbers separated by commas, 'all', or empty to cancel): ").strip()
            if not sets:
                return []
            sets = 'all' if sets.lower() == 'all' else [int(v) for v in sets.replace(',', ' ').split()]
        if isinstance(sets, (str, int)):
            sets = names if sets == 'all' else [sets]
        chosen = []
        for v in sets:
            if isinstance(v, (int, np.integer)):
                if not 1 <= v <= len(names):
                    raise ValueError(f'Set number {v} out of range (1-{len(names)}).')
                chosen.append(names[v - 1])
            else:
                match = [r for r in names if r.rstrip('/') == str(v).rstrip('/')]
                if not match:
                    raise ValueError(f"Set '{v}' not found. Available: {', '.join(r.rstrip('/') for r in names)}.")
                chosen.append(match[0])


    # #################### DOWNLOAD ####################
    folders = []
    for red in chosen:
        folder = casa_red + red.rstrip('/')
        if _set_folder(folder) is None and not _fetch_set(folder):
            print(f'WARNING ... no .tgz found for the set {red.rstrip("/")} in the grid share')
            continue
        folders.append(_set_folder(folder))
        print(f'   - {red.rstrip("/")}: {folders[-1]}')


    # #################### OUTPUT ####################
    return folders
