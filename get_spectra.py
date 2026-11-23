# Standard library
import shutil
import tarfile
from urllib.parse import parse_qs, unquote, urlparse

# Third-party packages
import pyvo
import requests
from astropy.time import Time
from astroquery.eso import Eso

from db import *

__all__ = ["get_mercator", "get_feros", "convert_vot_to_fits"]
SSA_URL = "https://mercatorvo.ster.kuleuven.be/hermes/q/ssa/ssap.xml"
SPECTRUM_EXTS = {".fits", ".fit", ".fts", ".vot", ".xml", ".gz", ".dat", ".txt"}

_ESO_CLIENT = None


# =============================================================================
# Shared helpers
# =============================================================================

def _is_coord_pair(obj):
    """Return True if obj is a single (ra, dec) numeric pair."""
    return (
        isinstance(obj, (tuple, list))
        and len(obj) == 2
        and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in obj)
    )


def _load_targets_input(targets):
    """
    Normalize input into a list of targets.
    Accepts:
      - Path/filename ending in .lst, .txt, or .fits (using db.findlist / db.findtable)
      - A single star name string (e.g. "V356 Per")
      - A single (ra_deg, dec_deg) tuple (e.g. (34.24046, 57.13020))
      - A list/iterable of star names and/or coordinate tuples
    """
    if isinstance(targets, (str, Path)):
        s = str(targets)
        if s.endswith((".lst", ".txt")):
            return list(findlist(s))
        if s.endswith(".fits"):
            tbl = findtable(s)
            colname = input("Enter the column name with target names: ")
            if colname not in tbl.colnames:
                raise KeyError(
                    f"Column '{colname}' not found in {s}. Available: {tbl.colnames}"
                )
            return list(tbl[colname])
        return [s]

    if _is_coord_pair(targets):
        return [targets]

    return list(targets)


def _get_eso_client():
    """Lazily initialize and authenticate the ESO client on first use."""
    global _ESO_CLIENT
    if _ESO_CLIENT is None:
        eso = Eso()
        eso_user = input("Provide the username of your ESO account: ")
        eso.login(username=eso_user, store_password=True)
        eso.ROW_LIMIT = -1
        _ESO_CLIENT = eso
    return _ESO_CLIENT


# =============================================================================
# Mercator / HERMES helpers & main function
# =============================================================================

def _get_row_value(row, key):
    """Safely retrieve a column value from a pyvo SSARecord (with case-insensitive fallback)."""
    try:
        val = row[key]
    except (KeyError, IndexError):
        val = None
        try:
            for k in row.keys():
                if k.lower() == key.lower():
                    val = row[k]
                    break
        except Exception:
            return None

    if val is None or np.ma.is_masked(val):
        return None
    return val


def _format_dateobs_for_filename(val):
    """Convert an ssa_dateobs value (MJD float or ISO string) to 'YYYYMMDD_HHMMSS'."""
    if val is None or np.ma.is_masked(val):
        return None

    if isinstance(val, (int, float, np.integer, np.floating)):
        if 40000 < float(val) < 80000:
            return Time(float(val), format="mjd").strftime("%Y%m%d_%H%M%S")
        return None

    s = str(val).strip()
    if not s:
        return None

    try:
        return Time(s).strftime("%Y%m%d_%H%M%S")
    except Exception:
        m = re.search(
            r"(\d{4})-?(\d{2})-?(\d{2})[T\s_]?(\d{2}):?(\d{2}):?(\d{2})", s
        )
        if m:
            y, mo, d, h, mi, sec = m.groups()
            return f"{y}{mo}{d}_{h}{mi}{sec}"
        return None


def _build_mercator_filename(row):
    """
    Build a Mercator/HERMES filename in the format:
    '<ssa_targname>_<YYYYMMDD>_<HHMMSS>_raw.fits'.
    Returns None if ssa_targname or ssa_dateobs is missing/unparseable.
    """
    targname = _get_row_value(row, "ssa_targname")
    dateobs = _get_row_value(row, "ssa_dateobs")
    if targname is None or dateobs is None:
        return None

    clean_target = re.sub(r"[^\w\-]", "_", str(targname).strip())
    clean_target = re.sub(r"_+", "", clean_target).strip("_")
    dt_str = _format_dateobs_for_filename(dateobs)
    if not clean_target or not dt_str:
        return None

    return f"{clean_target}_{dt_str}_raw.fits"


def _simplify_hermes_filename(raw_name, fallback_name="hermes_spectrum.fits"):
    """
    Simplify a Mercator/HERMES filename to '<YYYYMMDD>_<obsID>_<mode>.fits'.
    Example:
      'hermes_data_hermes_20201223_reduced_00985062_HRF_OBJ_ext_CosmicsRemoved_log_merged_cf.vot'
      -> '20201223_00985062_HRF.fits'
    """
    text = unquote(str(raw_name or "")).strip()

    # Match both flattened ('20201223_reduced_00985062_HRF') and path ('20201223/reduced/00985062_HRF')
    m = re.search(r"(\d{8})[/\\_]reduced[/\\_](\d+)_([A-Za-z0-9]+)", text)
    if m:
        date_str, obs_id, mode = m.groups()
        return f"{date_str}_{obs_id}_{mode}.fits"

    # Fallback: strip boilerplate prefixes/suffixes and force .fits extension
    stem = Path(text).stem if text else Path(fallback_name).stem
    stem = re.sub(r"^hermes_data_hermes_", "", stem)
    stem = re.sub(r"_OBJ_ext_CosmicsRemoved_log_merged_cf$", "", stem)
    stem = stem.replace("_reduced_", "_")
    return f"{stem}.fits"


def _extract_row_stems(row):
    """Extract candidate file stems (new format, simplified, and raw) from a row."""
    stems = []

    # 0. New target-date-time format (<ssa_targname>_<YYYYMMDD>_<HHMMSS>_raw)
    new_filename = _build_mercator_filename(row)
    if new_filename:
        stems.append(Path(new_filename).stem)

    candidates = []

    url = row.getdataurl()
    if url:
        parsed = urlparse(str(url))
        qs = parse_qs(parsed.query)
        for param in ("ID", "id", "pubdid"):
            if param in qs:
                candidates.extend(qs[param])

        path_str = unquote(parsed.path)
        if Path(path_str).suffix.lower() in SPECTRUM_EXTS:
            candidates.append(path_str)

    pubdid = _get_row_value(row, "ssa_pubdid")
    if pubdid:
        candidates.append(str(pubdid))

    for raw in candidates:
        s = unquote(str(raw)).strip()
        if not s:
            continue
        if "?" in s:
            s = s.split("?", 1)[1]
        s = s.split("&", 1)[0].split("#", 1)[0].strip("/")
        if not s:
            continue

        # 1. Legacy simplified stem (e.g. 20201223_00985062_HRF)
        simple_stem = Path(_simplify_hermes_filename(s)).stem
        if simple_stem and simple_stem not in stems:
            stems.append(simple_stem)

        # 2. DaCHS flattened stem (hermes_data_hermes_...)
        flat_stem = Path(s.replace("/", "_")).stem
        if flat_stem and flat_stem not in stems:
            stems.append(flat_stem)

        # 3. Original spectrum basename stem (00985062_HRF_OBJ_...)
        base_stem = Path(s).stem
        if base_stem and len(base_stem) >= 6 and base_stem not in stems:
            stems.append(base_stem)

    return stems


def _find_existing_file(row, existing_files):
    """Return the local Path if this spectrum already exists on disk, else None."""
    stems = _extract_row_stems(row)
    if not stems:
        return None

    for fpath in existing_files.values():
        local_stem = fpath.stem
        for key in stems:
            if local_stem == key or (len(key) >= 8 and local_stem.endswith(key)):
                return fpath
    return None


def _get_native_fits_url(row):
    """Check DataLink for a native FITS download URL before falling back to getdataurl()."""
    try:
        dl = row.getdatalink()
        for link in dl:
            url = str(link["access_url"] or "")
            ctype = str(link["content_type"] or "").lower()
            if url and ("fits" in ctype or url.lower().split("?")[0].endswith(".fits")):
                return url
    except Exception:
        pass
    return row.getdataurl()


def _download_mercator_spectrum(
    row, target_dir, fallback_name, existing_files, overwrite=False, timeout=60
):
    """
    Download a HERMES spectrum as '<ssa_targname>_<YYYYMMDD>_<HHMMSS>_raw.fits'.
    Prefers native FITS via DataLink; converts VOTable to FITS if the server returns .vot.
    """
    if not overwrite:
        existing_match = _find_existing_file(row, existing_files)
        if existing_match is not None:
            return existing_match, "skipped"

    filename = _build_mercator_filename(row)
    if filename is not None and not overwrite:
        if filename in existing_files:
            return existing_files[filename], "skipped"
        dest = target_dir / filename
        if dest.exists() and dest.stat().st_size > 0:
            existing_files[filename] = dest
            return dest, "skipped"

    url = _get_native_fits_url(row)
    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()

        if filename is None:
            cd = r.headers.get("Content-Disposition", "")
            match = re.search(r'filename="?([^";]+)"?', cd)
            if match:
                raw_name = Path(match.group(1)).name
            else:
                raw_name = unquote(url.split("?")[0])

            # Also pass ssa_pubdid if raw_name doesn't contain the date/reduced path
            name_source = raw_name
            if not re.search(r"\d{8}[/\\_]reduced[/\\_]\d+", name_source):
                pubdid = _get_row_value(row, "ssa_pubdid")
                if pubdid:
                    name_source = str(pubdid)

            filename = _simplify_hermes_filename(name_source, fallback_name=fallback_name)

        if not overwrite and filename in existing_files:
            return existing_files[filename], "skipped"

        dest = target_dir / filename
        if not overwrite and dest.exists() and dest.stat().st_size > 0:
            existing_files[filename] = dest
            return dest, "skipped"

        tmp_path = target_dir / f".{dest.stem}.tmp"
        with open(tmp_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)

    # Inspect the downloaded file header: native FITS starts with b"SIMPLE  ="
    with open(tmp_path, "rb") as f:
        magic = f.read(16)

    if magic.startswith(b"SIMPLE"):
        tmp_path.rename(dest)
    else:
        # Convert VOTable to FITS and remove temporary file
        try:
            vot_table = Table.read(tmp_path, format="votable")
            vot_table.write(dest, format="fits", overwrite=True)
        finally:
            if tmp_path.exists():
                tmp_path.unlink()

    existing_files[filename] = dest
    return dest, "downloaded"


def convert_vot_to_fits(folder, delete_vot=False):
    """
    Convert existing Mercator .vot files in `folder` (recursively) to simplified .fits files.

    Parameters
    ----------
    folder : str or Path
        Directory containing .vot files.
    delete_vot : bool, optional
        If True, deletes each .vot file after successful conversion (default: False).
    """
    root = Path(folder)
    vot_files = sorted(root.rglob("*.vot"))
    if not vot_files:
        print(f"No .vot files found in {root}")
        return

    for vot_path in vot_files:
        fits_name = _simplify_hermes_filename(vot_path.name)
        fits_path = vot_path.with_name(fits_name)
        try:
            tbl = Table.read(vot_path, format="votable")
            tbl.write(fits_path, format="fits", overwrite=True)
            print(f"Converted: {vot_path.name} -> {fits_path.name}")
            if delete_vot:
                vot_path.unlink()
        except Exception as e:
            print(f"Failed to convert {vot_path.name}: {e}")


def _resolve_target(target):
    """Resolve a (ra, dec) tuple, coordinate string, or star name to SkyCoord (ICRS)."""
    if _is_coord_pair(target):
        return SkyCoord(ra=target[0] * u.deg, dec=target[1] * u.deg, frame="icrs")

    if isinstance(target, str):
        for unit in ((u.hourangle, u.deg), (u.deg, u.deg)):
            try:
                return SkyCoord(target, unit=unit, frame="icrs")
            except ValueError:
                pass

        try:
            return SkyCoord.from_name(target, frame="icrs")
        except Exception:
            # Handle Durchmusterung names without a space (e.g. "BD+56473" -> "BD+56 473")
            spaced = re.sub(
                r"^(BD|CD|CPD)([+-]\d{2})(\d+)$",
                r"\1\2 \3",
                target.strip(),
                flags=re.IGNORECASE,
            )
            if spaced != target.strip():
                return SkyCoord.from_name(spaced, frame="icrs")
            raise

    raise ValueError(f"Unsupported target format: {target}")


def _target_folder_name(target):
    """Create a filesystem-safe folder name for a target."""
    if _is_coord_pair(target):
        return f"ra{target[0]:.5f}_dec{target[1]:.5f}"
    return re.sub(r"[^\w\-_+.]", "_", str(target).strip())


def _index_existing_files(search_dirs):
    """Recursively scan search_dirs and map filename -> Path for all non-empty files."""
    existing = {}
    for d in search_dirs:
        p = Path(d)
        if not p.exists():
            continue
        for f in p.rglob("*"):
            if f.is_file() and f.name != "metadata.csv" and f.stat().st_size > 0:
                existing.setdefault(f.name, f)
    return existing


def _save_metadata_csv(astropy_table, csv_path):
    """Save an Astropy Table to CSV, flattening any multidimensional VO columns."""
    clean = Table()
    for colname in astropy_table.colnames:
        col = astropy_table[colname]
        if col.ndim > 1:
            clean[colname] = [str(val) for val in col]
        else:
            clean[colname] = col
    clean.write(csv_path, format="csv", overwrite=True)


def _format_date(val):
    """Convert MJD floats to readable ISO timestamps; leave strings intact."""
    if np.ma.is_masked(val) or val is None:
        return "-"
    if isinstance(val, (int, float, np.integer, np.floating)):
        if 40000 < float(val) < 80000:
            return Time(float(val), format="mjd").iso[:19]
        return f"{float(val):.4f}"
    return str(val)[:19]


def _build_summary_table(result, raw_table, existing_files):
    """Extract key columns from the SSAP result into a compact, readable Table."""
    col_map = {c.lower(): c for c in raw_table.colnames}
    summary = Table()
    summary["#"] = np.arange(1, len(raw_table) + 1)

    if "ssa_targname" in col_map:
        summary["Target"] = [str(v) for v in raw_table[col_map["ssa_targname"]]]

    if "ssa_dateobs" in col_map:
        raw_dates = raw_table[col_map["ssa_dateobs"]]
        summary["Obs Date (UTC)"] = [_format_date(v) for v in raw_dates]

    if "ssa_timeext" in col_map:
        vals = raw_table[col_map["ssa_timeext"]]
        if not np.all(np.ma.getmaskarray(vals)):
            summary["Exp (s)"] = [
                "-" if np.ma.is_masked(v) else f"{float(v):.1f}" for v in vals
            ]

    if "ssa_snr" in col_map:
        vals = raw_table[col_map["ssa_snr"]]
        if not np.all(np.ma.getmaskarray(vals)):
            summary["SNR"] = [
                "-" if np.ma.is_masked(v) else f"{float(v):.1f}" for v in vals
            ]

    if "ssa_dstitle" in col_map:
        summary["Title"] = [str(v) for v in raw_table[col_map["ssa_dstitle"]]]
    elif "ssa_pubdid" in col_map:
        summary["ID"] = [str(v).split("/")[-1] for v in raw_table[col_map["ssa_pubdid"]]]

    summary["Local"] = [
        "Yes" if _find_existing_file(row, existing_files) is not None else "No"
        for row in result
    ]

    return summary


def get_mercator(
    targets,
    output_dir=None,
    radius_arcsec=15.0,
    download_spectra=True,
    existing_dirs=None,
    overwrite=False,
    return_table=False,
):
    """
    Search and download HERMES spectra from the Mercator Observatory SSAP archive.

    Parameters
    ----------
    targets : str, Path, tuple of (ra, dec), or list
        A star name (e.g. "V356 Per", "BD+56473"), a (ra_deg, dec_deg) tuple,
        a list mixing both, or a .lst/.txt/.fits file with target names.
    output_dir : str or Path, optional
        Output directory for metadata CSVs and downloaded spectra.
        Defaults to `<datadir>/HERMES/raw`.
    radius_arcsec : float, optional
        Cone search radius in arcseconds (default: 15.0).
    download_spectra : bool, optional
        If True, downloads the spectra into output_dir (default: True).
    existing_dirs : str, Path, or list of paths, optional
        Additional local directories to scan for already-downloaded spectra.
        `output_dir` is always scanned automatically.
    overwrite : bool, optional
        If True, re-download spectra even if they already exist locally (default: False).
    return_table : bool, optional
        If True, returns the full Astropy Table (single target) or dict of Tables
        (multiple targets). Default is False for clean terminal output.
    """
    if output_dir is None:
        out_path = Path(datadir).joinpath("HERMES/raw")
    else:
        out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    try:
        target_list = _load_targets_input(targets)
    except KeyError as e:
        print(e)
        return

    scan_dirs = [out_path]
    if existing_dirs is not None:
        if isinstance(existing_dirs, (str, Path)):
            scan_dirs.append(Path(existing_dirs))
        else:
            scan_dirs.extend(Path(d) for d in existing_dirs)
    existing_files = _index_existing_files(scan_dirs)

    service = pyvo.dal.SSAService(SSA_URL)
    diameter = (2.0 * radius_arcsec) * u.arcsec

    single_input = (
        (isinstance(targets, str) and not targets.endswith((".lst", ".txt", ".fits")))
        or _is_coord_pair(targets)
    )
    results_by_target = {}

    for target in target_list:
        if not _is_coord_pair(target):
            target = str(target).strip()

        print(f"\n=== {target} ===")
        try:
            coord = _resolve_target(target)
        except Exception as e:
            print(f"Could not resolve '{target}': {e}")
            continue

        ra_deg, dec_deg = coord.ra.deg, coord.dec.deg

        try:
            result = service.search(pos=coord, diameter=diameter)
        except Exception as e:
            print(f"SSAP query failed for '{target}': {e}")
            continue

        n_found = len(result)
        print(
            f"Found {n_found} spectra within {radius_arcsec}\" "
            f"of (RA={ra_deg:.5f}, Dec={dec_deg:.5f})"
        )

        table = result.to_table()
        results_by_target[target] = table

        if n_found == 0:
            continue

        summary = _build_summary_table(result, table, existing_files)
        summary.pprint_all()

        target_dir = out_path / _target_folder_name(target)
        target_dir.mkdir(parents=True, exist_ok=True)

        csv_path = target_dir / "metadata.csv"
        _save_metadata_csv(table, csv_path)
        print(f"\nSaved full metadata to {csv_path}")

        if download_spectra:
            downloaded_count = 0
            skipped_count = 0
            for idx, row in enumerate(result, start=1):
                if not row.getdataurl():
                    continue

                fallback_name = f"hermes_spectrum_{idx:03d}.vot"
                try:
                    path, status = _download_mercator_spectrum(
                        row,
                        target_dir,
                        fallback_name,
                        existing_files=existing_files,
                        overwrite=overwrite,
                    )
                    if status == "downloaded":
                        downloaded_count += 1
                        print(f"  [{idx}/{n_found}] Downloaded: {path.name}")
                    else:
                        skipped_count += 1
                        loc = "" if path.parent == target_dir else f" (in {path.parent})"
                        print(f"  [{idx}/{n_found}] Already exists, skipped: {path.name}{loc}")
                except Exception as e:
                    print(f"  [{idx}/{n_found}] Failed to download: {e}")

            print(
                f"Done: {downloaded_count} downloaded, {skipped_count} skipped "
                f"-> {target_dir}"
            )

    if return_table:
        return next(iter(results_by_target.values()), None) if single_input else results_by_target


# =============================================================================
# ESO / FEROS main function
# =============================================================================

def get_feros(fname, radius=2, program=None, limit=9999, output_dir=None):
    """
    Search and download FEROS spectra from the ESO archive.

    Parameters
    ----------
    fname : str, Path, or list
        Filename containing source names (.lst/.txt/.fits), a single star name,
        or a list of star names.
    radius : int/float, optional
        Max distance in arcmin to the queried result sources. Default is 2 arcmin.
    program : str, optional
        Program ID to search within the ESO archive. Default is None.
    limit : int, optional
        Maximum number of datasets to download per source. Default is 9999.
    output_dir : str or Path, optional
        Output directory for downloaded FITS files. Defaults to `<datadir>/FEROS/raw`.

    Returns
    -------
    None (downloads the FITS files and logs to downloaded_sources.txt).
    """
    try:
        sources = _load_targets_input(fname)
    except KeyError as e:
        return print(str(e))

    if len(sources) == 0:
        return msg.info("No sources provided.")

    eso = _get_eso_client()

    feros_dir = Path(output_dir) if output_dir is not None else Path(datadir).joinpath("FEROS/raw")
    feros_dir.mkdir(parents=True, exist_ok=True)

    msg.info(f"Files will be downloaded to {feros_dir}/")

    # Scale maxval by 100 per source so the bar can advance after every dataset
    maxval = len(sources) * 100
    bar = pb.ProgressBar(
        maxval=maxval,
        term_width=80,
        redirect_stdout=True,
        widgets=[pb.Bar("=", "[", "]"), " ", pb.Percentage()],
    )
    bar.start()

    with open(feros_dir.joinpath("downloaded_sources.txt"), "a") as db_log:
        for i, source in enumerate(sources):
            name = str(source).strip()
            print(f"\nSearching data for {name}...\n")

            query = eso.query_surveys(surveys="FEROS", cache=False, target=name)
            if query is None:
                db_log.write(f"{name}, 0, 0, Query returned no results.\n")
                bar.update((i + 1) * 100)
                continue

            if program is not None and "Run/Program ID" in query.colnames:
                query = query[query["Run/Program ID"] == program]

            simbad = query_Simbad(name)
            if simbad is None:
                db_log.write(f"{name}, 0, 0, Simbad return no results.\n")
                bar.update((i + 1) * 100)
                continue

            coords_S = SkyCoord(
                simbad["ra"], simbad["dec"], unit=(u.hourangle, u.deg), frame="icrs"
            )

            datasets = []
            for entry in query:
                try:
                    if re.sub("[ /_/]", "", str(entry["Object"])).replace("HD-", "HD") == name:
                        datasets.append(entry["ARCFILE"])
                    else:
                        coords_F = SkyCoord(
                            entry["RA"], entry["DEC"], unit=(u.deg, u.deg), frame="icrs"
                        )
                        if coords_S.separation(coords_F).arcmin[0] < radius:
                            datasets.append(entry["ARCFILE"])
                except Exception:
                    db_log.write(f"Source {name} had problems during the query.\n")
                    continue

            if len(datasets) == 0:
                db_log.write(f"{name}, 0, 0, No matches in FEROS database.\n")
                bar.update((i + 1) * 100)
                continue

            if len(datasets) > limit:
                db_log.write(f"{name}, 0, 0, >{limit} datasets for the source.\n")
                bar.update((i + 1) * 100)
                continue

            spectra = 0
            n_ds = len(datasets)
            for j, dataset in enumerate(datasets, start=1):
                miliseconds = round(float(dataset[-5:]) + 0.001, 3)

                if miliseconds == 1.0:
                    db_log.write(f"Please check: data from {name} is probably missing.\n")

                miliseconds_str = f"{miliseconds:.3f}"
                dataset_f = dataset[:22] + miliseconds_str

                print(f"[{j}/{n_ds}] {dataset_f} ({dataset})")
                if feros_dir.joinpath(dataset_f + ".fits").exists():
                    db_log.write(f"Dataset {dataset_f} already exists. Skipping download.\n")
                    spectra += 1
                    bar.update(min(i * 100 + int(100 * j / n_ds), maxval))
                    continue

                try:
                    eso.retrieve_data(dataset_f, destination=feros_dir)
                except Exception:
                    db_log.write(f"ERROR: Dataset {dataset_f} could not be retrieved.\n")
                    bar.update(min(i * 100 + int(100 * j / n_ds), maxval))
                    continue

                xml_file = feros_dir.joinpath(dataset_f + ".xml")
                if xml_file.exists():
                    xml_file.unlink()

                folder_name = feros_dir.joinpath(dataset_f)
                tgz_file = Path(str(folder_name) + ".tgz")

                with tarfile.open(tgz_file) as tar:
                    tar.extractall(folder_name)
                tgz_file.unlink()

                fits_files = [
                    f.name
                    for f in folder_name.iterdir()
                    if f.is_file() and f.suffix == ".fits"
                ]
                for fit_file in fits_files:
                    if not fit_file.endswith("1081.fits"):
                        folder_name.joinpath(fit_file).unlink()
                    else:
                        spectra += 1
                        print(f"Found a valid spectrum: {fit_file}")
                        shutil.move(
                            folder_name.joinpath(fit_file),
                            feros_dir.joinpath(fit_file),
                        )

                shutil.rmtree(folder_name)
                bar.update(min(i * 100 + int(100 * j / n_ds), maxval))

            db_log.write(f"{name}, {len(datasets)}, {spectra},--\n")

    bar.finish()
    return msg.g("\n Data retrieval completed. Check the folder.")