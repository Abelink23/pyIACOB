import re
from pathlib import Path
from astropy.table import Table
import astropy.units as u

UNIT_MAP = {
    "-": None,
    "--": None,
    "---": None,
    "yr": u.yr,
    "Msun": u.M_sun,
    "Lsun": u.L_sun,
    "K": u.K,
    "mass frac.": u.dimensionless_unscaled,
    "Ms/yr": u.M_sun / u.yr,
    "g/cm3": u.g / (u.cm**3),
    "s-1": u.s**-1,
    "km/s": u.km / u.s,
    "Msun/yr": u.M_sun / u.yr,
    "10^53 g cm2/s": 1e53 * u.g * (u.cm**2) / u.s,
}

def clean_unit(unit_val):
    """Strips brackets, replaces Msun with M_sun, and drops dashes."""
    if unit_val is None:
        return None

    # Remove brackets
    u_str = re.sub(r'[\[\]]', '', str(unit_val)).strip()

    # Drop dash variants
    if u_str in {'', '—', '--', '---'}:
        return None

    # Replace Msun
    return u_str.replace('Msun', 'M_sun')


def keep_columns(table, columns):
    cols_to_remove = [col for col in table.colnames if col not in columns]
    if len(cols_to_remove) > 0:
        table.remove_columns(cols_to_remove)
    return table


def cleanup_Brott11(source_dir: str, metallicity='014', output_subdir: str="processed_tracks"):

    columns =  ['t','Mass','Teff','logL','R','log(Mdot)','logg','Vsurf','Prot','Vcrit',
                'Ge','eps(H)','eps(He)','eps(C)','eps(N)','eps(O)','eps(F)','eps(Si)',
                'eps(Fe)','sH1','cH1','cHe3','cHe4']

    source_path = Path(source_dir)
    output_path = source_path / output_subdir
    output_path.mkdir(parents=True, exist_ok=True)

    PATTERN = re.compile(r"^f(\d+)-(\d+)\.mw\.fits$", re.IGNORECASE)

    for item in source_path.iterdir():
        if item.is_file() and not item.name.startswith("."):
            try:
                table = Table.read(item, format="fits")
            except Exception:
                continue

            match = PATTERN.match(item.name)
            if not match:
                continue

            raw_mass, raw_vini = match.groups()

            # Zero-pad mass and velocity to 3 digits (e.g. 7 -> 007, 0 -> 000, 284 -> 284)
            mass = f"{int(raw_mass):03d}"
            vini = f"{int(raw_vini):03d}"

            # Construct new output filename
            out_name = f"M{mass}Z{metallicity}V{vini}Av00.fits"
            out_file = output_path / out_name

            # 1. Clean columns
            table = keep_columns(table, columns)

            # 2. Clean column units
            for col_name in table.colnames:
                if table[col_name].unit is not None:
                    table[col_name].unit = clean_unit(table[col_name].unit)

            # 3. Save modified table
            table.write(out_file, format="fits", overwrite=True)
            print(f"Processed and saved: {out_file.name}")



def cleanup_MIST(source_dir: str, metallicity='14', vinivcrit='04', av='00', output_subdir: str="processed_tracks"):
    '''
    Parameters
    ----------
    source_dir : str
        Path to the directory containing raw MIST evolutionary track files.

    metallicity : str, optional
        Metallicity string identifier to include in the output filenames.
        Use only three digits. Default is '14' (z=0.014).

    vinivcrit : str, optional
        Initial rotation rate (v/vcrit) string identifier to include in the
        output filenames. Use two difits. Default is '04'.

    av : str, optional
        Extinction Av value of the computed tracks.
        Use two digits. Default is '00'.

    output_subdir : str, optional
        Subdirectory name inside `source_dir` where the resulting FITS
        files will be saved (default is 'processed_tracks').
    '''

    columns =  ["star_age","star_mass","star_mdot",
                "he_core_mass","c_core_mass","o_core_mass",
                "log_L","log_L_div_Ledd","log_LH","log_LHe","log_LZ",
                "log_Teff","log_abs_Lgrav","log_R","log_g","log_surf_cell_z",
                "surf_avg_omega","surf_avg_v_rot","surf_avg_omega_crit",
                "surf_avg_omega_div_omega_crit","surf_avg_v_crit","surf_avg_v_div_v_crit",
                "surface_h1","surface_h2","surface_he3","surface_he4",
                "surface_c12","surface_c13","surface_n13","surface_n14","surface_n15",
                "surface_o14","surface_o15","surface_o16","surface_o17","surface_o18",
                "surface_mg23","surface_mg24","surface_mg25","surface_mg26",
                "surface_si27","surface_si28","surface_si29","surface_si30",
                "log_center_T","log_center_Rho","pp","cno","tri_alfa","burn_c","burn_n","burn_o",
                "phase"]

    source_path = Path(source_dir)
    output_path = source_path / output_subdir
    output_path.mkdir(parents=True, exist_ok=True)

    if len(metallicity) != 3:
        print("\033[31m Input metallicity should have three digits! \033[0m")
        return None

    if len(vinivcrit) != 2:
        print("\033[31mInput vinivcrit should have two digits! \033[0m")
        return None

    if len(av) != 2:
        print("\033[31m Input Av should have two digits! \033[0m")
        return None

    for item in source_path.iterdir():
        if not item.is_file() or item.name.startswith("."):
            continue

        # Skip already converted FITS files
        if item.suffix.lower() == ".fits":
            continue

        # Extract the 3-character mass prefix (e.g., '009' from '0090000M.track.eep')
        match = re.match(r"^(\d{3})", item.name)
        if not match:
            continue
        mass = match.group(1)

        try:
            # Column headers are on line 12 (0-indexed 11); data starts on line 13 (0-indexed 12)
            table = Table.read(
                item,
                format="ascii",
                header_start=11,
                data_start=12,
                comment="#",
            )
        except Exception:
            continue

        # 1. Clean columns
        table = keep_columns(table, columns)

        # 2. Construct new filename and save modified table to FITS
        out_name = f"M{mass}Z{metallicity}V{vinivcrit}Av{av}.fits"
        out_file = output_path / out_name
        table.write(out_file, format="fits", overwrite=True)
        print(f"Processed and saved: {out_file.name}")



def cleanup_Sciarini26a(source_dir: str, output_subdir: str="processed_tracks"):

    columns = {
        1: 'model', 2: 'age', 3: 'mass', 4: 'logL', 5: 'logTeff', 6: 'sH1', 7: 'sHe4',
        8: 'sHe3', 9: 'sC12', 10: 'sC13', 11: 'sN14', 12: 'sO16', 13: 'sO17', 14: 'sO18',
        17: 'Mcc', 19: 'logMdot', 20: 'loggc', 21: 'logTc', 22: 'cH1', 23: 'cHe4', 24: 'cHe3',
        25: 'cC12', 26: 'cC13', 27: 'cN14', 28: 'cO16', 29: 'cO17', 30: 'cO18', 39: 'Vsurf/Vcrit',
        40: 'Vsurf', 41: 'Vcenter', 42: 'Rpol/Req', 53: 'Veq_before_conv', 58: 'Veq_after_conv',
        62: 'Mdot', 64: 'I', 65: 'Ltot', 111: 'MSper', 112: 'OB', 113: 'RSG', 114: 'WR',
    }

    source_path = Path(source_dir)
    output_path = source_path / output_subdir
    output_path.mkdir(parents=True, exist_ok=True)

    # Sorted list of 0-based column indices to extract
    target_indices = sorted(col_num - 1 for col_num in columns.keys())

    for item in source_path.glob("*.dat"):
        if item.name.startswith("."):
            continue

        try:
            # Read headerless whitespace-delimited ASCII file
            raw_table = Table.read(item, format="ascii.no_header")
        except Exception as e:
            print(f"Skipping {item.name}: failed to read ({e})")
            continue

        # Verify the table has enough columns
        max_needed_idx = target_indices[-1]
        if len(raw_table.colnames) <= max_needed_idx:
            print(
                f"Skipping {item.name}: expected at least {max_needed_idx + 1} columns, "
                f"found {len(raw_table.colnames)}")
            continue

        # Select only desired columns by their positional names ('col1', 'col2', etc.)
        selected_cols = [raw_table.colnames[idx] for idx in target_indices]
        filtered_table = raw_table[selected_cols]

        # Rename columns to their designated labels
        for idx in target_indices:
            orig_col = raw_table.colnames[idx]
            new_col =  columns[idx + 1]
            filtered_table.rename_column(orig_col, new_col)

        # Replace V4 by V04Av00 and V0 by V00Av00
        out_name = f"{item.stem}.fits"
        if 'V0' in out_name:
            out_name = out_name.replace('V0','V00Av00')
        elif 'V4' in out_name:
            out_name = out_name.replace('V4','V04Av00')

        # Save to FITS
        out_file = output_path / out_name
        filtered_table.write(out_file, format="fits", overwrite=True)
        print(f"Processed: {item.name} -> {out_file.name}")



def cleanup_Gene12(source_dir: str, output_subdir: str="processed_tracks"):

    remove_columns = ['20Ne_surf','22Ne_surf','26Al_surf','20Ne_cen','22Ne_cen','26Al_cen']

    source_path = Path(source_dir)
    output_path = source_path / output_subdir
    output_path.mkdir(parents=True, exist_ok=True)

    for item in source_path.glob("*.dat"):
        if item.name.startswith("."):
            continue

        # 1. Read header and units
        with open(item, "r") as f:
            header_line = f.readline().rstrip("\n")
            units_line = f.readline().rstrip("\n")

        # 2. Fix duplicated column names
        raw_names = header_line.split()
        seen = {}
        unique_names = []
        for name in raw_names:
            if name in seen:
                seen[name] += 1
                unique_names.append(f"{name}_{seen[name]}")
            else:
                seen[name] = 0
                unique_names.append(name)

        # 3. Split units
        pattern = r"10\^53\s+g\s+cm2/s|mass\s+frac\.|\S+"
        raw_units = re.findall(pattern, units_line)

        try:
            # Read headerless whitespace-delimited ASCII file
            table = Table.read(item,format="ascii",names=unique_names,data_start=2,guess=False,)
        except Exception as e:
            print(f"Skipping {item.name}: failed to read ({e})")
            continue

        # Remove unwanted columns
        table.remove_columns([i for i in remove_columns if i in table.colnames])

        # 5. Asignar las unidades a cada columna
        for col_name, u_str in zip(table.colnames, raw_units):
            if u_str in UNIT_MAP:
                table[col_name].unit = UNIT_MAP[u_str]
            else:
                try:
                    table[col_name].unit = u.Unit(u_str)
                except ValueError:
                    table[col_name].unit = None

        # Save to FITS
        out_file = output_path / f"{item.stem}.fits"
        table.write(out_file, format="fits", overwrite=True)
        print(f"Processed: {item.name} -> {out_file.name}")



def cleanup_Gene26(source_dir: str, output_subdir: str="processed_tracks"):

    remove_columns = ['20Ne_surf','22Ne_surf','26Al_surf','20Ne_cen','22Ne_cen','26Al_cen']

    source_path = Path(source_dir)
    output_path = source_path / output_subdir
    output_path.mkdir(parents=True, exist_ok=True)

    for item in source_path.glob("*.dat"):
        if item.name.startswith("."):
            continue

        # 1. Read header and units
        with open(item, "r") as f:
            header_line = f.readline().rstrip("\n")
            units_line = f.readline().rstrip("\n")

        # 2. Fix duplicated column names
        raw_names = header_line.split()
        seen = {}
        unique_names = []
        for name in raw_names:
            if name in seen:
                seen[name] += 1
                unique_names.append(f"{name}_{seen[name]}")
            else:
                seen[name] = 0
                unique_names.append(name)

        # 3. Split units
        pattern = r"10\^53\s+g\s+cm2/s|mass\s+frac\.|\S+"
        raw_units = re.findall(pattern, units_line)

        try:
            # Read headerless whitespace-delimited ASCII file
            table = Table.read(item,format="ascii",names=unique_names,data_start=2,guess=False,)
        except Exception as e:
            print(f"Skipping {item.name}: failed to read ({e})")
            continue

        # Remove unwanted columns
        table.remove_columns([i for i in remove_columns if i in table.colnames])

        # 5. Asignar las unidades a cada columna
        for col_name, u_str in zip(table.colnames, raw_units):
            if u_str in UNIT_MAP:
                table[col_name].unit = UNIT_MAP[u_str]
            else:
                try:
                    table[col_name].unit = u.Unit(u_str)
                except ValueError:
                    table[col_name].unit = None

        # Replace V4 by V04Av00 and V0 by V00Av00
        out_name = f"{item.stem}.fits"
        if '_hydro' in out_name:
            out_name = out_name.replace('_hydro','Av00_hydro')
        elif '_magn' in out_name:
            out_name = out_name.replace('_magn','Av00_magn')

        # Save to FITS
        out_file = output_path / out_name
        table.write(out_file, format="fits", overwrite=True)
        print(f"Processed: {item.name} -> {out_file.name}")
