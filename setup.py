import re
import sys
import shutil
import subprocess
from pathlib import Path

print('\033[96m\033[1m\033[4mWelcome to the pyIACOB package!\033[0m\n')
print('\033[96mDeveloped by Abel de Burgos and the IACOB team.\033[0m\n')
print('\033[96mThis setup script will help you configure the pyIACOB package focused on the IACOB-broad sub-module and its dependencies.\033[0m\n')

# Obtain the paths to the main and data mandatory directories.
while True:
    raw_input = input('Please provide the path to the main working directory: ').strip()
    if not raw_input:
        print('The main directory is mandatory, please provide the path to the main directory.')
        continue

    path_candidate = Path(raw_input).expanduser()
    if not path_candidate.exists():
        create = input('The main directory does not exist, do you want to create it? (y/n): ')
        if create.lower() == 'y':
            path_candidate.mkdir(parents=True, exist_ok=True)
            main_dir = path_candidate
            break
        else:
            print('Please provide a valid path to the main directory.')
    else:
        main_dir = path_candidate
        break

while True:
    raw_input = input('Please provide the path to the data directory: ').strip()
    if not raw_input:
        print('The data directory is mandatory, please provide the path to the data directory.')
        continue

    path_candidate = Path(raw_input).expanduser()
    if not path_candidate.exists():
        create = input('The data directory does not exist, do you want to create it? (y/n): ')
        if create.lower() == 'y':
            path_candidate.mkdir(parents=True, exist_ok=True)
            data_dir = path_candidate
            break
        else:
            print('Please provide a valid path to the data directory.')
    else:
        data_dir = path_candidate
        break

# Create the file 'paths.txt' which contains the paths to the directories.
with open('paths.txt', 'w') as f:
    f.write('main=' + str(main_dir) + '\n')
    f.write('data=' + str(data_dir) + '\n')

# Create the subfolders inside the main_dir directory
print('\nCreating the subfolders inside the main directory...\n')
subfolders = ['broadening', 'lists', 'lists/lines', 'plots', 'tables', 'tmp']
for subfolder in subfolders:
    target_folder = main_dir / subfolder
    if not target_folder.exists():
        print(f"Creating {subfolder} directory...")
        target_folder.mkdir(parents=True, exist_ok=True)

# Create a subfolder 'ASCII' inside the data_dir directory if it does not exist.
ascii_dir = data_dir / 'ASCII'
if not ascii_dir.exists():
    print("Creating ASCII directory...")
    ascii_dir.mkdir(parents=True, exist_ok=True)
    # Move the files in sample_spectra to the data_dir directory
    print("Moving the FITS files in sample_spectra to the data directory...")
    if Path('sample_spectra').exists():
        for file in Path('sample_spectra').iterdir():
            if file.is_file():
                shutil.move(str(file), str(data_dir / file.name))
                print(f"Moved {file.name} to {data_dir}")
    else:
        print("WARNING: sample_spectra directory does not exist. No files were moved.")

# Copy the file 'snr_gaps.txt' to the main_dir/lists directory.
print('Copying the file snr_gaps.txt to the main/lists directory...')
shutil.copyfile('snr_gaps.txt', main_dir / 'lists' / 'snr_gaps.txt')


# --- Environment Activation Confirmation & Manager Selection ---
print('\n' + '=' * 60)
print('DEPENDENCY INSTALLATION SETUP')
print('=' * 60)

# Acknowledge active environment
while True:
    ack = input(
        '\nHave you already activated the target environment where you want to install pyIACOB dependencies? (y/n): '
    ).strip().lower()
    if ack in ('y', 'yes'):
        break
    elif ack in ('n', 'no'):
        print('\nPlease activate your target virtual environment (Conda or venv) and re-run this setup script.')
        sys.exit(0)
    else:
        print("Please answer with 'y' or 'n'.")

# Choose environment manager
while True:
    env_choice = input(
        'Which environment manager are you using? (Enter "conda" or "venv"): '
    ).strip().lower()
    if env_choice in ('conda', 'venv'):
        break
    print("Invalid choice. Please enter either 'conda' or 'venv'.")

# Python version check
required_py = '3.11'
req_major, req_minor = map(int, required_py.split('.'))
if sys.version_info < (req_major, req_minor):
    print(f"\nError: Python >={required_py} is required, but you are running {sys.version.split()[0]}.")
    print("Please activate an environment with the required Python version and try again.")
    sys.exit(1)
else:
    print(f"\nPython >={required_py} requirement satisfied.")

# List of packages to check/install
packages = [
    'numpy==2.4.6', # 2.5.4 (py 3.13)
    'matplotlib==3.11.2',
    'scipy==1.17.1', # 1.18.1 (py 3.13)
    'astropy==8.0.2',
    'astroquery==0.4.11',
    'progressbar2==4.6.0', # 4.7.0 (py 3.13)
]

# Mapping of distribution package name -> Python import name
IMPORT_MAP = {
    'progressbar2': 'progressbar',
}

# Loop through packages and check if they are installed
print("\nChecking if all packages are installed...")
for package in packages:
    package_name = re.split(r'[=><]+', package)[0]
    module_to_import = IMPORT_MAP.get(package_name, package_name)

    try:
        __import__(module_to_import)
        print(f"✓ {package_name} is already installed.")
    except ImportError:
        print(f"✗ {package_name} is not installed.")

        confirm = input(f"Do you want to install {package} using {env_choice}? (y/n): ").strip().lower()
        if confirm not in ('y', 'yes'):
            print(f"Please install {package} manually and re-run the setup script.")
            sys.exit(1)

        print(f"Installing {package} using {env_choice}...")
        try:
            if env_choice == 'conda':
                try:
                    subprocess.check_call(['conda', 'install', '-y', package])
                except subprocess.CalledProcessError:
                    print(f"Failed to install {package_name} from default channels. Trying conda-forge...")
                    subprocess.check_call(['conda', 'install', '-c', 'conda-forge', '-y', package])
            else:  # venv
                # Using sys.executable guarantees package installation into the active Python interpreter
                subprocess.check_call([sys.executable, '-m', 'pip', 'install', package])
        except subprocess.CalledProcessError as err:
            print(f"Failed to install {package}: {err}")
            sys.exit(1)

print("\n\033[1m\033[4m\033[96mAll packages are installed!\033[0m")
