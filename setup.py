import os
import shutil
import subprocess

print('Welcome to the pyIACOB package!\n')

# Obtain the paths to the main and the data mandatoy directories.
while True:
    main_dir = input('Please provide the path to the main directory: ')
    if main_dir == '':
        print('The main directory is mandatory, please provide the path to the main directory.')
    else:
        break

while True:
    data_dir = input('Please provide the path to the data directory: ')
    if data_dir == '':
        print('The data directory is mandatory, please provide the path to the data directory.')
    else:
        break

# Obtain the paths to the optional directories.
maui_dir   = input('Please provide the path to the MAUI directory: ')
ib_dir     = input('Please provide the path to the IACOB-Broad directory: ')
models_dir = input('Please provide the path to the models directory: ')
tess_dir   = input('Please provide the path to the TESS directory: ')

# Create the file 'paths.txt' which contains the paths to the directories.
with open('paths.txt', 'w') as f:
    f.write('main=' + main_dir + '\n')
    f.write('data=' + data_dir + '\n')
    f.write('maui=' + maui_dir + '\n')
    f.write('ib=' + ib_dir + '\n')
    f.write('models=' + models_dir + '\n')
    f.write('tess=' + tess_dir)

# Create the subfolders inside the main_dir directory
subfolders = ['list', 'list/lines', 'plots', 'tables', 'tmp']
for subfolder in subfolders:
    if not os.path.exists(os.path.join(main_dir, subfolder)):
        os.mkdir(os.path.join(main_dir, subfolder))

# Copy the file 'snr_gaps.txt' to the main_dir/list directory.
shutil.copyfile('snr_gaps.txt', os.path.join(main_dir, 'list', 'snr_gaps.txt'))

# List of packages to check/install
packages = ['numpy', 'matplotlib', 'scipy', 'astropy', 'astroquery', 'progressbar', 'random', 'lightkurve']

# Loop through packages and check if they are installed
for package in packages:
    try:
        # Attempt to import the package
        __import__(package)
        print(f"{package} is already installed.")
    except ImportError:
        print(f"{package} is not installed. Installing...")
        # Install package using conda
        if package in ['astroquery', 'progressbar', 'lightkurve']:
            subprocess.call(['conda', 'install', '-c', 'conda-forge', package, '-y'])
        else:
            subprocess.call(['conda', 'install', package, '-y'])

print("All packages are installed!")
