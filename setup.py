import os
import sys
import re
import shutil
import subprocess

print('Welcome to the pyIACOB package!\n')

# Obtain the paths to the main and the data mandatoy directories.
while True:
    main_dir = input('Please provide the path to the main directory: ').strip()
    if main_dir == '':
        print('The main directory is mandatory, please provide the path to the main directory.')
    else:
        break

while True:
    data_dir = input('Please provide the path to the data directory: ').strip()
    if data_dir == '':
        print('The data directory is mandatory, please provide the path to the data directory.')
    else:
        break

# Obtain the paths to the optional directories.
maui_dir   = input('Please provide the path to the MAUI directory: ').strip()
ib_dir     = input('Please provide the path to the IACOB-Broad directory: ').strip()
models_dir = input('Please provide the path to the models directory: ').strip()
tess_dir   = input('Please provide the path to the TESS directory: ').strip()

# Create the file 'paths.txt' which contains the paths to the directories.
with open('paths.txt', 'w') as f:
    f.write('main=' + main_dir + '\n')
    f.write('data=' + data_dir + '\n')
    f.write('maui=' + maui_dir + '\n')
    f.write('ib=' + ib_dir + '\n')
    f.write('models=' + models_dir + '\n')
    f.write('tess=' + tess_dir + '\n')

# Create the subfolders inside the main_dir directory
print('Creating the subfolders inside the main directory...\n')
subfolders = ['lists', 'lists/lines', 'plots', 'tables', 'tmp']
for subfolder in subfolders:
    if not os.path.exists(os.path.join(main_dir, subfolder)):
        print(f"Creating {subfolder} directory...")
        os.mkdir(os.path.join(main_dir, subfolder))

# Create a subfolder 'ASCII' inside the data_dir directory if it does not exist.
if not os.path.exists(os.path.join(data_dir, 'ASCII')):
    print("Creating ASCII directory...")
    os.mkdir(os.path.join(data_dir, 'ASCII'))
    # move the files in test_data to the data_dir directory
    print("Moving the FITS files in test_data to the data directory...")
    for file in os.listdir('test_data'):
        shutil.move(os.path.join('test_data', file), os.path.join(data_dir, file))

# Copy the file 'snr_gaps.txt' to the main_dir/lists directory.
print('Copying the file snr_gaps.txt to the main/lists directory...\n')
shutil.copyfile('snr_gaps.txt', os.path.join(main_dir, 'lists', 'snr_gaps.txt'))

# Copy the 'atlas_lines' folder to the main_dir/lists/lines/ directory.
print('Copying the folder atlas_lines to the main/lists/lines directory...\n')
shutil.copytree('atlas_lines', os.path.join(main_dir, 'lists', 'lines', 'atlas_lines'), dirs_exist_ok=True)

# List of packages to check/install
packages = [
    'python>=3.12.2',
    'numpy==2.4.4',
    'matplotlib==3.10.8',
    'scipy==1.14.1',
    'astropy==7.2',
    'astroquery==0.4.11',
    'progressbar==2.5',
    ]

# Loop through packages and check if they are installed
print("Checking if all packages are installed...")
input_conda = ''
for package in packages:
    if package.startswith('python'):
        if sys.version_info >= (3, 12):
            print("python is already installed.")
            continue
        package_name = 'python'
    else:
        package_name = re.split(r'[=><]+', package)[0]

    try:
        if package_name != 'python':
            __import__(package_name)
            print(f"{package_name} is already installed.")
        else:
            raise ImportError
    except ImportError:
        print(f"{package_name} is not installed.")
        if input_conda == '':
            input_conda = input("Do you want to install the missing packages using conda? (y/n): ")
        if input_conda.lower() == 'y':
            try:
                subprocess.check_call(['conda', 'install', '-y', package])
            except subprocess.CalledProcessError:
                print(f"Failed to install {package_name} from default channels. Trying conda-forge...")
                subprocess.check_call(['conda', 'install', '-c', 'conda-forge', '-y', package])
        else:
            print(f"Please install {package_name} manually and run this setup script again.")
            exit()

print("\nAll packages are installed!")
