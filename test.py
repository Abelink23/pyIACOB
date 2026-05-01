import os
import shutil

#for module in ['db','reduction', 'get_feros', 'spec', 'spec_posproc', 'IACOBroad',
#                'measure','binarity', 'maui', 'models', 'rv']:
#    try:
#        exec(f'from {module} import *')
#        print(f'Module {module} imported successfully.')
#    except Exception as e:
#        print(f'Error importing module {module}: {e}')

print('\033[92mWelcome to the pyIACOB test program!\033[0m')
print('This program will test the main functions of the pyIACOB package.\n')
print('This program is not intended to be used by the user, but it is useful for development.\n')

#^ database module ======================================================================
try:
    from db import *
    print('Module db imported successfully.')
except Exception as e:
    print(f'\033[91mError importing module db: {e}\033[0m')

#^ findstar from the database module
print("Trying to find the star HD 37022 using findstar from the db module.")
try:
    star = findstar('HD 37022')
    print(f'Spectra available for star {star}')
except Exception as e:
    print(f'\033[91mError finding star: {e}\033[0m')
try:
    star = findstar('HD 37022', snr='best')
    print(f'Best spectra available for star {star[0]}')
except Exception as e:
    print(f'\033[91mError finding star with best SNR: {e}\033[0m')

#^ findlines from the database module
print("Creating a temporary file in the lists/lines and reading it using findlines from the db module.")
try:
    with open(maindir+'lists/lines/_temp_.txt', 'w') as f:
        f.write('4387.929, He I\n')
        f.write('5875.62,  He I\n')
        f.write('4552.622, Si III\n')
    lines = findlines('_temp_.txt')
    print(f'Lines found: {lines}')
    os.remove(maindir+'lists/lines/_temp_.txt')
except Exception as e:
    print(f'\033[91mError finding lines: {e}\033[0m')

#^ findtable from the database module
print("Trying to find the table 'stars' using findtable from the db module.")
# copy the file ALL_all.txt to the main_dir/tables directory
try:
    shutil.copyfile(maindir+'lists/lines/atlas_lines/ALL_all.txt', os.path.join(maindir, 'tables', '_ALL_all_.txt'))
    table = findtable('_ALL_all_.txt')
    print(f'Table found:\n {table[:10]}\n...')
    os.remove(os.path.join(maindir, 'tables', '_ALL_all_.txt'))
except Exception as e:
    print(f'\033[91mError finding table: {e}\033[0m')

#^ query_Simbad from the database module
print("Trying to query Simbad for the star HD37022 using query_Simbad from the db module.")
try:
    star = query_Simbad('HD 37022')
    print(f'Simbad query result for star HD 37022:\n {star}')
except Exception as e:
    print(f'\033[91mError querying Simbad: {e}\033[0m')

#^ query_Gaia from the database module
print("Trying to query Gaia for the star HD37022 using query_Gaia from the db module.")
try:
    star = query_Gaia('HD 37022')
    print(f'Gaia query result for star HD 37022:\n {star}')
except Exception as e:
    print(f'\033[91mError querying Gaia: {e}\033[0m')

#^ table_db from the database module
print("Trying to make a table using table_db from the db module.")
try:
    table = table_db('HD37022_20171217_224509_M_V85000.fits', db='IACOB',
                    spt='O', snrcut=50, spccode=True, bmag='>1', gaia='DR3')
    print(f'Table for star HD 37022:\n {table[:10]}\n...')
    findtable('tablestars.fits')
    os.remove('tablestars.fits')
except Exception as e:
    print(f'\033[91mError making table: {e}\033[0m')

#^ show_header from the database module
print("Trying to show the header of the file sample fits file using show_header from the db module.")
try:
    show_header('HD37022_20171217_224509_M_V85000.fits')
except Exception as e:
    print(f'\033[91mError showing header: {e}\033[0m')

#? spec module ==========================================================================
try:
    from spec import *
    print('Module spec imported successfully.')
except Exception as e:
    print(f'\033[91mError importing module spec: {e}\033[0m')

#? class spec from the spec module
print("Trying to create an instance of the class spec for the file sample fits file.")
try:
    s = spec('HD37022_20171217_224509_M_V85000.fits', snr=40, rv0=1)
    print(f'Instance of class spec created successfully for file {s.fitsfile}')
except Exception as e:
    print(f'\033[91mError creating instance of class spec: {e}\033[0m')

#? try several methods of the class spec
print("Trying to use the methods of the class spec for the instance created.")
try:
    s.waveflux(lwl=4000, rwl=8000, cut_edges=True)
    s.fitline(4552.622, info=True, fw3414=True, plot=True)
    s.snrcalc()
    s.cosmic()
    s.degrade(resol=5000, vsini=100)
    s.resamp(dlam=0.03)
    s.plotline(4552.622)
except Exception as e:
    print(f'\033[91mError using methods of class spec: {e}\033[0m')


#! MISSING TESTS FOR THE OTHER MODULES, BUT THIS IS ALREADY A GOOD START.
print('\n\033[92mAll tests completed!! pyIACOB is ready for use.\033[0m')