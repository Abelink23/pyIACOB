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
print("Trying to find the star HD 37042 using findstar from the db module.")
try:
    star = findstar('HD 37042')
    print(f'\033[92mSpectra available for star {star}\033[0m')
except Exception as e:
    print(f'\033[91mError finding star: {e}\033[0m')
try:
    star = findstar('HD 37042', snr='best')
    print(f'Best spectra available for star {star[0]}')
except Exception as e:
    print(f'\033[91mError finding star with best SNR: {e}\033[0m')

#^ findlines from the database module
print("Creating a temporary file in the lists/lines and reading it using findlines from the db module.")
tmp_file_path = maindir / 'lists/lines/_temp_.txt'
try:
    with open(tmp_file_path, 'w') as f:
        f.write('4387.929, He I\n')
        f.write('5875.62,  He I\n')
        f.write('4552.622, Si III\n')
    lines = findlines('_temp_.txt')
    print(f'\033[92mLines found: {lines}\033[0m')
    tmp_file_path.unlink()
except Exception as e:
    print(f'\033[91mError finding lines: {e}\033[0m')

#^ findtable from the database module
print("Trying to find the table 'stars' using findtable from the db module.")
# copy the file all.txt to the main_dir/tables directory
try:
    shutil.copyfile(maindir/'lists/lines/atlas_lines/all.txt', maindir/'tables'/'_ALL_all_.txt')
    table = findtable('_ALL_all_.txt')
    print(f'\033[92mTable found:\n {table[:10]}\n...\033[0m')
    Path(maindir / 'tables' / '_ALL_all_.txt').unlink()
except Exception as e:
    print(f'\033[91mError finding table: {e}\033[0m')

#^ query_Simbad from the database module
print("Trying to query Simbad for the star HD37042 using query_Simbad from the db module.")
try:
    star = query_Simbad('HD 37042')
    print(f'\033[92mSimbad query result for star HD 37042:\n {star}\033[0m')
except Exception as e:
    print(f'\033[91mError querying Simbad: {e}\033[0m')

#^ query_Gaia from the database module
print("Trying to query Gaia for the star HD37042 using query_Gaia from the db module.")
try:
    star = query_Gaia('HD 37042')
    print(f'\033[92mGaia query result for star HD 37042:\n {star}\033[0m')
except Exception as e:
    print(f'\033[91mError querying Gaia: {e}\033[0m')

#^ table_db from the database module
print("Trying to make a table using table_db from the db module.")
try:
    table = table_db('HD37042_20171030_024023_M_V85000.fits', db='IACOB',
                    spt='B', snrcut=50, spccode=True, bmag='>1', gaia='DR3')
    print(f'\033[92mTable for star HD 37042:\n {table[:10]}\n...\033[0m')
    findtable('tablestars.fits')
    Path(maindir / 'tables' / 'tablestars.fits').unlink()
except Exception as e:
    print(f'\033[91mError making table: {e}\033[0m')


#? spec module ==========================================================================
try:
    from spec import *
    print(f'\033[92mModule spec imported successfully.\033[0m')
except Exception as e:
    print(f'\033[91mError importing module spec: {e}\033[0m')

#? class spec from the spec module
print("Trying to create an instance of the class spec for the file sample fits file.")
try:
    s = spec('HD37042_20171030_024023_M_V85000.fits', snr=40, rv0=1)
    print(f'\033[92mInstance of class spec created successfully for file {s.filename}\033[0m')
except Exception as e:
    print(f'\033[91mError creating instance of class spec: {e}\033[0m')

#? try several methods of the class spec
print("Trying to use the methods of the class spec for the instance created.")
try:
    s.waveflux(lwl=4000, rwl=8000, cut_edges=True)
    print('Close the plot window to continue...')
    s.fitline(4552.622, info=True, fw3414=True, plot=True)
    plt.show(block=True)
    s.snrcalc()
    s.cosmic()
    s.convolution(resol=5000, vsini=100)
    s.resamp(dlam=0.03)
    print(f'\033[92mMethods of class spec used successfully.\033[0m')
except Exception as e:
    print(f'\033[91mError using methods of class spec: {e}\033[0m')


#! MISSING TESTS FOR THE OTHER MODULES, BUT THIS IS ALREADY A GOOD START.
print('\n\033[92mAll tests completed!!\033[0m')