from db import *

from astroquery.nist import Nist
from astroquery.atomic import AtomicLineList


def atline(lwl, rwl, unit, elements=None, source='ALL'):

    '''
    Function to retrieve spectral lines from the Atomic Line databases.
    Input and output is given in Air conditions.

    Example: atline(lwl=4550,rwl=4576,unit='aa',elements='Si III',source='ALL')

    Parameters
    ----------
    lwl : float/int
        Sets the start wavelength.

    rwl : float/int
        Sets the end wavelength.

    unit : str
        Units of the input wavelengths in lwl and rwl (nm/angstrom).

    elements : str, optional
        Enter a string with the elements associated to the lines.
        If 'OB' it will pick the most characteristic lines of OB stars.

    source : str, optional
        Choose between 'ALL' and 'NIST' databases. Default is 'ALL'

    Note
    ----
    The wavelength type is set to 'Air' and the wavelength accuracy to 20.

    Returns
    -------
    List with the spectral lines.
    '''

    if lwl > rwl:
        print('ERROR: left wavelength must be smaller than right wavelength. Exiting...')
        return None

    if unit in ['angstrom','aa','A']:
        wavelength_range = (lwl*u.AA,rwl*u.AA)
    elif unit in ['nanometer','nm']:
        wavelength_range = (lwl*u.nm,rwl*u.nm)
    else:
        print('ERROR: wrong input units for the wavelength range. Exitting...')
        return None

    source = source.upper()

    if elements == 'OB':
        elements = 'H I, He I-II, O I-IV, C I-IV, Ne I-III, Fe I-III, N I-IV, \
        Si I-IV, Mg I-IV, S I-IV, V I-II, Cr I-II, Ni I-II, Sc I-II, Ti I-II, \
        Ca I-II, Na I-II'

    if source == 'ALL':
        columns = ('spec', 'type', 'conf', 'term', 'angm', 'prob')
        query = AtomicLineList.query_object(wavelength_range=wavelength_range, wavelength_type='Air',
                wavelength_accuracy=20, element_spectrum=elements, output_columns=columns,
                cache=False)

        if 'A_ki' in query.colnames and query['A_ki'].dtype != 'float64':
            query['A_ki'] = [float(i.replace('None','nan')) for i in query['A_ki']]

    elif source == 'NIST':
        Nist.TIMEOUT = 60
        query = Nist.query(minwav=wavelength_range[0], maxwav=wavelength_range[1],
                wavelength_type='vac+air', linename=elements)

    return query


def find_configurations(table, source, element, ion, configuration=None, term=None, lwl=None, rwl=None):
    '''
    Function to find the configurations of a given element and ion.
    Example: find_configurations(table,'Si','III','3s.4p-3s.4d')

    Example: find_configurations(table=tabla,source='ALL',element='Si',ion='III',
                                    configuration='3s.4s-3s.4p',lwl=4550,rwl=4576)


    Parameters
    ----------
    table : astropy.table.Table
        Table with the atomic lines.

    source : str
        Source of the atomic lines, 'NIST' or 'ALL'.

    element : str
        Element of the atomic lines.

    ion : str
        Ion of the atomic lines.

    configuration : str, optional
        Configuration of the atomic lines.

    term : str, optional
        Term of the atomic lines.

    lwl : float/int, optional
        Sets the start wavelength. Default is None.

    rwl : float/int, optional
        Sets the end wavelength. Default is None.

    Returns
    -------
    List with the atomic lines.
    '''

    if source == 'NIST':
        table.rename_column('Ritz','wl')
        table.rename_column('Spectrum','spc')
        table['config'] = [i['Lower level'].split('|')[0].strip() + '-' + i['Upper level'].split('|')[0].strip() for i in table]
        table['term'] = [i['Lower level'].split('|')[1].strip().replace('*','o') + '-' + i['Upper level'].split('|')[1].strip().replace('*','o') for i in table]
    elif source == 'ALL':
        table.rename_column('LAMBDA AIR ANG','wl')
        table.rename_column('SPECTRUM','spc')
        table.rename_column('CONFIGURATION','config')
        table.rename_column('TERM','term')

    table['elem'] = [i.split(' ')[0] for i in table['spc']]
    table['ion'] = [i.split(' ')[1] for i in table['spc']]
    if configuration != None and term == None:
        sub_table = table[(table['elem'] == element) & (table['ion'] == ion) & (table['config'] == configuration)]
    elif term != None and configuration == None:
        sub_table = table[(table['elem'] == element) & (table['ion'] == ion) & (table['term'] == term)]
    elif term != None and configuration != None:
        sub_table = table[(table['elem'] == element) & (table['ion'] == ion) & (table['config'] == configuration) & (table['term'] == term)]

    sub_table.sort('term')
    terms = list(set(sub_table['term']))
    configs = list(set(sub_table['config']))

    if lwl != None and rwl != None:
        if lwl > rwl:
            print('ERROR: left wavelength must be smaller than right wavelength. Exiting...')
            return None
        sub_table = sub_table[(sub_table['wl'] > lwl) & (sub_table['wl'] < rwl)]
        if len(sub_table) == 0:
            print('WARNING: no lines found in the given wavelength range.')

    return sub_table, terms, configs


