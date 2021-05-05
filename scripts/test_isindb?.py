import sys
sys.path.append('../')

from db import *

table = findtable('SUGUS-SIMBAD-SOUTH-V9-B0-9-I-II.csv')

to_search = []
for source in table:
    query_simbad = Simbad.query_objectids(source['identifier'])

    names = [i.replace(' ','') for i in query_simbad['ID'] if
    'HD' in i or 'BD' in i or 'ADS' in i or 'CDS' in i or 'CD' in i]
    print(names)

    if any([findstar(name) is not None for name in names]) is False:
        to_search.append(source['identifier'])

print(to_search)
