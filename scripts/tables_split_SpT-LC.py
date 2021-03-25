'''=============================================================================
        Script to split SpC column into SpT and LC from any csv table
============================================================================='''

import os
import csv

dir = raw_input('Enter the path where the file is located: (~/Desktop)')
if dir == '': dir = "/home/abelink/Desktop/"

fn = raw_input('Enter the name of the file: (*.csv)')
if fn == '': print 'No file was selectec'
elif fn[-4:] != ".csv": fn = fn + ".csv"

for root, dirs, files in os.walk(dir):
    for file in files:
        if file == fn:
            fn_dir = os.path.join(root, file)
            with open(fn_dir, 'r') as csv_file:
                data = csv.reader(csv_file, delimiter = ',')

                header = [i.strip(' ') for i in data.next()]
                header_spt = header.index('SPTYPE')
                header.insert(header_spt, 'SPT')
                header.insert(header_spt + 1, 'LC')
                header.remove('SPTYPE')

                with open(fn_dir[:-4]+'_new.csv','w') as output:
                        new_fn = csv.writer(output, delimiter = ',', lineterminator = '\n')

                        new_fn.writerow(header)
                        for row in data:
                            if row == []: new_fn.writerow(row); continue
                            elif 'I' in row[header_spt]:
                                spt = row[header_spt][:row[header_spt].index('I')]
                                lc  = row[header_spt][row[header_spt].index('I'):]
                            elif 'V' in row[header_spt]:
                                spt = row[header_spt][:row[header_spt].index('V')]
                                lc  = row[header_spt][row[header_spt].index('V'):]
                            else: continue
                            row.insert(header_spt, spt)
                            row.insert(header_spt +1, lc)
                            row.remove(row[header_spt +2])
                            new_fn.writerow(row)
csv_file.close(); output.close()
