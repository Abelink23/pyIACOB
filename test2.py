import re
import time
import os.path
import platform

import warnings; warnings.filterwarnings("ignore")

import numpy as np

import astropy.units as u
from astropy.io import fits
from astropy.table import Table, join, setdiff, vstack, hstack
from astropy.coordinates import SkyCoord

f_dir = []; max = 0; min = 0
for root, dirs, files in os.walk("F:/backup_SSD/Documents/ML/dataset/spectra/"):
    for file in files:
        if file.endswith('.ascii'):
            f_dir = os.path.join(root,file)
            table =Table.read(f_dir,format='ascii',delimiter=' ')
            if max < table['col2'].max(): max = table['col2'].max()
            if min > table['col2'].min(): min = table['col2'].min()

table['col2'].std()
