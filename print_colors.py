# This file contains the class color, which is use in the
# pyIACOB package to print colored text in the terminal.

class color:
    p = '\033[95m'
    c = '\033[96m'
    b = '\033[94m'
    g = '\033[92m'
    y = '\033[93m'
    r = '\033[91m'
    bold = '\033[1m'
    end = '\033[0m'
    _ = '\033[4m'
    error = r + 'ERROR: '
    warn = y + 'WARNING: '