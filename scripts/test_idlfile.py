import sys
sys.path.append('../')

from db import *

from scipy.io.idl import readsav

idldata = readsav('/Users/abelink/Documents/MAUI/SOLUTION/emulated_solution_mcmc_sqexp_mat1_HD2905_20091109_191245_N_V46000_2021-01-10.idl')
for i in idldata.keys():
    try: print(i,len(idldata[i]))
    except: print(i)

filename = idldata.aa[0][3].decode()
id = filename.split('_')[0]
resolution = int(filename.split('_V')[-1][0:5])
gridname = idldata.modelgridname.decode()
synwave = idldata.xx_mod
synflux = idldata.spec_prim
dx = (idldata.xx_mod[-1]-idldata.xx_mod[0])/len(idldata.xx_mod)

vsini = idldata.obsdat.spectrum[0].VSINI[0]
vmac = idldata.obsdat.spectrum[0].MACRO[0]

[i for i in idldata.solution.dtype.names]

idldata.solution[0].var_label
idldata.solution[0].sol

data_rows = []
data_row = []
data_row.extend([id])
parameters = [i.decode() for i in idldata.solution.var_label[0]]
for i in range(len(parameters)):
    # Los 11/13 parametros
    median = idldata.solution[0].sol[i]
    sol_max = idldata.solution[0].sol_max[0][i] # maximum of distribution

    err_sim = idldata.solution[0].err_sol[0][i] # symetric error of distribution at 67%
    err_dw = abs(sol_max-idldata.solution[0].hpd_interval[0][i])
    err_up = abs(sol_max-idldata.solution[0].hpd_interval[1][i])

    if err_up < err_sim/2: label1 = '<'
    elif err_up > err_sim/2: label1 = '>'
    elif err_up == err_sim/2: label = '='
    if err_dw < err_sim/2: label2 = '<'
    elif err_dw < err_sim/2: label2 = '>'
    elif err_dw == err_sim/2: label = '='
    print(label1,label2)
    data_row.extend([sol_max,err_dw,err_up])

if len(parameters) == 11:
    data_row.extend([np.nan,np.nan,np.nan,np.nan,np.nan,np.nan])
data_rows.append(tuple(data_row))
