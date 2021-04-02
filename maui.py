from spec import *
from scipy.io.idl import readsav

grids = {
'nlte_10.1.6_SOLAR_expoclump_2019-10-24': ('BSgs_CNOSiMg_old',1),
'nlte_10.1.6_bdwarfs_SOLAR_2020-01-29': ('BDws_CNOSIMg_old',2),
'nlte_10.4.7_OB.Sg_SOLAR_2021-01-23': ('OBSgs_hot_NOSi_new',3),
'nlte_10.4.7_late.bsgs_SOLAR_expoclump_NOSi.djl_2021-02-06': ('BSgs_cool_NOSi_new',4)}

class idl():
    def __init__(self, idlfile):
        '''
        Parameters
        ----------

        idlfile : str
            Enter the input spectrum full path to the .idl file.

        NOTE: poner tooodos los errores y toooodas las cosas en la clase
         y luego cada vez que encuentre un .idl que sea fácil para ingestar una
         fila en la tabla por idl. De el programa input_maui podría escupir una
         lista o tabla que use aquí en gen_table para crear la tabla summary (py)
        '''

        idldata = readsav(idlfile)
        #for i in idldata.keys():
        #    try: print(i,len(idldata[i]))
        #    except: print(i)

        self.file_name = idldata.aa[0][3].decode()
        self.name_star = self.file_name.split('_')[0]
        self.resolution = int(self.file_name.split('_V')[-1][0:5])
        self.gridname = idldata.modelgridname.decode()
        self.synwave = idldata.xx_mod
        self.synflux = idldata.spec_prim
        self.dx = (idldata.xx_mod[-1]-idldata.xx_mod[0])/len(idldata.xx_mod)

        self.vsini = idldata.obsdat.spectrum[0].VSINI[0]
        self.vmac = idldata.obsdat.spectrum[0].MACRO[0]


def gen_ascii(solution_dir=maindir+'tables',
              save_dir=datadir+'ASCII/Synthetic_MAUI/',
              lwl=3900,rwl=5080):
    '''
    Parameters
    ----------

    solution_dir : str
        Enter the directory where to search for all the .idl file.

    save_dir : str
        Enter the directory where to save all the output spectra in ascii format.

    lwl : float, optional
        Sets the start wavelenght of the output spectrum.

    rwl : float, optional
        Sets the end wavelenght of the output spectrum.
    '''

    if solution_dir == 'def':
        solution_dir = '/net/nas/proyectos/hots/masblue/maui2021/RESULTS_BSGS_202101/SOLUTION/'
    if save_dir == 'def':
        save_dir = '/net/nas/proyectos/hots/masblue/obs_iac/spec_opt/IACOB_DB/ASCII/SYNTHETIC/'

    for file in os.listdir(solution_dir):
        if not file.startswith('._') and file.endswith('.idl') and not 'HERMES' in file:
            idlspec = idl(solution_dir+file)

            star_db = spec(idlspec.name_star,SNR='best')

            if idlspec.file_name != star_db.file_name[:-5]:
                print('\nWARNING: %s does not match with best spectrum available.' % idlspec.file_name[:-5])

            else:
                #idlspec.file_name = idlspec.file_name.replace(str(idlspec.resolution),'85000')
                new_idlspec = '%s_red%i.dat' % (idlspec.file_name,grids[idlspec.gridname][1])
                np.savetxt(save_dir+new_idlspec,np.c_[idlspec.synwave,idlspec.synflux],
                           fmt=('%.4f','%.6f'))

                star_idl = spec(new_idlspec,txt=True)
                star_idl.txtwaveflux(lwl,rwl)
                #plt.plot(star_idl.wave,star_idl.flux,'r',lw=.3) # plot to check
                star_idl.degrade(profile='rotmac',vsini=idlspec.vsini,vmac=idlspec.vmac)
                #plt.plot(star_idl.wave,star_idl.flux,'g',lw=.3) # plot to check

                np.savetxt(save_dir+new_idlspec,np.c_[star_idl.wave,star_idl.flux],
                           fmt=('%.4f','%.6f'))


def gen_table(solution_dir='def', tables_dir='def', input_table='MAUI_ver9.txt'):
    '''
    Parameters
    ----------

    solution_dir : str
        Enter the directory where to search for all the .idl file.

    tables_dir : str
        Enter the directory where to locate the input_table and output fits.

    input_table : str
        Name of the input table contaning the list of .idl files to search.
    '''

    if solution_dir == 'def':
        solution_dir = '/net/nas/proyectos/hots/masblue/maui2021/RESULTS_BSGS_202101/SOLUTION/'
    if tables_dir == 'def':
        tables_dir = '/net/nas/proyectos/hots/adeburgos/tables/'

    table = findtable(input_table)

    bar = pb.ProgressBar(maxval=len(table),
                         widgets=[pb.Bar('=','[',']'),' ',pb.Percentage()])
    bar.start()

    data_rows = []
    for row,i in zip(table,range(len(table))):
        file_name = row['filename']

        match = []
        for file in os.listdir(solution_dir):
            if file.endswith('.idl') and \
               file.split('_sqexp_mat1_')[1][:-14]+'RV.ascii' == file_name:
                   match.append(solution_dir+file)

        if len(match) == 0:
            print('\nWARNING: No .idl file found for %s. Continuing...' % file_name)
            continue

        elif len(match) > 1:
            match = sorted(match, key=lambda x: int(x[-14:-4].replace('-','')),reverse=True)

        idldata = readsav(match[0])

        file_name = idldata.aa[0][3].decode()
        name_star = file_name.split('_')[0]

        star_db = spec(name_star,SNR='best')

        if file_name != star_db.file_name[:-5]:
            print('\nWARNING: %s does not match with best spectrum available.' % file_name[:-5])

        data_row = []; data_row.extend([name_star])
        parameters = [j.decode() for j in idldata.solution.var_label[0]]
        for j in range(len(parameters)):
            # Los 11/13 parametros
            median = idldata.solution[0].sol[j]
            sol_max = idldata.solution[0].sol_max[0][j] # maximum of distribution

            err_sim = idldata.solution[0].err_sol[0][j] # symetric error of distribution at 67%
            err_dw = abs(sol_max-idldata.solution[0].hpd_interval[0][j])
            err_up = abs(sol_max-idldata.solution[0].hpd_interval[1][j])

            if err_up < err_sim/2: label1 = '<'
            elif err_up > err_sim/2: label1 = '>'
            elif err_up == err_sim/2: label = '='
            if err_dw < err_sim/2: label2 = '<'
            elif err_dw < err_sim/2: label2 = '>'
            elif err_dw == err_sim/2: label = '='
            data_row.extend([sol_max,err_dw,err_up])

        if len(parameters) == 11:
            data_row.extend([np.nan,np.nan,np.nan,np.nan,np.nan,np.nan])
        data_rows.append(tuple(data_row))

        bar.update(i)
        time.sleep(0.1)

    bar.finish()

    output = Table(
    rows=data_rows,
    names=('Name',
    'Teff','Teff_eUP','Teff_eDW',
    'lgf','lgf_eUP','lgf_eDW',
    'He','He_eUP','He_eDW',
    'Micro','Micro_eUP','Micro_eDW',
    'logQs','logQs_eUP','logQs_eDW',
    'beta','beta_eUP','beta_eDW',
    'C','C_eUP','C_eDW',
    'N','N_eUP','N_eDW',
    'O','O_eUP','O_eDW',
    'Mg','Mg_eUP','Mg_eDW',
    'Si','Si_eUP','Si_eDW',
    'fcl','fcl_eUP','fcl_eDW',
    'vcl','vcl_eUP','vcl_eDW'))

    hdu = fits.BinTableHDU(data=output.filled(np.nan))
    hdu.writeto(tables_dir+'MAUI_results.fits',overwrite=True)

    return('DONE')
