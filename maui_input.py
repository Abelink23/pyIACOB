import sys
sys.path.append('../')

import matplotlib.path as mpath

from spec import *
from RV import *
from mist import *

table = findtable('IACOB_O9BAs_SNR20.fits') # file where quality flags are
table_REF = findtable('O9BAs_RVEWFWs.fits') # file where RVs, EWs and FWs are
table_IB = findtable('IB_results_ver5.txt') # file where vsini and vmac are
results = findtable('zsummary_results.txt') # file with output from MAUI

def input_maui(table='IACOB_O9BAs_SNR20.fits',output_name='MAUI_verX',RV0tol=200,ascii_0='y'):

    table = findtable(table)

    maui_txt = open(maindir+'lists/%s.txt' % output_name,'a')
    maui_txt.write(
            "{:<40}".format('starname')+"{:<48}".format('filename')+"{:<6}".format('vrad')+\
            "{:<8}".format('vsini')+"{:<7}".format('evsini')+"{:<8}".format('rt_mac')+\
            "{:<7}".format('emac')+"{:<10}".format('R')+"{:<4}".format('SNR')+\
            ' ;# SpC       FW34-14 lgf lTf SiIII SiII\n')

    quit = ''
    for row in table:

        ascii = ascii_0

        if quit == 'quit': break

        name = row['Name'].strip()

        # Filer based on properties from the main table:
        if 'SB2' in row['SB']: continue
        if 'Em' in row['CHb']:
            if not 'Em(p)' in row['CHb']: continue
        if 'PCyg' in row['CHb']: continue
        if row['QIB'] < 2: continue

        match_REF = table_REF[[i.strip()==name for i in table_REF['Name']]]
        match_IB = table_IB[table_IB['Name']==name]

        if len(match_REF) == 0 or len(match_IB) == 0: continue

        # Filter based on Si lines properties:
        if row['QSiII']<3 or match_REF['EWSiII']<50 or np.isnan(match_REF['EWSiII']) or\
            match_REF['depSiII']<3/match_REF['SNR_B']: SiIIFG = 0
        else: SiIIFG = 1

        if row['QSiIII']<3 or match_REF['EWSiIII1']<50 or np.isnan(match_REF['EWSiIII1']) or\
            match_REF['depSiIII1']<3/match_REF['SNR_B']: SiIIIFG = 0
        else: SiIIIFG = 1

        if  SiIIIFG == 0:# or SiIIFG == 0:
            print('No SiIII found for %s\n' % name); continue

        star = spec(name,SNR='best')

        if match_IB['filename'][0] != star.file_name:
            print('Warning: Different files from best SNR and from IB results for %s' % name)
            print(star.file_name,' vs ',match_IB['filename'][0])

        # Extra information appended to the end of each row:
        match_results = results[results['Name']==name]
        if len(match_results) == 0: loggf = logTf = 0
        else: loggf = 5.39-match_results['lgf'][0]; logTf = 4+np.log10(match_results['Teff'][0])
        # --------------------------------------------------------------------------

        if ascii_0 == 'y' and not type(search(match_IB['filename'][0][:-5]+'_RV.ascii',\
        os.path.expanduser('~')+'/Documents/MAUI/ASCII/')) == type(None): ascii = 'n'

        if ascii == 'y':

            # If RV SiIII is good enough, it uses it for the offset:
            if abs(match_REF['RVSiIII1']-match_REF['RVHb']) < 10:
                star.offset = float(match_REF['RVSiIII1'])*float(4500)*1000/cte.c
                star.waveflux(); print('RV good enough.\n')

            # Otherwise it calculates the offset with the RV0 program:
            else:
                next = 'n'
                while next == 'n':

                    skip = input('%s - Hit return to continue, type "s" to skip: ' % name)
                    if skip == 's': break

                    if match_REF['SpT_code']<2.6: star.plotspec(4530,4590)
                    else: star.plotspec(6337.11,6357.11)

                    SpT = '-'
                    while SpT not in ['O','B','A']:
                        SpT = input('Choose SpT for the RV0 list of lines (default is B): ')
                        if SpT == '' or SpT == 'B': SpT = 'B'; spt_list = 'rv_Bs.lst'
                        elif SpT == 'A': SpT = 'A'; spt_list = 'rv_As.lst'
                        elif SpT == 'O': SpT = 'O'; spt_list = 'rv_Os.lst'
                    fun = '-'
                    while fun not in ['g','l','v','r','vr']:
                        fun = input('Choose function to fit between g/l/v/r/vr (default is g): ')
                        if fun == '': fun = 'g'
                    wid = '-'
                    while type(wid) is not float:
                        wid = input('Choose the initial width in angstroms (default is 15): ')
                        if wid == '': wid = 15.
                        else: wid = float(wid)

                    plt.close()

                    star.offset = RV0(spt_list,star.spectrum,func=fun,ewcut=30,tol=RV0tol)
                    star.waveflux() # Applies the offset
                    #star.cosmic(sigclip=0.005)

                    if match_REF['SpT_code'] <= 2.5: star.plotspec(4530,4590,poslines='OB')
                    else: star.plotspec(6337.11,6357.11,poslines='OB')

                    input(); plt.close('all')

                    next = input('Type "n" to repeat, hit return to move to the next star. ')

            star.export(tail='_RV',extension='.ascii')

        # MAUI input last modifications:
        if star.resolution == 67000: star.resolution = 85000

        if match_IB['vsini'] < 5:
            match_IB['vsini'] = 0
            match_IB['evsini'] = 0

        if match_IB['rt_mac'][0] < 10:
            match_IB['rt_mac'] = 0
            match_IB['emac'] = 0

        if match_REF['SNR_B'][0] > 200:
            match_REF['SNR_B'] = 200

        maui_txt.write(
            "{:<40}".format(star.file_name[:-5])+\
            "{:<48}".format(star.file_name[:-5]+'_RV.ascii')+"{:<6}".format('0.0d0')+\
            "{:<8}".format(str(int(round(match_IB['vsini'][0],0))))+\
            "{:<7}".format(str(int(round(match_IB['evsini'][0]))))+\
            "{:<8}".format(str(int(round(match_IB['rt_mac'][0]))))+\
            "{:<7}".format(str(int(round(match_IB['emac'][0]))))+\
            "{:<10}".format(str(star.resolution)+'.')+\
            "{:<5}".format(str(int(round(match_REF['SNR_B'][0],0))))+';# '+\
            "{:<12}".format(row['SpC'].strip().replace(' ',''))+' '+\
            "{:<6}".format(str(round(match_REF['FW34Hb'][0]-match_REF['FW14Hb'][0],2)))+\
            "{:<5}".format(str(round(loggf,2)))+"{:<5}".format(str(round(logTf,2)))+\
            str(SiIIIFG)+' '+str(SiIIFG)+\
            '\n')

    maui_txt.close()

    if ascii == 'y': print('Remember to move the new ascii into "ASCII_ARCHIVE" folder')


def shdr(gentable='y'):

    plt.close('all'); plt.figure(figsize=(5,8))

    #mass_list = [.8,.9,1.0,1.1,1.2,1.3,1.5,1.7,2,2.5,3,4,5,7,9,12,15,20,25,32,40,60,85,120]
    #mass_list = [.8,.9,1.,1.1,1.2,1.3,1.5,1.7,2,2.5,3,4,5,7,9,12,15,20,25,32,40,60,85]
    mass_list = [7,9,12,15,20,25,32,40,60,85]

    # Information of each box:
    names = ['BSgs_CNOSiMg_old','BDws_CNOSIMg_old','OBSgs_hot_NOSi_new','BSgs_cool_NOSi_new']
    # nlte_10.1.6_SOLAR_expoclump_2019-10-24.idl
    box1 = [[4.190,4.477,4.477,4.190,4.190],[3.785,3.785,4.391,4.391,3.785]]
    # nlte_10.1.6_bdwarfs_SOLAR_2020-01-29.idl
    box2 = [[4.290,4.543,4.543,4.290,4.290],[2.391,2.391,3.889,3.889,2.391]]
    # nlte_10.4.7_OB.Sg_SOLAR_2021-01-23.idl
    box3 = [[4.399,4.544,4.544,4.399,4.399],[3.488,3.488,4.386,4.386,3.488]]
    # nlte_10.4.7_late.bsgs_SOLAR_expoclump_NOSi.djl_2021-02-06.idl
    box4 = [[4.146,4.322,4.322,4.146,4.146],[3.092,3.092,4.391,4.391,3.092]]

    plt.plot(box1[0],box1[1],lw=1); plt.plot(box2[0],box2[1],lw=1)
    plt.plot(box3[0],box3[1],lw=1); plt.plot(box4[0],box4[1],lw=1)

    # Plot the tracks first:
    for i in mass_list:
        mist = trackmist(mass=i,vr=0.0)
        mist = mist[mist['phase']<=4]
        log_Teff = mist['log_Teff']
        log_LLsol = 4*mist['log_Teff']-mist['log_g']-10.61
        plt.scatter(log_Teff,log_LLsol,s=.3,c=mist['surface_he4'],cmap='gnuplot')
        plt.text(log_Teff[0]+.05,log_LLsol[0]-.07,str(i),fontsize=7)

    # Plot the results from MAUI creating the sub-tables:
    log_Teff = np.asarray(4+np.log10(results['Teff']))
    log_LLsol = np.asarray(5.39-results['lgf'])

    points = np.column_stack([log_Teff,log_LLsol])
    for name,box in zip(names,[box1,box2,box3,box4]):
        verts = np.array([box[0],box[1]]).T
        path = mpath.Path(verts)
        inout = path.contains_points(points)
        log_Teff_in,log_LLsol_in = points[inout].T

        results_in = results[path.contains_points(points)]['Name']
        table_red = table[[i['Name'].strip() in results_in for i in table]]

        if gentable=='y': input_maui(table=table_red,output_name=name,ascii_0='n')

        plt.scatter(log_Teff_in,log_LLsol_in,s=6,label=name)

    #plt.colorbar(shrink=0.75)

    plt.tick_params(direction='in',top='on')
    plt.gca().invert_xaxis()
    plt.xlabel(r"log(T$_{eff})\,[K]$",size=13)
    plt.ylabel(r"log($\mathcal{L}$/$\mathcal{L}_{\odot}$)",size=13)
    plt.xlim(4.8,4)
    plt.ylim(2,5)
    plt.tight_layout()

    plt.legend(); plt.show(block=False)
