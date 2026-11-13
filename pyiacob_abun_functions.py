"""
@author: A. Casasbuenas Corral

GOAL:
    Quantitative spectroscopic analysis of B-type stars, mainly within the IACOB database. 
"""

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import numpy.ma as ma

from astropy.io.ascii import read
from math import ceil
from scipy.interpolate import splrep, BSpline
from scipy.ndimage import gaussian_filter1d,label, binary_dilation #  median_filter, 
import scipy.io 

import glob as glob
import time 
from astropy.table import unique
import os

# FUNCIONES DE ABEL DE BURGOS 
from db import * 
from measure import * 
from spec import *

import warnings
warnings.filterwarnings("ignore")

maindir = str(maindir)+ "/"
fastwinddir = maindir + "FASTWIND/"

tex=False
font="serif"
dpi=180
plt.rcParams.update({
    "font.size": 14,
    "font.family": font,
    "text.usetex": tex,
    "figure.subplot.top": 0.9,
    "figure.subplot.right": 0.9,
    "figure.subplot.left": 0.15,
    "figure.subplot.bottom": 0.12,
    "figure.subplot.hspace": 0.4,
    "savefig.dpi": dpi,
    "savefig.format": "png",
    "axes.titlesize": 16,
    "axes.labelsize": 14,
    "axes.axisbelow": True,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 12,
    "xtick.minor.size": 6,
    "xtick.major.pad": 12,
    "xtick.minor.pad": 6,
    "ytick.major.pad": 12,
    "ytick.minor.pad": 6,
    "ytick.major.size": 12,
    "ytick.minor.size": 6,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 14,
    "legend.framealpha": 1,
    "figure.titlesize": 18,
    "lines.linewidth": 2,
        })
np.set_printoptions(threshold=np.inf)




def FASTWIND_parameters(path,grid_name,savepath,path_evolutionarytracks,
                        el_input='Si',T='*',g='*',He='*',Q='*',b='*'):
    '''
    Function that 
    '''
    
    if el_input=='O':
        abun_centr='8.75'
    elif el_input=='N':
        abun_centr='7.80'
    elif el_input=='Mg':
        abun_centr='7.55'
    elif el_input=='H_He':
        abun_centr='7.80'
    elif el_input=='C':
        abun_centr='8.35'
    elif el_input=='Si':
        abun_centr='7.50'
    else:
        print('Please choose a valid element between O, N, Mg, H_He, C, Si')
    
    # if lines_file != None:
    #     wavelist=findlines(lines_file)[0]
    #     linelambda=findlines(lines_file)[1]
        
    folders=sorted(set(glob.glob(path+'T*g*'+'He'+He+'Q'+Q+'b'+b+'*')))
    folders=[f[:-4] for f in folders if f.endswith('.xdr')]  
    
    T=np.empty(len(folders),dtype=int)
    g=np.empty(len(folders),dtype=float)
    He=np.empty(len(folders),dtype=int)
    Q=np.empty(len(folders),dtype=float)
    b=np.empty(len(folders),dtype=int)
    ab=np.empty(len(folders),dtype=float)
    
    for i,source in zip(range(len(folders)),folders):
        print(source)
        folder=source.split(os.sep)[-1]
        T[i]=int(folder.split('T')[1].split('g')[0])*100
        g[i]=float(folder.split('g')[1].split('He')[0])/100
        He[i]=int(folder.split('He')[1].split('Q')[0] )
        Q[i]=-1*float(folder.split('Q')[1].split('b')[0])/100
        b[i]=int(folder.split('b')[1][:2])
        if folder[-4]=='p':
            abun=(float(abun_centr)+float(folder[-3:])/100)
            formatter = '{:.' + '{}'.format(2) + 'f}'
            ab[i]=float(formatter.format(round(abun,2)))
        elif folder[-4]=='m':
            abun=(float(abun_centr)-float(folder[-3:])/100)
            formatter = '{:.' + '{}'.format(2) + 'f}'
            ab[i]=float(formatter.format(round(abun,2)))
        if i<1:
            Proof=scipy.io.readsav(source+".xdr")
            Models=Proof['l']
            lines=sorted(set([Models[jj][0].decode("utf-8") for jj in range(len(Models)) ]))
            micros=sorted(set([Models[jj][1].decode("utf-8") for jj in range(len(Models)) ] ))
    
    print('#################################################')
    print('Teff 10^3:', sorted(set(np.array(T/1000,dtype=int))))
    print('logg:', sorted(set(g)))    
    print('Q:', sorted(set(Q)))
    print('He:', sorted(set(He)))
    print('b:', sorted(set(b)))
    print('ab:', sorted(set(ab)))
    print("micros: ", sorted(set(micros)))
    print("lines: ", sorted(set(lines)))
    print('#################################################')
    
    Lsp=4*np.log10(T)-g-10.61

    #FICHEROS TRAZAS EVOLUTIVAS
    
    M=[7,9,12,15,20,25,32,40,60,85,120]
    archivos_v0=['M007Z14V0.dat', 'M009Z14V0.dat','M012Z14V0.dat',
                 'M015Z14V0.dat','M020Z14V0.dat','M025Z14V0.dat','M032Z14V0.dat',
                 'M040Z14V0.dat','M060Z14V0.dat','M085Z14V0.dat','M120Z14V0.dat']

    tablas_v0=[]

    G=6.674e-11 #m^3kg^-1s^-2
    Rsun=6.987e8 #m
    Msun=1.988e30 #kg
    gsun=(G*Msun/(Rsun)**2)*10**2 #cm s^-2
    Tsun=5777 #K
    
    for i in range(len(archivos_v0)):
        fichero=path_evolutionarytracks+archivos_v0[i] 
        header=read(fichero, data_start=0, data_end=1)
        list=header[0]
        header=[]
        for ii in range(len(list)):
            header.append(list[ii])
        header[4]='lg(Teff)1'
        data=read(fichero, data_start=2, names=header)
        data['R']=(10**(data['lg(L)']))**(1/2)*((10**(data['lg(Teff)']))/Tsun)**(-2) #Rsun
        data['lgg']=np.log10((data['mass']/data['R']**2)*gsun)
        data['lg(Lsp)']=4*data['lg(Teff)']-data['lgg']-10.61   
        tablas_v0.append(data)

    eps=1e-3
    logTeff_ZAMS_v0=[]
    logL_ZAMS_v0=[]
    logLsp_ZAMS_v0=[]
    logg_ZAMS_v0=[]

    for i in range(len(archivos_v0)):
        logTeff_ZAMS_v0.append((tablas_v0[i])['lg(Teff)'][0])
        logL_ZAMS_v0.append((tablas_v0[i])['lg(L)'][0])
        logLsp_ZAMS_v0.append((tablas_v0[i])['lg(Lsp)'][0])
        logg_ZAMS_v0.append((tablas_v0[i])['lgg'][0])        
        
        
    fig=plt.figure(figsize=(10,14),constrained_layout=True)
    gs = fig.add_gridspec(2,1)
    
    #SPECTROSCOPIC HR DIAGRAM
    ax=fig.add_subplot(gs[0,0])
    for i in range(len(archivos_v0)):
        if M[i]<30:
            ax.plot((tablas_v0[i])['lg(Teff)'],(tablas_v0[i])['lg(Lsp)'],'--r',alpha=0.5)
        else:
            ax.plot((tablas_v0[i])['lg(Teff)'][tablas_v0[i]['1H_cen']>eps],
                    (tablas_v0[i])['lg(Lsp)'][tablas_v0[i]['1H_cen']>eps],'--r',alpha=0.5)
        ax.text(logTeff_ZAMS_v0[i]+0.05,logLsp_ZAMS_v0[i],M[i],zorder=200)

    ax.plot(logTeff_ZAMS_v0,logLsp_ZAMS_v0,'-k',alpha=0.5)
    ax.plot(np.log10(T),Lsp,'g*',label=grid_name)

    #ax.plot([np.log10(33500),np.log10(29500)],[3.35,4.3], '-', c='k', lw=2, zorder=99,alpha=0.5)

    ax.set_ylabel('log $\mathcal{L}/\mathcal{L}_\odot$')
    ax.set_xlabel('$T_{\mathrm{eff}}$ [$10^3$K]')

    ax.set_xlim(4.8,4.1)
    ax.set_ylim(2.2,4.7)

    listofxticks=np.array(sorted(np.arange(4.1,4.8,0.1),reverse=True))
    labelofxticks=['%d'%i for i in 10**(listofxticks)/1000]
    ax.set_xticks(listofxticks)
    ax.set_xticklabels(labelofxticks)
    
    ax.legend(fancybox=True, loc='lower left',fontsize=14)

    ax=fig.add_subplot(gs[1,0])

    for i in range(len(archivos_v0)):
        if M[i]<30:
            ax.plot(((tablas_v0[i])['lg(Teff)']),(tablas_v0[i])['lgg'],'--r')
        else:
            ax.plot(((tablas_v0[i])['lg(Teff)'])[tablas_v0[i]['1H_cen']>eps],
                     (tablas_v0[i])['lgg'][tablas_v0[i]['1H_cen']>eps],'--r')

    ax.plot(np.log10(T),g,'g*',label=grid_name)

    ax.set_xlim(4.8,4.1)
    ax.set_ylim(4.5,1.5)
    listofxticks=np.array(sorted(np.arange(4.1,4.8,0.1),reverse=True))
    labelofxticks=['%d'%i for i in 10**(listofxticks)/1000]
    ax.set_xticks(listofxticks)
    ax.set_xticklabels(labelofxticks)
    ax.set_ylabel('log $g$')
    ax.set_xlabel('$T_{\mathrm{eff}}$ [$10^3$K]')
    
    ax.legend(fancybox=True, framealpha=0.5,loc='upper left',fontsize=14) 
    
    fig.savefig(savepath+grid_name+'_parameters.png')
    plt.show(block=False)

    return T, g, lines

def where_in_grid(path,grid_name,path_evolutionarytracks,
                  T_fixed,g_fixed,el_input='Si'):

    T='*'
    g='*'
    He='*'
    Q='*'
    b='*'
    
    if el_input=='O':
        abun_centr='8.75'
    elif el_input=='N':
        abun_centr='7.80'
    elif el_input=='Mg':
        abun_centr='7.55'
    elif el_input=='H_He':
        abun_centr='7.80'
    elif el_input=='C':
        abun_centr='8.35'
    elif el_input=='Si':
        abun_centr='7.50'
    else:
        print('Please choose a valid element between O, N, Mg, H_He, C, Si')
    
    # if lines_file != None:
    #     wavelist=findlines(lines_file)[0]
    #     linelambda=findlines(lines_file)[1]
        
    folders=sorted(set(glob.glob(path+'T'+T+'g'+g+'He'+He+'Q'+Q+'b'+b+'*')))
    folders=[f[:-4] for f in folders if f.endswith('.xdr')]  
    
    T=np.empty(len(folders),dtype=int)
    g=np.empty(len(folders),dtype=float)
    He=np.empty(len(folders),dtype=int)
    Q=np.empty(len(folders),dtype=float)
    b=np.empty(len(folders),dtype=int)
    # ab=np.empty(len(folders),dtype=float)
    
    
    for i,source in zip(range(len(folders)),folders):
        folder=source.split(os.sep)[-1]
        T[i]=int(folder.split('T')[1].split('g')[0])*100
        g[i]=float(folder.split('g')[1].split('He')[0])/100
        He[i]=int(folder.split('He')[1].split('Q')[0] )
        Q[i]=-1*float(folder.split('Q')[1].split('b')[0])/100
        b[i]=int(folder.split('b')[1][:2])
    #     if folder[-4]=='p':
    #         abun=(float(abun_centr)+float(folder[-3:])/100)
    #         formatter = '{:.' + '{}'.format(2) + 'f}'
    #         ab[i]=float(formatter.format(round(abun,2)))
    #     elif folder[-4]=='m':
    #         abun=(float(abun_centr)-float(folder[-3:])/100)
    #         formatter = '{:.' + '{}'.format(2) + 'f}'
    #         ab[i]=float(formatter.format(round(abun,2)))
    
    # print('#################################################')
    # print('Teff 10^3:', sorted(set(np.array(T/1000,dtype=int))))
    # print('logg:', sorted(set(g)))    
    # print('Q:', sorted(set(Q)))
    # print('He:', sorted(set(He)))
    # print('b:', sorted(set(b)))
    # print('ab:', sorted(set(ab)))
    # print('#################################################')
    
    Lsp=4*np.log10(T)-g-10.61
    print(T_fixed, type(T_fixed))
    Lsp_fixed=4*np.log10(T_fixed)-g_fixed-10.61

    #FICHEROS TRAZAS EVOLUTIVAS
    
    M=[7,9,12,15,20,25,32,40,60,85,120]
    archivos_v0=['M007Z14V0.dat', 'M009Z14V0.dat','M012Z14V0.dat',
                 'M015Z14V0.dat','M020Z14V0.dat','M025Z14V0.dat','M032Z14V0.dat',
                 'M040Z14V0.dat','M060Z14V0.dat','M085Z14V0.dat','M120Z14V0.dat']

    tablas_v0=[]

    G=6.674e-11 #m^3kg^-1s^-2
    Rsun=6.987e8 #m
    Msun=1.988e30 #kg
    gsun=(G*Msun/(Rsun)**2)*10**2 #cm s^-2
    Tsun=5777 #K
    
    for i in range(len(archivos_v0)):
        fichero=path_evolutionarytracks+archivos_v0[i] 
        header=read(fichero, data_start=0, data_end=1)
        list=header[0]
        header=[]
        for ii in range(len(list)):
            header.append(list[ii])
        header[4]='lg(Teff)1'
        data=read(fichero, data_start=2, names=header)
        data['R']=(10**(data['lg(L)']))**(1/2)*((10**(data['lg(Teff)']))/Tsun)**(-2) #Rsun
        data['lgg']=np.log10((data['mass']/data['R']**2)*gsun)
        data['lg(Lsp)']=4*data['lg(Teff)']-data['lgg']-10.61   
        tablas_v0.append(data)

    eps=1e-3
    logTeff_ZAMS_v0=[]
    logL_ZAMS_v0=[]
    logLsp_ZAMS_v0=[]
    logg_ZAMS_v0=[]

    for i in range(len(archivos_v0)):
        logTeff_ZAMS_v0.append((tablas_v0[i])['lg(Teff)'][0])
        logL_ZAMS_v0.append((tablas_v0[i])['lg(L)'][0])
        logLsp_ZAMS_v0.append((tablas_v0[i])['lg(Lsp)'][0])
        logg_ZAMS_v0.append((tablas_v0[i])['lgg'][0])        
        
        
    fig=plt.figure(figsize=(10,14),constrained_layout=True)
    gs = fig.add_gridspec(2,1)
    
    #SPECTROSCOPIC HR DIAGRAM
    ax=fig.add_subplot(gs[0,0])
    for i in range(len(archivos_v0)):
        if M[i]<30:
            ax.plot((tablas_v0[i])['lg(Teff)'],(tablas_v0[i])['lg(Lsp)'],'--r',alpha=0.5)
        else:
            ax.plot((tablas_v0[i])['lg(Teff)'][tablas_v0[i]['1H_cen']>eps],
                    (tablas_v0[i])['lg(Lsp)'][tablas_v0[i]['1H_cen']>eps],'--r',alpha=0.5)
        ax.text(logTeff_ZAMS_v0[i]+0.05,logLsp_ZAMS_v0[i],M[i],zorder=200)

    ax.plot(logTeff_ZAMS_v0,logLsp_ZAMS_v0,'-k',alpha=0.5)
        
    ax.plot(np.log10(T),Lsp,'g*',label=grid_name)
    ax.scatter(np.log10(T_fixed),Lsp_fixed,s=100,marker='o',color='red')

    ax.set_ylabel('log $\mathcal{L}/\mathcal{L}_\odot$')
    ax.set_xlabel('$T_{\mathrm{eff}}$ [$10^3$K]')

    ax.set_xlim(4.8,4.1)
    ax.set_ylim(2.2,4.7)

    listofxticks=np.array(sorted(np.arange(4.1,4.8,0.1),reverse=True))
    labelofxticks=['%d'%i for i in 10**(listofxticks)/1000]
    ax.set_xticks(listofxticks)
    ax.set_xticklabels(labelofxticks)
    

    ax.legend(fancybox=True, loc='lower left',fontsize=14)

    ax=fig.add_subplot(gs[1,0])

    for i in range(len(archivos_v0)):
        if M[i]<30:
            ax.plot(((tablas_v0[i])['lg(Teff)']),(tablas_v0[i])['lgg'],'--r')
        else:
            ax.plot(((tablas_v0[i])['lg(Teff)'])[tablas_v0[i]['1H_cen']>eps],
                     (tablas_v0[i])['lgg'][tablas_v0[i]['1H_cen']>eps],'--r')

    ax.plot(np.log10(T),g,'g*',label=grid_name)
    ax.scatter(np.log10(T_fixed),g_fixed,s=100,marker='o',color='red')

    ax.set_xlim(4.8,4.1)
    ax.set_ylim(4.5,1.5)
    listofxticks=np.array(sorted(np.arange(4.1,4.8,0.1),reverse=True))
    labelofxticks=['%d'%i for i in 10**(listofxticks)/1000]
    ax.set_xticks(listofxticks)
    ax.set_xticklabels(labelofxticks)
    ax.set_ylabel('log $g$')
    ax.set_xlabel('$T_{\mathrm{eff}}$ [$10^3$K]')
    
    ax.legend(fancybox=True, framealpha=0.5,loc='upper left',fontsize=14) 
    
    #fig.savefig(savepath+grid_name+'_parameters.png')
    plt.show(block=False)
    
    return None

def table_db_FASTWIND_EW(path,grid_name,el_input='Si',T='*',g='*',He='*',Q='*',b='*',join_table=None):

    '''
    Function that generates a FITS table with the EW information of the lines
    of a certain element for a path with FASTWIND models with a certain set of
    parameters T, g, He, b, in order to calculate the element abundances with
    the curve of growth method
    
    Parameters
    ----------
    path : str
        path of the folder with the FASTWIND models. Inside, there should be folders with
        a name structure similar to T310g400He10Q1350b10CNOp000 and in each of these folders,
        files with a name structure similar to OUT.C1173bl_VT030. 
        E.g. maindir+'FASTWIND/iac_bmsgrid_HHeSi_Z100_He10_p1/'

    grid_name: str
        Name of the grid, for reference. 
        E.g. iac_bmsgrid_HHeSi_Z100
    
    el_input : str 
        Element consdiered for the abundance calculation 
        E.g. write "N" for considering only the Nitrogen lines.
        The valid elements so far are O, N, Mg, H_He, C, Si
    
    T : str 
        Effective temperature of the model to be taken into account.
        E.g. write "310" if you want the models with T=31000K
        If '*', it will look for all available T inside the folder
        Default is '*'
        
    g : str 
        Surface gravity of the model to be taken into account, 
        E.g. write "400" if you want the models with log g=4.00
        If '*', it will look for all available g inside the folder
        Default is '*'        
        
    He : str 
        Helium abundance of the model to be taken into account, 
        E.g. write "10" if you want the models with ab(He)=1.0
        If '*', it will look for all available He inside the folder
        Default is '*'
        
    Q : str 
        Wind strength to be taken into account.
        E.g. write "1350" if you want the models with log Q=-13.5 dex
        If '*', it will look for all available Q inside the folder
        Default is '*'        
        
    b : str 
        Wind to be taken into account. 
        E.g. write "10" if you want the models with b=10
        If '*', it will look for all available b inside the folder
        Default is '*'        
    
    lines_file: str
        If specified, the program will collect the EW information for the lines
        indicated in that file.
        If None, the program will create a file with the FASTWIND lines available.
        Default is None. 

    join_table: str
        If "none", creates a new file
        If not "none", joins the table to another existing table. 

    Returns
    -------
  
    Nothing, but saves the list of FASTWIND lines in:
        lists/lines/*el_input*_lines_FASTWIND.lst
        
    And also saves the table with the EW information of these lines in:
        tables/measure_FASTWIND_*grid_name*_*el_input*.fits 
    '''
    

    start_time=time.time()
    folders=sorted(set(glob.glob(path+'T'+T+'g'+g+'He'+He+'Q'+Q+'b'+b+'*')))
    
    folders=[f for f in folders if f.endswith('.xdr')]   
    
    bar = pb.ProgressBar(maxval=len(folders),widgets=[pb.Bar('=','[',']'),' ',pb.Percentage()])
    bar.start()       
    tab_index=0
    
    for source,ii in zip(folders,range(len(folders))):
        
        print(ii,'/',len(folders))
        folder=source #.split(os.sep)[-1] #nombre del xdr
        Proof=scipy.io.readsav(folder)
        Param=Proof['m']
        Models=Proof['l']
        lines=[Models[jj][0].decode("utf-8") for jj in range(len(Models)) ]
        micros=[Models[jj][1].decode("utf-8") for jj in range(len(Models)) ] 
        EWs=[Models[jj][2] for jj in range(len(Models)) ]
        WL=[Models[jj][3][80] for jj in range(len(Models)) ]

        # WL = [(Models[jj][3])[np.nanargmin(Models[jj][4])] for jj in range(len(Models)) ]

        if ii<1:

            # CREAR UNA LISTA DE LÍNEAS, SOLO UNA VEZ
            linemask=(np.array(((np.array(WL,dtype=float))>3700)) & np.array(((np.array(WL,dtype=float)))<9000))
            wavelist=np.array(WL)[linemask]
            linelist=np.array(lines)[linemask]
            tab_lines=Table([wavelist,linelist],names=['wave','line'],dtype=['float64','S20'])
            tab_lines=unique(tab_lines, keys='line')
            tab_lines=tab_lines[ [I for I,J in zip(range(len(tab_lines['line'])),tab_lines['line']) if J.startswith(el_input)] ] 
            tab_lines.sort('line')
            
            tab_lines['wave','line'].write(maindir+'lists/lines/'+el_input+'_lines_FASTWIND.lst',
                format='ascii.commented_header',delimiter=',',overwrite=True)    
                                        
            wavelist=tab_lines['wave']
            linelist=tab_lines['line']

            # CREAR LA TABLA, SOLO UNA VEZ 
            
            empty_str=np.empty(len(folders)*len(Models),dtype=str)
            empty_int=np.empty(len(folders)*len(Models),dtype=int)
            empty_float=np.empty(len(folders)*len(Models),dtype=float)
            
            if el_input=='O':
                abunname='abun_O'
            elif el_input=='N':
                abunname='abun_N'
            elif el_input=='Mg':
                abunname='abun_Mg'
            elif el_input=='C':
                abunname='abun_C'
            elif el_input=='Si':
                abunname='abun_Si'
            else:
                print('Please choose a valid element between O, N, Mg, H_He, C, Si')
                    
            tab_final=Table([empty_str,empty_str,empty_str,empty_str,empty_str,empty_str,empty_str,empty_str,empty_str],
                names=['Folder','T','g','He','Q','b',abunname,'micro','Grid_name'],
                dtype=['S200','S10','S10','S10','S10','S10','S10','S10','S50'])
            for linesew in ['EW_'+lineii for lineii in linelist]:
                tab_final[linesew]=empty_float
            
        for microjj in sorted(set(micros)):
            tab_final['Folder'][tab_index]=(source.split(os.sep)[-1])[:-4]
            tab_final['T'][tab_index]='%d'%(Param.teff/100)
            tab_final['g'][tab_index]='%d'%(np.round(Param.logg*100,-1))

            # if len('%d'%(Param.yhe*100)) < 2:
            #     yHe = '0%d'%(Param.yhe*100)
            # else:
            #     yHe = '%d'%(Param.yhe*100) 
            # tab_final['He'][tab_index]= yHe

            tab_final['He'][tab_index]= '%d'%(Param.yhe*100) 
            tab_final['Q'][tab_index]='%d'%(np.round(Param.Q*10)*-10)
            tab_final['b'][tab_index]='%d'%(Param.beta*10)

            if el_input=='O':
                tab_final['abun_O'][tab_index]='%.2f'%(Param.abun_o+12)
            elif el_input=='N':
                tab_final['abun_N'][tab_index]='%.2f'%(Param.abun_n+12)
            elif el_input=='Mg':
                tab_final['abun_Mg'][tab_index]='%.2f'%(Param.abun_mg+12)
            elif el_input=='C':
                tab_final['abun_C'][tab_index]='%.2f'%(Param.abun_c+12)
            elif el_input=='Si':
                tab_final['abun_Si'][tab_index]='%.2f'%(Param.abun_si+12)
            else:
                print('Please choose a valid element between O, N, Mg, H_He, C, Si')

            tab_final['micro'][tab_index]=microjj 
            for linejj in sorted(set(lines)):
                if linejj in linelist:
                    mask=( (np.array(lines)==linejj) & (np.array(micros)==microjj) ) 
                    try:
                    #print(mask)                        
                        lindex=(np.argwhere(mask>0))[0][0]
                        tab_final['EW_'+linejj][tab_index]=Models[lindex][2]*(-1000)
                        if tab_final['EW_'+linejj][tab_index]<0:
                            tab_final['EW_'+linejj][tab_index] = +0.0
                    except:
                        tab_final['EW_'+linejj][tab_index]=-999
                    
            tab_final['Grid_name'][tab_index]=grid_name
            tab_index=tab_index+1

        end_time=time.time()
        elapsed_time=end_time-start_time
        bar.update(ii)

        print(' %d seconds'%(elapsed_time))
    
    tab_final=tab_final[tab_final['micro']!='']
    print('\n')
    bar.finish()  

    tab_final.sort(['Folder','micro','abun_'+el_input])

    if join_table is not None:
        join_tab = findtable(join_table)
        if join_tab is not None:
            tab_final = vstack([join_tab, tab_final])     
        tab_final.write(maindir+'tables/' +join_table,format='fits',overwrite=True)
        print('\nFiles saved as: \n\tlists/lines/'+el_input+'_lines_FASTWIND.lst','\n\ttables/' + join_table )    

    else:
        tab_final.write(maindir+'tables/measure_FASTWIND_'+grid_name+'_'+el_input+'.fits',
                              format='fits',overwrite=True)                     
        print('\nFiles saved as: \n\tlists/lines/'+el_input+'_lines_FASTWIND.lst','\n\ttables/measure_FASTWIND_'+
         	grid_name+'_'+el_input+'.fits' )        
          
    return None


def FASTWIND_tab_interpolate(T_new,grid_table,g,He,Q,b,el_input='Si',n=3):
    '''
    Parameters
    ----------
    T_new : str
        New temperature, in 10^3 K, for which to interpolate the EW for each of the lines. 
        
    grid_table : str
        Name of the table generated with the function table_db_FASTWIND_abundances 
        with the desired lines and EWs for the temperatura range.
        
    n: int
        Number of T point with which the program will try to interpolate.
        Minimum number of points is two. Default is 4. 
        
    el_input:
        Element for which to create the table.
        
        
    Returns
    -------
    None, but saves the table with the EW information of these lines in the
    same table:
        grid_table 

    '''
    
    FW_table_og=findtable(grid_table,path=maindir+'tables/')
    FW_table_og = FW_table_og[(FW_table_og['He']==He)]
    FW_table_og=FW_table_og[FW_table_og['Grid_name']==FW_table_og['Grid_name'][0]]
    FW_table=FW_table_og[(FW_table_og['g']==g) & 
                         (FW_table_og['Q']==Q) & (FW_table_og['b']==b)] 
    list_of_tables=[]
    
    T_list=sorted(set(FW_table['T']))
    
    if len(T_list)<n:
        print('ERROR: Try with different parameters: not enough T points in the table '+
              'to interpolate with ' + str(n) + ' points.')
        return None
    
    Tarr=[]
    for j in range(n):
        nearest=T_list[min(range(len(T_list)),key=lambda i : abs(int(T_list[i])-int(T_new)))]
        if j<1:
            new_tab=FW_table[FW_table['T']==nearest]

        list_of_tables.append(FW_table[FW_table['T']==nearest])
        Tarr.append(nearest)
        T_list.remove(nearest)
    
    Tarr=sorted(Tarr)
    
    new_tab['T']=T_new
    
    bar = pb.ProgressBar(maxval=len(new_tab),widgets=[pb.Bar('=','[',']'),' ',pb.Percentage()])
    bar.start()
    
    new_tab['Folder']=new_tab['Folder'].astype('60str')
    
    if np.sum(np.array(Tarr,dtype=int)<int(T_new))!=np.sum(np.array(Tarr,dtype=int)>int(T_new)):
        print('WARNING: asymetric interpolation')
        new_tab['Grid_name']='asymetric-interp'
    else:
        print('Symetric interpolation')
        new_tab['Grid_name']='symetric-interp'
    
    ind=0

    for ab_i in list(sorted(set(FW_table['abun_'+el_input]))):
        for micro_i in list(sorted(set(FW_table['micro']))):
            mask=((FW_table['micro']==micro_i) & (FW_table['abun_'+el_input]==ab_i))
            masknewtab=((new_tab['micro']==micro_i) & (new_tab['abun_'+el_input]==ab_i)) 
            for colname in new_tab.colnames:
                if 'EW' in colname:
                    EW_arr=np.empty(len(Tarr),dtype=float)
                    for ii in range(len(Tarr)):
                        try:
                            EW_arr[ii]=FW_table[colname][(FW_table['T']==Tarr[ii]) & mask][0]
                            cont='y'
                        except:
                            print('ERROR in', ab_i,micro_i,Tarr[ii], g)
                            EW_arr[ii]=np.nan
                            cont='n'
                    if cont=='y':
                        tck=splrep(np.array(Tarr,dtype=int),np.array(EW_arr,dtype=float),s=50)
                        new_tab[colname][masknewtab]=BSpline(*tck)(np.array(T_new,dtype=int))
                        new_tab['Folder'][masknewtab]='T(int)g'+(FW_table['Folder'][(FW_table['T']==Tarr[ii]) & mask][0]).split('g')[1]
                        
                    else:
                        print('NO MODEL WITH AB=%s, MICRO=%s, T=%s'%(ab_i,micro_i,Tarr[ii]))
            ind=ind+1
            # bar.update(ind)
            
    print('\n')
    bar.finish()
    
    new_tab=vstack([FW_table_og,new_tab])
    
    #PARA QUE VALGA PARA CUALQUIER ELEMENTO, COPIAR EL FINAL Y EL PRINCIPIO DE LOS ARCHIVOS
    #ANTES DE T Y DESPUÉS DE G, Y SOLO CAMBIAR EL TNEW (y asi no tener que especificar elinput)
    
    new_tab.write(maindir+'tables/'+grid_table,format='fits',overwrite=True) 
        
    print('File saved as: tables/'+grid_table)
    
    return None

#%%

def plot_lines(ref_file,spc_star,orig,
               star_measure,list_of_lines,savepath,
               list_of_models=None,broadening=None):
    '''
    Parameters
    ----------
    star : str
        Name of the star or spectra file to 
        
    spc_star: str
    
    orig: str
        Supports 'IACOB' or 'txt'
    
    star_measure_table : TYPE
        DESCRIPTION.
        
    list_of_lines : TYPE
        DESCRIPTION.
        
    savepath: TYPE
        DESCRIPTION
        
    list_of_models : TYPE
        List of the parameters of the models:
            [[grid_path,el_input,T,g,He,Q,b,micro,abun],[grid_path,el_input,T2,g2,He2,Q2,b2,micro2,abun2]]
        
    broadening: TYPE
        DESCRIPTION

    Returns
    -------
    None.

    '''

    star_measure=findtable(star_measure,path=maindir+'tables/')
    
    star_measure=star_measure[star_measure['Ref_file']==ref_file]
    star_name=star_measure['ID'][0]    

    SpC=spc_star
         
    rv0=star_measure['RV0']
    star=spec(ref_file, orig=orig,rv0=rv0)
    wave=star.wave
    flux=star.flux
        
    lines=findlines(list_of_lines)
    linewaves,linenames= lines[0], lines[1]
    
    '''
    ARREGLAR PINTAR UN MODELO ENCIMA
    (los modelos de fastwind son para cada línea...)
    Igual tengo que eliminar la función que busca los modelos y buscarlos aquí.
    + Permitir pintar varios modelos (llamarlos a, b, c.... 
    y luego poner la info de ellos abajo de subtítulo si se puede?)
    '''
    
    if list_of_models is not None:
        grid_path,el_input,T,g,He,Q,b,micro,abun=list_of_models
        vsini,vmac,R=broadening
        model_wave=[]
        model_flux=[]
        for linei in linenames:
            modeli=find_FW_model(grid_path,el_input,T,g,He,Q,b,abun)
            if modeli is not None:
                m=broadened_FASTWIND_model(modeli,micro,linei,int(vsini),int(vmac),int(R))
                if m is not None:
                    wavei,fluxi=m[0],m[1]
                    model_wave.append(list(wavei))
                    model_flux.append(list(fluxi))
                else:
                    model_wave.append(list([0]))
                    model_flux.append(list([0]))
            else:
                print('Model not found: ', list_of_models)
                model_wave.append(list([0]))
                model_flux.append(list([0]))
    rows=ceil(len(linewaves)/4)
    cols=60
    ncols=[4,3,3,4,5]
        
    fig=plt.figure(figsize=(15,20),constrained_layout=True)
    gs = fig.add_gridspec(rows, cols)
    count=0
    for i in range(rows):
        start=0
        step=int(cols/ncols[i])
        for j in range(ncols[i]):
            xrange=(linewaves[count]-20, linewaves[count]+20)
            
            ax=fig.add_subplot(gs[i,start:start+step])
            ax.plot(wave,flux,label=linenames[count],color='navy',alpha=0.75)
            
            if list_of_models is not None:
                ax.plot(model_wave[count],model_flux[count],'--',color='green',label='model',alpha=0.75)
            ax.axvline(linewaves[count], color='red',ls='--', alpha=0.75)
            ax.set_xlim(xrange[0],xrange[1])
            ax.set_ylim(0.4,1.25)
            ax.set_xlabel(r'$\lambda$ ($\AA$)')
            ax.legend(loc='upper right',fancybox=True, framealpha=0.7)
            if start<1:
                ax.set_ylabel('Norm.  flux')
            else:
                ax.set_yticks([])
            count=int(count+1)
            start=start+step
            
    fig.suptitle('%s (%s)'%(ref_file,SpC),fontsize='25') 
    plt.savefig(savepath+ref_file.split(".")[0]+'_Tglines.png') 
    plt.show(block=False) 
    
#%%

def gaussian_reconstruct(RV,FWHM,DEP,lam_rest):
    
    # # --- parameters you already have ---
    # RV = -32.5              # km/s
    # FWHM = 0.35             # Å
    # DEP = 0.42              # depth (relative to continuum)
    # lam_rest = 6562.8       # H-alpha for example, Å

    # --- reconstruct Gaussian ---
    c = 299792.458          # km/s
    RV = 0 #realmente ya está corregido en RV
    lam0 = lam_rest * (1 + RV / c)
    sigma = FWHM / (2*np.sqrt(2*np.log(2)))

    wave = np.linspace(lam0 - 5*FWHM, lam0 + 5*FWHM, 2000)

    flux = 1 - DEP * np.exp(- (wave - lam0)**2 / (2*sigma**2))

    # flux now contains your reconstructed fitted line    

    return wave, flux   


#%%

def plot_best_model_element(ref_file,el_input,star_db,star_measure,
                    linestab,FW_EW_measure,parameterstable,savepath):

    print("Working on plotting the best fitting model to each %s line..."%(el_input)) 
    
    star_dbtab=findtable(star_db) 
    star_dbtab = star_dbtab[star_dbtab['Ref_file']==ref_file]

    star_name=star_dbtab['ID'][0]    
    if "Class." in star_dbtab.colnames:
        spc_star=star_dbtab['Class.'][0]   
    elif "SpType" in findtable(star_db).colnames:
        spc_star=findtable(star_db)['SpType'][findtable(star_db)['Ref_file']==ref_file][0]         
    else:
        spc_star = "-"
    
    if "vsini" in star_dbtab.colnames: 
        vsini_star=star_dbtab['vsini'][0]
    elif "vsini_GF" in star_dbtab.colnames: 
        vsini_star=star_dbtab['vsini_GF'][0]
    if "rtmac" in star_dbtab.colnames: 
        vmac_star=star_dbtab['rtmac'][0]
    elif "vmac" in star_dbtab.colnames: 
        vmac_star=star_dbtab['vmac'][0]
    elif "vmac_GF" in star_dbtab.colnames: 
        vmac_star=star_dbtab['vmac_GF'][0]

    R_star=star_dbtab['Resol'][0]
    star_measuretab=findtable(star_measure)
    FW_EW_measuretab=findtable(FW_EW_measure)
    FW_EW_measuretab = FW_EW_measuretab[(FW_EW_measuretab["Grid_name"] != "symetric-interp") &
                                        (FW_EW_measuretab["Grid_name"] != "asymetric-interp")]
    # FW_EW_measuretab=FW_EW_measuretab[FW_EW_measuretab['Grid_name']==FW_EW_measuretab['Grid_name'][0]]

    star_measuretab=star_measuretab[star_measuretab['Ref_file']==ref_file]
    rv0=star_measuretab['RV0'][0]

    if '.ascii' in ref_file:
        orig='txt'
    elif '.txt' in ref_file:
        orig='txt'
    elif 'fit' in ref_file:
        orig='IACOB'
        
    star=spec(ref_file, snr='bestHF', orig=orig,rv0=rv0)

    wave=star.wave
    flux=star.flux        

    fulllinestab = findtable(linestab,format='csv',path=savepath)
    linewaves = fulllinestab ['WL']
    linenames = fulllinestab['label']
    nn = fulllinestab[ref_file]
    flagtab = fulllinestab['flag']
    
    format_i = "csv" if ".csv" in parameterstable else "ecsv" if ".ecsv" in parameterstable else "ascii"
    delimiter_i = "," if ".csv" in parameterstable else "\s" if ".ecsv" in parameterstable else "\s"    
    parameters=findtable(parameterstable,path=savepath,format=format_i,delimiter=delimiter_i)

    parameters = parameters[parameters["Ref_file"]==ref_file]

    flag=parameters['FLAG'][0]

    fixed_T=parameters['T'][0]
    fixed_g=parameters['g'][0]
    fixed_He=parameters['He'][0]
    fixed_Q=parameters['Q'][0]
    fixed_b=parameters['b'][0]
    abun_final=parameters['ab'+el_input][0]
    micro_center=parameters['micro'][0]
    microlist=sorted(set(FW_EW_measuretab['micro']))
    micro_final=microlist[min(range(len(microlist)),key=lambda i : abs(int(microlist[i])-float(micro_center)))]

    FW_EW_measuretab = FW_EW_measuretab[(FW_EW_measuretab["Grid_name"] != "symetric-interp") &
                                        (FW_EW_measuretab["Grid_name"] != "asymetric-interp")]
    # FW_EW_measuretab=FW_EW_measuretab[FW_EW_measuretab['Grid_name']==FW_EW_measuretab['Grid_name'][0]]
    
    gFW=sorted(set(FW_EW_measuretab['g']))
    nearest_g=gFW[min(range(len(gFW)),key=lambda i : abs(int(gFW[i])-int(fixed_g)))] #NEAREST
    
    TFW=sorted(set(FW_EW_measuretab['T'][FW_EW_measuretab['g']==str(int(nearest_g))]))
    nearest_T=TFW[min(range(len(TFW)),key=lambda i : abs(int(TFW[i])-int(fixed_T/100)))] #NEAREST
    
    abunFW=sorted(set(FW_EW_measuretab['abun_'+el_input][(FW_EW_measuretab['g']==str(int(nearest_g))) & 
                                                         (FW_EW_measuretab['T']==str(int(nearest_T)))]))   
    nearest_abun=abunFW[min(range(len(abunFW)),key=lambda i : abs(float(abunFW[i])-float(abun_final)))] #NEAREST
    
    modeli=find_FW_model(FW_EW_measuretab['Grid_name'][FW_EW_measuretab["He"] == str(fixed_He)][0], 
                         el_input, str('%d'%(int(nearest_T))), str(int(nearest_g)),
                  str(fixed_He), str(fixed_Q), str(fixed_b),  str(nearest_abun))

    nrows, ncols = even_plot(len(linenames))
    fig,ax = plt.subplots(nrows, ncols, figsize=(16,16), layout = "constrained")
    axs =  ax.flatten()

    start=0
    
    fluxmin=1
    fluxmax=1.3

    for i in range(len(linenames)):

        wavei=linewaves[i]
        linei=linenames[i]   
        flagi = int(nn[i])

        if wavei >=np.amax(wave)-5:
            continue
        elif wavei <= np.amin(wave)+5:
            continue

        xrange=(wavei-5, wavei+5)
        wave_st=wave[(wave>xrange[0]) & (wave<xrange[1])]
        flux_st=flux[(wave>xrange[0]) & (wave<xrange[1])]

        wavemodi,fluxmodi=broadened_FASTWIND_model(modeli,micro_final, linei ,int(vsini_star),
                                        int(vmac_star),R=int(R_star)) 
                
        axs[i].axvline(wavei, color='tab:blue',ls='--', alpha=0.5)
        axs[i].plot(wave_st,flux_st,label=linei,c='navy')

        axs[i].set_xlim(np.amin(wave_st),np.amax(wave_st))
        axs[i].set_xlabel(r'$\lambda$ ($\AA$)')

        if flagi <0.5:
            axs[i].set_title(linei, color = "tab:green")
            axs[i].plot(wavemodi,fluxmodi,alpha=0.75,c='tab:green')
        
        elif flagi <1.5: 
            axs[i].plot(wavemodi,fluxmodi,alpha=0.75,c='grey')
            axs[i].set_title(linei, color = "grey")
        
        else: 
            axs[i].plot(wavemodi,fluxmodi,alpha=0.75,c='darkred')
            axs[i].set_title(linei, color = "darkred")

    fig.suptitle(r"%s (%s) log g=%.1f dex, Teff=%i K, log Q=-%.1f dex, $\beta$=%.1f, $\xi_t$=%d km/s, $\epsilon$(%s)=%.2f dex" 
                  %(star_name,spc_star,int(nearest_g)/100,int(nearest_T)*100, 
                  float(fixed_Q)/100,float(fixed_b)/100,int(micro_final),el_input,float(nearest_abun)),
                  fontsize=25)
    
    print("Saving the image...")
    fig.savefig(savepath+ref_file.split(".")[0]+'_bestfit_%s.png'%el_input, bbox_inches = "tight")    

    stop = input("Press return to finish with the analysis.")

    plt.show(block=False)

    print("DONE!")

    return wave, flux
    

def plot_best_model(ref_file,el_input,star_db,star_measure,
                    linelist,FW_EW_measure,parameterstable,savepath):

    ######################################################################
    ############# PLOT DE LAS LÍNEAS DE SILICIO! #########################
    ######################################################################
    star_dbtab=findtable(star_db) 
    star_dbtab = star_dbtab[star_dbtab['Ref_file']==ref_file]

    star_name=star_dbtab['ID'][0]    
    if "Class." in star_dbtab.colnames:
        spc_star=star_dbtab['Class.'][0]    
    elif "SpType" in findtable(star_db).colnames:
        spc_star=findtable(star_db)['SpType'][findtable(star_db)['Ref_file']==ref_file][0]        
    else:
        spc_star = "-"

    if "vsini" in star_dbtab.colnames: 
        vsini_star=star_dbtab['vsini'][0]
    elif "vsini_GF" in star_dbtab.colnames: 
        vsini_star=star_dbtab['vsini_GF'][0]
    if "rtmac" in star_dbtab.colnames: 
        vmac_star=star_dbtab['rtmac'][0]
    elif "vmac" in star_dbtab.colnames: 
        vmac_star=star_dbtab['vmac'][0]
    elif "vmac_GF" in star_dbtab.colnames: 
        vmac_star=star_dbtab['vmac_GF'][0]

    R_star=star_dbtab['Resol'][0]
    star_measuretab=findtable(star_measure)
    FW_EW_measuretab=findtable(FW_EW_measure)
    FW_EW_measuretab = FW_EW_measuretab[(FW_EW_measuretab["Grid_name"] != "symetric-interp") &
                                        (FW_EW_measuretab["Grid_name"] != "asymetric-interp")]    
    # FW_EW_measuretab=FW_EW_measuretab[FW_EW_measuretab['Grid_name']==FW_EW_measuretab['Grid_name'][0]]

    star_measuretab=star_measuretab[star_measuretab['Ref_file']==ref_file]
    rv0=star_measuretab['RV0'][0]

    if '.ascii' in ref_file:
        orig='txt'
    elif '.txt' in ref_file:
        orig='txt'
    elif '.dat' in ref_file:
        orig = "txt"
    elif 'fit' in ref_file:
        orig='IACOB'
        
    star=spec(ref_file, snr='bestHF', orig=orig,rv0=rv0)

    wave=star.wave
    flux=star.flux        

    linenames=findlines(linelist)[1]
    linewaves=findlines(linelist)[0]
    
    format_i = "csv" if ".csv" in parameterstable else "ecsv" if ".ecsv" in parameterstable else "ascii"
    delimiter_i = "," if ".csv" in parameterstable else "\s" if ".ecsv" in parameterstable else "\s"    
    parameters=findtable(parameterstable,path=savepath,format=format_i,delimiter=delimiter_i)

    parameters = parameters[parameters["Ref_file"]==ref_file]

    flag=parameters['FLAG'][0]

    fixed_T=parameters['T'][0]
    fixed_g=parameters['g'][0]
    fixed_He=parameters['He'][0]
    fixed_Q=parameters['Q'][0]
    fixed_b=parameters['b'][0]
    abun_final=parameters['ab'+el_input][0]
    micro_center=parameters['micro'][0]
    microlist=sorted(set(FW_EW_measuretab['micro']))
    micro_final=microlist[min(range(len(microlist)),key=lambda i : abs(int(microlist[i])-float(micro_center)))]

    # FW_EW_measuretab=FW_EW_measuretab[FW_EW_measuretab['Grid_name']==FW_EW_measuretab['Grid_name'][0]]

    FW_EW_measuretab = FW_EW_measuretab[(FW_EW_measuretab["Grid_name"] != "symetric-interp") &
                                        (FW_EW_measuretab["Grid_name"] != "asymetric-interp")]    
    
    gFW=sorted(set(FW_EW_measuretab['g']))
    nearest_g=gFW[min(range(len(gFW)),key=lambda i : abs(int(gFW[i])-int(fixed_g)))] #NEAREST
    
    TFW=sorted(set(FW_EW_measuretab['T'][FW_EW_measuretab['g']==str(int(nearest_g))]))
    nearest_T=TFW[min(range(len(TFW)),key=lambda i : abs(int(TFW[i])-int(fixed_T/100)))] #NEAREST
    
    abunFW=sorted(set(FW_EW_measuretab['abun_'+el_input][(FW_EW_measuretab['g']==str(int(nearest_g))) & 
                                                         (FW_EW_measuretab['T']==str(int(nearest_T)))]))   
    nearest_abun=abunFW[min(range(len(abunFW)),key=lambda i : abs(float(abunFW[i])-float(abun_final)))] #NEAREST
      
    #H
    Balmer_lineindex=[index for index, i in enumerate(linenames) if (i[0]=='H') and (i[1:3]!='EI')]
    #HE
    He_lineindex=[index for index, i in enumerate(linenames) if (i[0]=='H') and (i[1:3]=='EI')]
    #SI/el_input
    El_lineindex=[index for index, i in enumerate(linenames) if (i[0:len(el_input)]==el_input)]
    
    rows0=3
    cols0=len(He_lineindex)*len(Balmer_lineindex)*len(El_lineindex)
    
    fig0=plt.figure(figsize=(20,10),constrained_layout=True)    
    
    
    modeli=find_FW_model(FW_EW_measuretab['Grid_name'][FW_EW_measuretab["He"] == str(fixed_He)][0], 
                         el_input, str('%d'%(int(nearest_T))), str(int(nearest_g)),
                  str(int(fixed_He)), str(int(fixed_Q)), str(int(fixed_b)),  str("%.2f"%float(nearest_abun)))
    
    print(modeli)

    print(FW_EW_measuretab['Grid_name'][FW_EW_measuretab["He"] == str(fixed_He)][0], 
          el_input, str('%d'%(int(nearest_T))), str(int(nearest_g)),
                  str(fixed_He), str(fixed_Q), str(fixed_b),  str(nearest_abun), el_input)
    
    print("Working on plotting the spectrum vs the best fitting model...")
    
    #FILA 1: BALMER

    print("H LINES")
    
    row0=0
    gs0=fig0.add_gridspec(rows0,cols0)

    start=0
    step=int(cols0/len(Balmer_lineindex))
    fluxmin=1
    fluxmax=1.3
    for c,ind in zip(range(len(Balmer_lineindex)),Balmer_lineindex):
        wavei=linewaves[ind]
        xrange=(wavei-15, wavei+15)
        flux_st=flux[(wave>xrange[0]) & (wave<xrange[1])]
        if np.nanmin(flux_st)<fluxmin:
            fluxmin=np.nanmin(flux_st)-0.05
        if np.nanmax(flux_st)>fluxmax:
            fluxmax=np.nanmax(flux_st)+0.05
            
    for c,ind in zip(range(len(Balmer_lineindex)),Balmer_lineindex):
        wavei=linewaves[ind]
        linei=linenames[ind]
        
        xrange=(wavei-15, wavei+15)
        wave_st=wave[(wave>xrange[0]) & (wave<xrange[1])]
        flux_st=flux[(wave>xrange[0]) & (wave<xrange[1])]

        wavemodi,fluxmodi=broadened_FASTWIND_model(modeli,micro_final, linei ,int(vsini_star),
                                        int(vmac_star),R=int(R_star)) 
        
        ax=fig0.add_subplot(gs0[row0,int(start):int(start+step)])
        ax.axvline(wavei, color='red',ls='--', alpha=0.75)
        ax.plot(wave_st,flux_st,label=linei,c='navy')
        ax.plot(wavemodi,fluxmodi,alpha=0.75,c='green')
        
        ax.set_xlim(np.amin(wave_st),np.amax(wave_st))
        ax.set_ylim(fluxmin,fluxmax)
        ax.set_xlabel(r'$\lambda$ ($\AA$)')
        ax.legend(loc='upper right',fancybox=True, framealpha=0.7)
        
        if start<1:
            ax.set_ylabel('Normalized flux')
        else:
            ax.set_yticks([])

        start=start+step

    #FILA 2: HELIO
    print("HE LINES")

    row0=1
    start=0
    step=int(cols0/len(He_lineindex))
    fluxmin=1
    for c,ind in zip(range(len(He_lineindex)),He_lineindex):
        wavei=linewaves[ind]
        xrange=(wavei-10, wavei+10)
        flux_st=flux[(wave>xrange[0]) & (wave<xrange[1])]
        if np.nanmin(flux_st)<fluxmin:
            fluxmin=np.nanmin(flux_st)-0.05
            
    for c,ind in zip(range(len(He_lineindex)),He_lineindex):

        wavei=linewaves[ind]
        linei=linenames[ind]        
        xrange=(wavei-5, wavei+5)
        wave_st=wave[(wave>xrange[0]) & (wave<xrange[1])]
        flux_st=flux[(wave>xrange[0]) & (wave<xrange[1])]
        
        wavemodi,fluxmodi=broadened_FASTWIND_model(modeli,micro_final, linei ,int(vsini_star),
                                        int(vmac_star),R=int(R_star))
        
        ax=fig0.add_subplot(gs0[row0,int(start):int(start+step)])
        ax.axvline(wavei, color='red',ls='--', alpha=0.75)
        ax.plot(wave_st,flux_st,label=linei,c='navy')
        ax.plot(wavemodi,fluxmodi,alpha=0.75,c='green')
        
        ax.set_xlim(np.amin(wave_st),np.amax(wave_st))
        ax.set_ylim(fluxmin,1.3)
        ax.set_xlabel(r'$\lambda$ ($\AA$)')
        ax.legend(loc='upper right',fancybox=True, framealpha=0.7)
        
        if start<1:
            ax.set_ylabel('Normalized flux')
        else:
            ax.set_yticks([])            
        start=start+step

    #FILA 3: ELEMENTO

    print("%s LINES"%el_input)
        
    row0=2  
    start=0
    step=int(cols0/len(El_lineindex))
    fluxmin=1
    for c,ind in zip(range(len(El_lineindex)),El_lineindex):
        wavei=linewaves[ind]
        xrange=(wavei-5, wavei+5)
        flux_st=flux[(wave>xrange[0]) & (wave<xrange[1])]
        if np.nanmin(flux_st)<fluxmin:
            fluxmin=np.nanmin(flux_st)-0.05
            
    for c,ind in zip(range(len(El_lineindex)),El_lineindex):
        wavei=linewaves[ind]
        linei=linenames[ind]
        
        xrange=(wavei-5, wavei+5)
        wave_st=wave[(wave>xrange[0]) & (wave<xrange[1])]
        flux_st=flux[(wave>xrange[0]) & (wave<xrange[1])]
        
        wavemodi,fluxmodi=broadened_FASTWIND_model(modeli,micro_final, linei ,int(vsini_star),
                                        int(vmac_star),R=int(R_star))
        
        ax=fig0.add_subplot(gs0[row0,int(start):int(start+step)])
        ax.axvline(wavei, color='red',ls='--', alpha=0.75)
        ax.plot(wave_st,flux_st,label=linei,c='navy')
        ax.plot(wavemodi,fluxmodi,alpha=0.75,c='green')
        
        ax.set_xlim(np.amin(wave_st),np.amax(wave_st))
        ax.set_ylim(fluxmin,1.3)
        ax.set_xlabel(r'$\lambda$ ($\AA$)')
        ax.legend(loc='upper right',fancybox=True, framealpha=0.7)

        if start<1:
            ax.set_ylabel('Normalized flux')
        else:
            ax.set_yticks([])         

        #reconstructed gaussian fit

        RVi = star_measuretab["RV_%s"%linei][0]
        FWHMi = star_measuretab["FW_%s"%linei][0]
        DEPi = star_measuretab["dep_%s"%linei][0]

        wave_gaussian, flux_gaussian = gaussian_reconstruct(RVi,FWHMi,DEPi,wavei)

        ax.plot(wave_gaussian[(wave_gaussian>xrange[0]) & (wave_gaussian<xrange[1])],
                flux_gaussian[(wave_gaussian>xrange[0]) & (wave_gaussian<xrange[1])],
                "--",c="purple",alpha=0.9, linewidth=0.5)
        
        #!! igual se puede poner de colores diferntes si se mide "bien"  o no (como abel en sus plots de measure...)

        start=start+step    
        
    fig0.suptitle(r"%s (%s) log g=%.1f dex, Teff=%i K, log Q=-%.1f dex, $\beta$=%.1f, $\xi_t$=%d km/s, $\epsilon$(%s)=%.2f dex" 
                  %(star_name,spc_star,int(nearest_g)/100,int(nearest_T)*100, 
                  float(fixed_Q)/100,float(fixed_b)/100,int(micro_final),el_input,float(nearest_abun)),
                  fontsize=25)
    
    print("Saving the image...")
    fig0.savefig(savepath+ref_file.split(".")[0]+'_bestfit.png')    
    
    plt.show(block=False)

    return None


def find_FW_model(grid_path,el_input,T,g,He,Q,b,abun):
    '''
    Given a series of parameters, and the path of the original grid, this program 
    finds the FASTWIND OUT. file path with the information about that line.
    
    Parameters
    ----------
    grid_path : TYPE
        DESCRIPTION.
    el_input : TYPE
        DESCRIPTION.
    T : TYPE
        DESCRIPTION.
    g : TYPE
        DESCRIPTION.
    He : TYPE
        DESCRIPTION.
    Q : TYPE
        DESCRIPTION.
    b : TYPE
        DESCRIPTION.
    abun : TYPE
        DESCRIPTION.
        
    Returns
    -------
    model_path : TYPE
        DESCRIPTION.

    '''
    
    if el_input=='O':
        abun_centr='8.75'
    elif el_input=='N':
        abun_centr='7.80'
    elif el_input=='Mg':
        abun_centr='7.55'
    elif el_input=='H_He':
        abun_centr='7.80'
    elif el_input=='C':
        abun_centr='8.35'
    elif el_input=='Si':
        abun_centr='7.50'
    else:
        print('Please choose a valid element between O, N, Mg, H_He, C, Si')
        
    # abun=-3.4
    # abun_centr=-3.25
    
    abun=float(abun)
    abun_centr=float(abun_centr)
    listofabs=[-15,-10,-5,0,10,15]
    Aab=int(abun*100-abun_centr*100)
    Aab=listofabs[min(range(len(listofabs)),key=lambda i : abs(int(listofabs[i])-int(Aab)))]

    if Aab<0:
        Aab='m'+'0'*(3-len(str(abs(Aab))))+str(abs(Aab))
    elif Aab>0:
        Aab='p'+'0'*(3-len(str(abs(Aab))))+str(abs(Aab))
    else:
        Aab='p'+'0'*(3)

    he_i = "%d"%int(He)
    if len(he_i) < 2:
        he_i = "0%d"%int(He)

    path=fastwinddir+grid_path+'/'
    
    folder=glob.glob(path+'T'+str(int(T))+'g'+str(int(g))+'He'+he_i+'Q'+str(int(Q))+'b'+str(int(b))+'*'+Aab+'*')
    
    if len(folder)>1:
        print('ERROR: more than one model with these parameters')
        print(path+'T'+str(int(T))+'g'+str(int(g))+'He'+he_i+'Q'+str(int(Q))+'b'+str(int(b))+'*'+Aab+'*')
        return None
    elif len(folder)<1:
        print('ERROR: no model found with these parameters')
        print(path+'T'+str(int(T))+'g'+str(int(g))+'He'+he_i+'Q'+str(int(Q))+'b'+str(int(b))+'*'+Aab+'*')
        return None
    else:
        folder=folder[0] 

    return folder

class spec_FASTWIND_line():
    def __init__(self, path_model, micro, line):
        '''
        Parameters
        ----------
        path_model: str
            full path to the FASTWIND model in xdr format.
        micro: str
            microturbulence value as recorded in FASTWIND e.g. "001", "003", "025", etc.
            It needs to be a value that exists in the model in path_model. 
        line: str
            label of the line to read inside the FASTWIND model in xdr format. 
            It needs to be a label that exists in the model in path_model. 
        '''

        self.fullpath = path_model
        self.modelname = Path(self.fullpath).name
        self.micro = micro
        self.line = line

        #we could also add the rest of the parameters, but they are in the model name in self.path 

        Proof=scipy.io.readsav(path_model)
        Models=Proof['l']

        lines=[Models[jj][0].decode("utf-8") for jj in range(len(Models)) ]
        micros=[Models[jj][1].decode("utf-8") for jj in range(len(Models)) ] 

        mask=( (np.array(lines)==line) & (np.array(micros)==micro) ) 

        index=(np.argwhere(mask>0))[0][0]

        wave=Models[index][3]
        flux=Models[index][4]
 
        ind=int(len(wave) / 2)
        l0=wave[ind]       

        self.wave = wave
        self.flux = flux 
        self.l0 = l0
        self.dlam = np.mean(self.wave[1::] - self.wave[0:-1])
        self.resolution = None #should be >=0 or None. Try with np.nan

def broadened_FASTWIND_model(path_model,micro,line,vsini,vmac,R,beta=1.5):
        s = spec_FASTWIND_line(path_model, micro, line)

        stp = s.l0 / (R * 10)      

        f1 = interp1d(s.wave,s.flux, kind='linear')  
        wv = np.arange(start=s.wave[0],stop=s.wave[-2],step=stp)
        wv = wv[wv <= s.wave[-2]]
        flx = f1(wv) 

        s.wave = wv
        s.flux = flx
        s.dlam = np.mean(s.wave[1::] - s.wave[0:-1])

        spec.convolution(s, vsini = vsini, vmac = vmac, resol = R, quiet = "y")

        return s.wave, s.flux 

def read_OUT_lines(gridpath, folder, atom, micro):

    '''
    Reads the EW of a folder with FASTWIND models from the OUT.line files.
    Given the path where the files are, the atom of the lines you want, and
    the desired microturbulence, it gives a table. 
    '''

    folder = gridpath + folder

    lines = "OUT.%s*_VT00%d"%(atom,micro)

    folders = glob.glob(folder+lines)
    ew = np.zeros(len(folders), dtype=float)
    linelabel = np.empty(len(folders), dtype=object)
    ion =  np.zeros(len(folders), dtype=int)
    color = np.empty(len(folders), dtype=object)
    shape = np.empty(len(folders), dtype=object)

    for filei,i in zip(folders, range(len(folders))):
        linelabel[i] = ((filei.split("OUT.")[-1]).split("_VT")[0])
        if "IV" in linelabel[i]:
            ion[i] = 4
            color[i] = "purple"
            shape[i] = "s"
        elif "III" in linelabel[i]:
            ion[i] = 3
            color[i] = "green"
            shape[i] = "^"         
        elif "II" in linelabel[i]:
            ion[i] = 2
            color[i] = "brown"
            shape[i] = "d"

        with open(filei,'r') as file:
            file_list = file.readlines()
            try:
                ew[i] = float(file_list[-1])*-1000
            except:
                ew[i] = float(file_list[-2])*-1000

    t = Table([linelabel, ew, ion, color, shape],
              names= ["linelabel","ew", "ion", "color", "shape"])

    return t   


def calc_Tg(ref_file,spc_star,orig,el_input,star_measure,syn_measure,lines,
            vtable,savepath,fixed_Q='1400',fixed_He='10',fixed_b='10',
            fixed_micro='007',fixed_abun='centr', tol=50,save='y', show="y"):
    
    '''
    Programa que calcula las mejores parejas de T, g para una estrella según 
    una red de T, g, con el método de anchuras equivalentes.
    
    Parameters
    ----------
    star_name : str
        El nombre de la estrella a analizar. 
        
    spc_star: str

    orig: str
        Supports 'IACOB' and 'txt'        
        
    star_measure: str
         Nombre de la tabla en la que se encuentran los ficheros con las medidas de EW de la estrella.
    
    lines: str
        Enter the wavelenght(s) of the line(s) to show, either in a coma-separated
        string, or in a .txt/.lst file containing the lines.
        The default is lines_T_g.lst
 
    star_grid: str
        Grid de espectros sintéticos con la que se trabajará para hacer el análisis.
        Actualmente las opciones que existen son:
            'BDw': grid 1, óptima para B Dwarfs.
            'BSg': grid 2, óptima para B Super-Giants.
        Si el input no es ninguna de estas dos opciones, dará a elegir manualmente entre ellas.
        
    spc_star: str
        Tipo espectral de la estrella.
    
    vtable: str
        Nombre de la tabla en la que están medidas las velocidades de rotación y de
        macroturbulencia de la estrella.
        
    savepath: str
        Path donde se guardan las figuras y tablas generadas. 
        
    orig: str
        'IACOB' para tratar espectros de la IACOB
        'synthetic' para tratar espectros sintéticos.
        
    tol: int
        Sets the tolerance [km/s] to shifting the spectrum in order to fit the lines 
        after being corrected from radial velocity.
        Default is 50. 
        
    logQ_ref: float
        Valor del parámetro del viento logQ medido con MAUI.
        Default is -13.5
        
    save: str
        If 'y', saves all the plots and tables. 
        If 'n', saves nothing, just gives the parameters as output

    Returns
    -------    
    output: Table 
        Tabla con las parejas de T, g que cumplen el equilibrio de ionización para la estrella dada.
        Esta tabla también se guarda en star_name+'_Tg_EWSi.ecsv'. 
    
    La función también guarda un archivo zcrea_Tg.pro para ejecutar en IDL y descargar con MAUI
    la red de espectros sintéticos necesaria para el siguiente paso. 
    '''

    #Primero visualizamos tres líneas de Silicio  de cada tipo.
    flag = "NO_ANALYSIS"

    syn_measure=findtable(syn_measure,path=maindir+'tables/')

    syn_measure = syn_measure[(syn_measure["Grid_name"] != "symetric-interp") &
                                        (syn_measure["Grid_name"] != "asymetric-interp")]  
    # syn_measure=syn_measure[syn_measure['Grid_name']==syn_measure['Grid_name'][0]]

    star_measure=findtable(star_measure,path=maindir+'tables/') 
    star_measure=star_measure[star_measure['Ref_file']==ref_file][0]

    star_name = star_measure["ID"]

    vtable=findtable(vtable,path=maindir+'tables/')
        
    if len(vtable[vtable['Ref_file']==ref_file])<1:
        print('\nNo vsini or vmac info known about the star.')
        vsini_star=50
        vmac_star=0
    else:
        if 'vsini_GF' in vtable.colnames and 'vmac_GF' in vtable.colnames:    
            vsini_star=vtable['vsini_GF'][vtable['Ref_file']==ref_file][0]
            vmac_star=vtable['vmac_GF'][vtable['Ref_file']==ref_file][0]
        elif 'vsini' in vtable.colnames and 'rtmac' in vtable.colnames:
            vsini_star=vtable['vsini'][vtable['Ref_file']==ref_file][0]
            vmac_star=vtable['rtmac'][vtable['Ref_file']==ref_file][0]
        elif 'vsini' in vtable.colnames and 'vmac' in vtable.colnames:
            vsini_star=vtable['vsini'][vtable['Ref_file']==ref_file][0]
            vmac_star=vtable['vmac'][vtable['Ref_file']==ref_file][0]            
        else:
            print('Please check the column names in vtable for vsini and rtmac')
            print('\nNo vsini or vmac info known about the star.')
            vsini_star=50
            vmac_star=0  
            
        print('\nvsini=%.2f, vmac=%.2f'%(vsini_star,vmac_star))

    star=spec(ref_file, snr='bestHF', orig=orig,rv0=star_measure['RV0']) 
    
    wave=star.wave
    flux=star.flux
    
    lines=findlines(lines)
    linewaves,linenames= lines[0], lines[1]
    
    Si_lineindex=[index for index, i in enumerate(linenames) if i[0] == 'S']

    print("Using the following lines for the ionization equilibrium:")
    for ind in Si_lineindex:
        print(linenames[ind])

    rows=1
    cols=len(Si_lineindex)
    
    fig,axs=plt.subplots(rows,cols,figsize=(20,6))
    EWline_array=np.empty(len(Si_lineindex),dtype=float)
        
    for i,ind in zip(range(len(Si_lineindex)),Si_lineindex):
        xrange=(linewaves[ind]-30, linewaves[ind]+30)
        ylow=1-5*np.std(flux[(wave>xrange[0]) & (wave<xrange[1])]) 
        yhigh=1+2*np.std(flux[(wave>xrange[0]) & (wave<xrange[1])]) 
        yrange=(ylow,yhigh)
        
        axs[i].plot(wave,flux,label=linenames[ind])
        axs[i].axvline(linewaves[ind], color='red',ls='--', alpha=0.75)
        axs[i].set_xlim(xrange[0],xrange[1])
        axs[i].set_ylim(0.5,1.2) 
        axs[i].legend(loc='upper right')
        axs[i].set_xlabel(r'$\lambda$ ($\AA$)')
        
        if (~np.ma.getmask(star_measure['EW_'+str(linenames[ind])])):
            EWline=star_measure['EW_'+str(linenames[ind])][(star_measure['EW_'+str(linenames[ind])]>5) & 
                                                      (star_measure['dep_'+str(linenames[ind])]>0.005) &
                                                     (np.abs(star_measure['RV_'+str(linenames[ind])]-star_measure['RV0'])<0.5*tol)]
        
            if len(EWline)<1:
                EWline=0
            else: 
                EWline = EWline[0]

        else:
            EWline=0

        EWline_array[i]=EWline
        
        axs[i].text(linewaves[ind],0.6,'EW=%.1f' %EWline)
    
    axs[0].set_ylabel('Norm.  flux')      
    fig.suptitle('%s (%s) \n %s'%(star_name,spc_star, ref_file),fontsize='30')        
    plt.tight_layout()
    plt.show(block=False)
    fig.savefig(savepath+star_name+'_Silines.png')

    stop = input("Press return to continue with the analysis.")

    si_selection=''
        
    # tol_FW=1.5
    while si_selection not in['1','2','3']:
        if ( (EWline_array[0]>5) and (EWline_array[1]>5) and (EWline_array[2]<5) ):
            si_selection='2'
        elif ( (EWline_array[0]>5) and (EWline_array[1]>5) and (EWline_array[2]>5) and (EWline_array[0]>EWline_array[2])): 
            si_selection='2'

        elif ( (EWline_array[0]<5) and (EWline_array[1]>5) and (EWline_array[2]>5) ):
            si_selection='1'
        elif ( (EWline_array[0]>5) and (EWline_array[1]>5) and (EWline_array[2]>5) and (EWline_array[2]>=EWline_array[0])): 
            si_selection='1'
        
        else:
            si_selection = '0'
            # si_selection=input('Choose between SiIV4116/SiIII4552 (1) or SiIII4552/SiII6347 (2): ')
            print("Chosen lines are not strong enough to calculate the Teff, logg at the ionization equilibrium for this star. ")
            raise Exception("Chosen lines are not strong enough to calculate the ionization equilibrium for this star.") 

        if si_selection=='1':
            if show=="y": print('Analyzing SiIV/SiIII')
            line2='EW_SiIV4116'
            # line2 = "EW_SiIV4631"
            line1 = 'EW_SiIII4552'

        elif si_selection=='2':
            if show=="y": print('Analyzing SiIII/SiII')
            line2='EW_SiIII4552'
            line1='EW_SiII6371'

        else:
            print('Please press 1 or 2')

    line3= 'EW_SiIII4567'
    line4= 'EW_SiIII4574'

    T_lines = "%s/%s"%(line2.split("_")[-1],line1.split("_")[-1])

    EW2_star=float(star_measure[line2])
    EW1_star=float(star_measure[line1])
    EW3_star=float(star_measure[line3])
    EW4_star=float(star_measure[line4])
    
    #cambiar errores a 5 -> CONSISTENCIA. 
    EW2_star_err=5
    EW1_star_err=5

    EWstar=EW2_star/EW1_star
    EWstar_err=np.sqrt((EW2_star_err/EW1_star)**2+(EW2_star*EW1_star_err/(EW1_star)**2)**2)    

    ##### PRUEBA 2025-08-29 ###########
    # log = False
    # if log == True: 
    #     EWstar = np.log(EWstar)
    #     EWstar_err = (EWstar_err)/EWstar  #! COMPROBAR
        
    abun=[x for x in syn_measure.colnames if 'abun' in x][0]
    
    if fixed_abun=='centr':
        if el_input=='O':
            fixed_abun='8.75'
        elif el_input=='N':
            fixed_abun='7.80'
        elif el_input=='Mg':
            fixed_abun='7.55'
        elif el_input=='H_He':
            fixed_abun='7.80'
        elif el_input=='C':
            fixed_abun='8.35'
        elif el_input=='Si':
            fixed_abun='7.50'
            
    else:
        
        fixed_abun=str(fixed_abun)
   
    mask=((syn_measure['He']==str(fixed_He)) & (syn_measure['Q']==str(fixed_Q)) & 
      (syn_measure[abun]==fixed_abun) & (syn_measure['micro']==str(fixed_micro))
      & (syn_measure['b']==str(fixed_b) ) )     

    syn_measure = syn_measure[mask]

    T=np.array(sorted(set(syn_measure['T'])),dtype=int)*100
    g=sorted(set(syn_measure['g']))
    Ag=10 #error en g

    fig = plt.figure (figsize = (10,14),constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[3,2,2])
    ax0=fig.add_subplot(gs[0,:])
    ax1 = fig.add_subplot(gs[1,0])
    ax2 = fig.add_subplot(gs[2,0])
    ax3 = fig.add_subplot(gs[1,1])
    ax4 = fig.add_subplot(gs[2,1])

    colors = {"black": "#000000", "orange" : "#E69F00", "sky blue": "#56B4E9", "bluish green": "#009E73", 
            "yellow": "#F0E442", "blue": "#0072B2", "vermillion": "#D55E00", "reddish purple": "#CC79A7", 
            "light orange": "#ffd881", "light green" : "#00f6b3"}    

    color_list = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00',
                  '#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']
                 
    EWline_array=np.empty(len(Si_lineindex),dtype=float)

    ylow_all = 1
    yhigh_all = 1
        
    for i,ind in zip(range(len(Si_lineindex)),Si_lineindex):
        
        if (~np.ma.getmask(star_measure['EW_'+str(linenames[ind])])):
            EWline=star_measure['EW_'+str(linenames[ind])][(star_measure['EW_'+str(linenames[ind])]>5) & 
                                                      (star_measure['dep_'+str(linenames[ind])]>0.005) &
                                                     (np.abs(star_measure['RV_'+str(linenames[ind])]-star_measure['RV0'])<0.5*tol)]
            if len(EWline)<1:
                EWline=0
            else: 
                EWline = EWline[0]
        else:
            EWline=0

        EWline_array[i]=EWline    
        
        if EWline > 0 : 
            RV_line = star_measure['RV_'+str(linenames[ind])]
            RV_star= star_measure['RV0']
            v = RV_line - RV_star


            c = 299792.458  # km/s

            lambda0 = linewaves[ind] * (1 + v / c)

            sigma = star_measure['FW_'+str(linenames[ind])] / (2.35482)
            depth = (star_measure['EW_'+str(linenames[ind])]/1000) / (sigma * np.sqrt(2*np.pi))

            xrange=(linewaves[ind]-5, linewaves[ind]+5)
            wave_gaussian = np.linspace(lambda0-5, lambda0+5, 500)
            flux_gaussian = 1 - depth * np.exp(-(wave_gaussian - lambda0)**2 / (2*sigma**2))

        else: 
            xrange=(linewaves[ind]-5, linewaves[ind]+5)   
        
        ylow=1-5*np.std(flux[(wave>xrange[0]) & (wave<xrange[1])]) 
        yhigh=1+2*np.std(flux[(wave>xrange[0]) & (wave<xrange[1])]) 

        if ylow < ylow_all: 
            ylow_all = ylow
        if yhigh > yhigh_all:
            yhigh_all = yhigh

    EWmin=1000
    EWmax=15

    log = False
    # log = True
    if log == True:
        EWmax = np.log(EWmax)


    Tmin=np.amax(T)
    Tmax=np.amin(T)

    flag1 = np.array(["OK"]*len(g))
    flag2 = np.array(["OK"]*len(g))

    flag = "OK"

    mask_ind=np.round(np.linspace(0,len(g)-1,10)).astype(int)

    count=1

    for i,gi in zip(range(len(g)),g):

        EW_2=syn_measure[line2][(syn_measure['g']==gi)]
        EW_1=syn_measure[line1][(syn_measure['g']==gi)]
        EW_3 = syn_measure[line3][(syn_measure['g']==gi)]
        EW_4 = syn_measure[line4][(syn_measure['g']==gi)]
        Ti=np.array((syn_measure['T'][(syn_measure['g']==gi)]),dtype=int)*100
               
        maski= ((EW_2>5) & (EW_1>5))
        EW_2=EW_2[maski]
        EW_1=EW_1[maski]
        EW_3 = EW_3[maski]
        EW_4 = EW_4[maski]
        Ti=Ti[maski]

        if len(EW_2/EW_1)<2:
            flag1[i] = "X"
            flag2[i] = "X"
            continue
        elif min(EW_2/EW_1)>EWstar+EWstar_err:
            flag1[i] = "X"
            flag2[i] = "X"
            continue
        elif max(EW_2/EW_1)<EWstar+EWstar_err:
            flag1[i] = "X"
            flag2[i] = "X"            
            continue    
        else:

            if int(gi) in np.array(g, dtype=int)[mask_ind]: 

                if log == True:
                    ax0.plot(Ti,np.log(EW_2/EW_1),'.-',label='log g=%.1f' %(float(gi)/100), c = color_list[count])

                    count = count+1
                    if min(np.log(EW_2/EW_1))<EWmin:
                        EWmin=min(np.log(EW_2/EW_1))
                        EWmax=EWstar+np.abs(EWstar-EWmin)              
                else:
                    ax0.plot(Ti,EW_2/EW_1,'.-',label='log g=%.1f' %(float(gi)/100))
             
                                                      
                    count = count+1
                    if min(EW_2/EW_1)<EWmin:
                        EWmin=min(EW_2/EW_1)
                        EWmax=EWstar+np.abs(EWstar-EWmin)                  
                            
                ax1.plot(Ti,EW_2,'.-',label='log g=%.1f' %(float(gi)/100))
                ax1.plot(np.linspace(Tmin-1000,Tmax+1000,100),EW2_star*np.ones(100),'r--')
                ax2.plot(Ti, EW_1,'.-',label='log g=%.1f' %(float(gi)/100))  
                ax2.plot(np.linspace(Tmin-1000,Tmax+1000,100),EW1_star*np.ones(100),'r--')
                ax3.plot(Ti,EW_3,'.-',label='log g=%.1f' %(float(gi)/100))
                ax3.plot(np.linspace(Tmin-1000,Tmax+1000,100),EW3_star*np.ones(100),'r--')
                ax4.plot(Ti, EW_4,'.-',label='log g=%.1f' %(float(gi)/100))  
                ax4.plot(np.linspace(Tmin-1000,Tmax+1000,100),EW4_star*np.ones(100),'r--')                                    
                    
                Timin=Ti[np.nanargmin(np.abs(np.array(EW_2/EW_1)-EWmax))]
                Timin=np.min(Ti)
                Timax=np.max(Ti)
                if Timin<Tmin:
                    Tmin=Timin
                if Timax>Tmax:
                    Tmax=Timax

            if np.sum(np.array(EW_1, dtype=float) > float(EW1_star)) < 1: 
                flag1[i] = "EXTRAPOL"
            if np.sum(np.array(EW_2, dtype=float) > float(EW2_star)) < 1:
                flag2[i] = "EXTRAPOL"   

    Delta_Ti = np.abs(int(sorted(set(syn_measure["T"]))[1])-int(sorted(set(syn_measure["T"]))[0]))

    if np.sum(np.array(flag1) =="OK") <1:
        flag = "EXTRAPOL"
    if np.sum(np.array(flag2)=="OK") <1:
        flag = "EXTRAPOL"           
    
    ax0.plot(np.linspace(Tmin-1000,Tmax+1000,100),EWstar*np.ones(100),'r--')
    ax0.fill_between(np.linspace(Tmin-1000,Tmax+1000,100), (EWstar+EWstar_err)*np.ones(100), 
                    (EWstar-EWstar_err)*np.ones(100), alpha=.25,color='yellow')  
    
    ax0.set_xlim(Tmin-100,Tmax+100)
    ax0.set_ylim(EWmin,EWmax)  #PLOT SIMÉTRICO POR ARRIBA Y POR ABAJO:
    ax0.set_xlabel('Teff (K)')
    if log == True:
        ax0.set_ylabel(" log(%s / %s)" %(line2,line1))
    else:
        ax0.set_ylabel("%s / %s" %(line2,line1))
    ax1.set_ylabel(line2)
    ax2.set_ylabel(line1)
    ax3.set_ylabel(line3)
    ax4.set_ylabel(line4)    
    ax2.set_xlabel('Teff (K)')
    ax4.set_xlabel('Teff (K)')

    ax0.legend(fancybox=True, framealpha=0.7,fontsize="15") #,loc='upper left'
    ax0.set_title("He=%d, Q=%d, b=%d, abSi=%.2f, micro=%d" %(int(fixed_He),int(fixed_Q),int(fixed_b),float(fixed_abun),int(fixed_micro)))

    fig.suptitle('%s (%s) \n %s'%(star_name,spc_star, ref_file),fontsize=26)
        
    output = Table(
            names = (['g','Ag','Teff','ATeff+','ATeff-']),
            dtype = (['int64','int64','int64','int64','int64'] ) ) 
    
    for i,gi in zip(range(len(g)),g):

        EW_2=syn_measure[line2][(syn_measure['g']==gi) ]
        EW_1=syn_measure[line1][(syn_measure['g']==gi) ]
        Ti=np.array((syn_measure['T'][(syn_measure['g']==gi) ]),dtype=int)*100
            
        maski= ((EW_2>5) & (EW_1>5))
        EW_2=EW_2[maski]
        EW_1=EW_1[maski]
        Ti=Ti[maski]  
        
        if len(EW_2/EW_1)<4:
            continue
        elif min(EW_2/EW_1)>EWstar: 
            continue
        elif max(EW_2/EW_1)<EWstar: 
            continue        

        else:
            # normi=np.abs((EW_2/EW_1)**2-(EWstar)**2)

            if log == True: 
                normi=np.abs(np.log(EW_2/EW_1)-(EWstar))
                T0=Ti[np.nanargmin(normi)]
                EW0=np.log(EW_2[np.nanargmin(normi)]/EW_1[np.nanargmin(normi)])

                EW_2 = np.delete(EW_2,np.nanargmin(normi))
                EW_1 = np.delete(EW_1,np.nanargmin(normi))
                Ti = np.delete(Ti,np.nanargmin(normi))

                T1=Ti[min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0)))]
                EW1 = np.log(EW_2[min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0)))]/EW_1[min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0)))])
                EW_2 = np.delete(EW_2,min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0))))
                EW_1 = np.delete(EW_1,min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0))))
                Ti = np.delete(Ti,min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0))))

                T2=Ti[min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0)))]
                EW2 = np.log(EW_2[min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0)))]/EW_1[min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0)))])
                EW_2 = np.delete(EW_2,min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0))))
                EW_1 = np.delete(EW_1,min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0))))
                Ti = np.delete(Ti,min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0))))

                T3=Ti[min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0)))]
                EW3 = np.log(EW_2[min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0)))]/EW_1[min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0)))])
                EW_2 = np.delete(EW_2,min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0))))
                EW_1 = np.delete(EW_1,min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0))))
                Ti = np.delete(Ti,min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0))))

            else: 
                normi=np.abs((EW_2/EW_1)-(EWstar))
                T0=Ti[np.nanargmin(normi)]
                EW0=(EW_2[np.nanargmin(normi)]/EW_1[np.nanargmin(normi)])

                EW_2 = np.delete(EW_2,np.nanargmin(normi))
                EW_1 = np.delete(EW_1,np.nanargmin(normi))
                Ti = np.delete(Ti,np.nanargmin(normi))

                T1=Ti[min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0)))]
                EW1 = (EW_2[min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0)))]/EW_1[min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0)))])
                EW_2 = np.delete(EW_2,min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0))))
                EW_1 = np.delete(EW_1,min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0))))
                Ti = np.delete(Ti,min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0))))

                T2=Ti[min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0)))]
                EW2 = (EW_2[min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0)))]/EW_1[min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0)))])
                EW_2 = np.delete(EW_2,min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0))))
                EW_1 = np.delete(EW_1,min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0))))
                Ti = np.delete(Ti,min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0))))

                T3=Ti[min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0)))]
                EW3 = (EW_2[min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0)))]/EW_1[min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0)))])
                EW_2 = np.delete(EW_2,min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0))))
                EW_1 = np.delete(EW_1,min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0))))
                Ti = np.delete(Ti,min(range(len(Ti)),key=lambda i : abs(int(Ti[i])-int(T0))))                

            list1=[T0,T1,T2,T3]
            list2=[EW0,EW1,EW2,EW3]
            list2, list1 = (list(t) for t in zip(*sorted(zip(list2, list1))))
       
            T_star=np.interp(EWstar, list2, list1)
            
            T_up=np.interp(EWstar+EWstar_err, list2, list1)
            T_down=np.interp(EWstar-EWstar_err, list2, list1)
            if int(gi) in np.array(g, dtype=int)[mask_ind]:
                ax0.plot(list1, list2, 'k.-', alpha=0.5, zorder=400)
                ax0.errorbar(T_star, EWstar, marker ="*", zorder=500, alpha=0.7, mfc="black", mec="black", ms = 15)

            ###################################################################

            ATeff_up = np.abs(T_up-T_star) if np.abs(T_up-T_star)<Delta_Ti else Delta_Ti
            ATeff_down = np.abs(T_down-T_star) if np.abs(T_down-T_star)<Delta_Ti else Delta_Ti
            
            T_source = Table([[gi],[Ag],[T_star],[np.abs(T_up-T_star)],[np.abs(T_down-T_star)]],
                             names=(['g','Ag','Teff','ATeff+','ATeff-']),
                             dtype = (['int64','int64','int64','int64','int64']))
            
            output = vstack([output,T_source])

    if save=='y':
        fig.savefig(savepath+star_name+'EW2EW1_vs_Teff.png')

    plt.show(block=False)       
                
    if show=="y":
        print('Possible g, T, combinations: \n', output)

    if len(output)<1:
        flag = "GRID"

    return output, flag, T_lines

#%%

def best_Tg(ref_file,spc_star,orig,table_calcTg,
        star_measure_file,FW_measure,vtable_file,parameterstable,savepath,
        fixed_Q='1400',fixed_He='10',fixed_b='10',
        fixed_micro='007',fixed_abun='centr',beta_fixed = False, save='y', show="y"):
    
    el_input='Si'
    
    if isinstance(table_calcTg, str):
        table_calcTg=findtable(table_calcTg,path=savepath,format='ecsv',delimiter='\s')
    
    FW_measure=findtable(FW_measure)
    
    # FW_measure=FW_measure[FW_measure['Grid_name']==FW_measure['Grid_name'][0]]
    FW_measure = FW_measure[(FW_measure["Grid_name"] != "symetric-interp") &
                                        (FW_measure["Grid_name"] != "asymetric-interp")]
    
    if len(sorted(set(FW_measure['b'])))<2:
        beta_fixed = True

    if beta_fixed: 
        FW_measure = FW_measure[(FW_measure['b']==fixed_b)]
    else:
        FW_measure = FW_measure[(FW_measure['b']=="10") | (FW_measure['b']=="20")]

    star_measure=findtable(star_measure_file,path=maindir+'tables/')
    star_measure=star_measure[star_measure['Ref_file']==ref_file][0]
    star_name = star_measure["ID"]

    vtable=findtable(vtable_file,path=maindir+'tables/')
    
    if len(vtable[vtable['Ref_file']==ref_file])<1:
        print('\nNo vsini or vmac info known about the star.')
        print('Assuming vsini=50, vmac=0')
        vsini_star=50 
        vmac_star=0 
    else:
        if 'vsini_GF' in vtable.colnames and 'vmac_GF' in vtable.colnames:    
            vsini_star=vtable['vsini_GF'][vtable['Ref_file']==ref_file][0]
            vmac_star=vtable['vmac_GF'][vtable['Ref_file']==ref_file][0]
        elif 'vsini' in vtable.colnames and 'rtmac' in vtable.colnames:
            vsini_star=vtable['vsini'][vtable['Ref_file']==ref_file][0]
            vmac_star=vtable['rtmac'][vtable['Ref_file']==ref_file][0]
        elif 'vsini' in vtable.colnames and 'vmac' in vtable.colnames:
            vsini_star=vtable['vsini'][vtable['Ref_file']==ref_file][0]
            vmac_star=vtable['vmac'][vtable['Ref_file']==ref_file][0]                 
        else:
            print('Please check the column names in vtable for vsini and rtmac')
            print('\nNo vsini or vmac info known about the star.')
            print('Assuming vsini=50, vmac=0')
            vsini_star=50
            vmac_star=0  
            
        if show == "y": print('vsini=%.2f, vmac=%.2f'%(vsini_star,vmac_star))
    
    rv0=star_measure['RV0']  

    lines=([6562.8,4861.325,4340.463,4101.735, 3970.072, #! PROBANDO HEPS (2025-10-27)
            4387.929,4471.474,4713.139,4199.832,4541.591,5411.52],
           ['HALPHA','HBETA','HGAMMA','HDELTA', 'HEPS',
            'HEI4387','HEI4471','HEI4713','HEII4200','HEII4541','HEII5411'])

    linewaves,linenames= lines[0], lines[1]

    # #INCLUDING HEPS:
    Balmer_lineindex=[index for index, i in enumerate(linenames) if (i[0]=='H') and (i[1:3]!='EI')]
    He_lineindex=[index for index, i in enumerate(linenames) if (i[0]=='H') and (i[1:3]=='EI')]
    
    window_Balmer = {'HALPHA':[6515.00 ,6605.25],
                     'HBETA':[4845.40, 4875.30],
                     'HGAMMA':[4321.27, 4360.02],
                     'HDELTA':[4082.02, 4113.80],
                     'HEPS':[3950.50,3990.50]}

    #WE FIX THE PARAMETERS  
    garr=np.array(table_calcTg['g'],dtype=int)
    Tarr=np.zeros(len(table_calcTg['Teff']),dtype=int)
    for j in range(len(table_calcTg['Teff'])):
        T_new=table_calcTg['Teff'][j] 
        T_list=np.array(sorted(set(FW_measure['T'][FW_measure['g']==str(garr[j])])),dtype=int)*100
        nearest=T_list[min(range(len(T_list)),key=lambda i : abs(int(T_list[i])-int(T_new)))]
        Tarr[j]=str(int(nearest/100))
    
    if fixed_abun=='centr':
        if el_input=='O':
            fixed_abun='8.75'
        elif el_input=='N':
            fixed_abun='7.80'
        elif el_input=='Mg':
            fixed_abun='7.55'
        elif el_input=='H_He':
            fixed_abun='7.80'
        elif el_input=='C':
            fixed_abun='8.35'
        elif el_input=='Si':
            fixed_abun='7.50'
    else:
        fixed_abun=str(fixed_abun)

    #!! WE DO SIGMA CLIPPING TO TAKE AWAY THE PEAKS THAT DO NOT BELONG TO THE BALMER LINES. 

    format_i = "csv" if ".csv" in parameterstable else "ecsv" if ".ecsv" in parameterstable else "ascii"
    delimiter_i = "," if ".csv" in parameterstable else "\s" if ".ecsv" in parameterstable else "\s"    
    final_parameters=findtable(parameterstable,path=savepath,format=format_i,delimiter=delimiter_i)

    if final_parameters is not None:
        if ref_file in final_parameters["Ref_file"]:
            final_parameters = final_parameters[final_parameters["Ref_file"] == ref_file]
            T_temp=final_parameters['T'][-1]
            g_temp=final_parameters['g'][-1]
        else:
            T_temp = table_calcTg["Teff"][int((len(table_calcTg) - 1)/2)]
            g_temp = table_calcTg["g"][int((len(table_calcTg) - 1)/2)]
    else:    
        T_temp = table_calcTg["Teff"][int((len(table_calcTg) - 1)/2)]
        g_temp = table_calcTg["g"][int((len(table_calcTg) - 1)/2)]     

    T_list=np.array(sorted(set(FW_measure['T'][FW_measure['g']==str(g_temp)])),dtype=int)*100
    nearest=T_list[min(range(len(T_list)),key=lambda i : abs(int(T_list[i])-int(T_temp)))]
    T_temp = str(int(nearest/100))
    g_temp = str(int(g_temp))        

    # READING THE SPECTRA OF THE STAR: 
    star=spec(ref_file, snr='bestHF', orig=orig,rv0=rv0)
    wave=star.wave
    flux=star.flux

    ###################################################################################
    #### FUNCTIONS FOR MASKING THE SPECTRUM TO ONLY KEEP THE BALMER LINES #############
    ###################################################################################
    
    def process_balmer_peaks(ref_file, orig, rv0, snr_key, Balmer_lineindex,
                            linewaves, linenames, FW_measure, el_input,
                            T_temp, g_temp, fixed_He, fixed_Q, fixed_b,
                            fixed_abun, fixed_micro, vsini_star, vmac_star,
                            savepath, star_name, spec, find_FW_model,
                            broadened_FASTWIND_model,
                            nsigma=3.0, max_iter=8, R=85000, beta=1.5,
                            kernel_mult=4.0, min_width_frac=0.6, pad_frac=0.5,
                            min_median_kernel_aa=1.0, min_min_group_aa=0.05,
                            min_pad_aa=0.1):
        
        print("Working on removing the contaminating lines in the Balmer line profiles...")

        star = spec(ref_file, snr=snr_key, orig=orig, rv0=rv0)
        wave = star.wave
        flux = star.flux

        mask_peakproc = np.ones(len(wave), dtype=bool)
    
        show_plot = "n" #change to "y" whenever relevant. 

        rows = 2
        cols = len(Balmer_lineindex)

        if show_plot == "n": plt.ioff()  # turn off interactive mode

        fig_sc, axs_sc = plt.subplots(rows, cols, figsize=(9 * cols, 6 * rows))
        if cols == 1:
            axs_sc = axs_sc.reshape(rows, 1)
    
        for i, ind in zip(range(len(Balmer_lineindex)), Balmer_lineindex):
    
            # exclude the line core itself from the "other lines" search
            mask_peakproc[(wave > linewaves[ind] - 1.5) &
                        (wave < linewaves[ind] + 1.5)] = False
            
            
    
            mask_st = (wave > linewaves[ind] - 30) & (wave < linewaves[ind] + 30)
            wave_st = wave[mask_st]
            flux_st = flux[mask_st]
            pixel_scale = np.median(np.diff(wave_st))


            modeli = find_FW_model(FW_measure['Grid_name'][FW_measure["He"] == str(fixed_He)][0], 
                                   el_input, T_temp,
                                    g_temp, fixed_He, fixed_Q, fixed_b, fixed_abun)

            wavemodi, fluxmodi = broadened_FASTWIND_model(
                modeli, fixed_micro, linenames[ind], int(vsini_star),
                int(vmac_star), R=R)
            fit_syn = np.interp(wave_st, wavemodi, fluxmodi)

            local_mask = mask_peakproc[mask_st]
            local_mask[((fit_syn > 0.99) & (fit_syn < 1.01))] = False
    
            FW_cols = [coli for coli in star_measure.colnames if "FW" in coli] 

            FW_vals = np.ma.masked_array(list(star_measure[FW_cols]), dtype=float)
            W = np.nanmin(FW_vals.filled(np.nan))

            W =np.nanmin( np.ma.masked_array(list(star_measure[FW_cols])).filled(np.nan) )
            median_kernel_aa_line = max(min_median_kernel_aa, kernel_mult * W)
            min_group_aa_line = max(min_min_group_aa, min_width_frac * W)
            pad_aa_line = max(min_pad_aa, pad_frac * W)

            median_kernel_px = max(3, int(round(median_kernel_aa_line / pixel_scale)))
            if median_kernel_px % 2 == 0:
                median_kernel_px += 1  # odd kernel size

            smoothed = gaussian_filter1d(flux_st, sigma=10)  # light denoise only
            resid_smooth = fit_syn - flux_st
    
            axs_sc[0, i].plot(wave_st, resid_smooth, c='blue', lw=0.7)

            med0 = np.median(resid_smooth)
            sig0 = 1.4826 * np.median(np.abs(resid_smooth - med0))

            for k in (1, 2, 3):
                axs_sc[0, i].axhline(med0 + k * sig0, color='k', ls='--', alpha=0.2 / k * 2)
                axs_sc[0, i].axhline(med0 - k * sig0, color='k', ls='--', alpha=0.2 / k * 2)

            axs_sc[0, i].axhline(med0, color='k')
            axs_sc[0, i].set_title("%s  (W=%.3f AA)" % (linenames[ind], W))
            
            if pixel_scale and pixel_scale > 0:
                pad_px = max(1, int(round(pad_aa_line / pixel_scale)))
                if min_group_aa_line is not None:
                    min_group = max(1, int(round(min_group_aa_line / pixel_scale)))
            else:
                pad_px = max(1, int(round(pad_aa_line)))
        
            for n_1 in range(1, max_iter + 1):
                good = resid_smooth[local_mask]
                if good.size < 5:
                    break
        
                med = np.median(good)
                sigma = 1.4826 * np.median(np.abs(good - med))
                # sigma = np.std(good) 
                if sigma == 0:
                    break
        
                outlier = ((resid_smooth > med + nsigma * sigma) |
                        (resid_smooth < med - nsigma * sigma))
                
                candidate_bad = outlier & local_mask
                if not candidate_bad.any():
                    break

                #keep only contiguous runs long enough to be a real feature: 

                labels, n_groups = label(candidate_bad)
                confirmed_bad = np.zeros_like(candidate_bad)

                for k in range(1, n_groups + 1):
                    run = (labels == k)
                    if run.sum() >= min_group:
                        confirmed_bad |= run

                if not confirmed_bad.any():
                    break

                if pad_px > 0:
                    confirmed_bad = binary_dilation(confirmed_bad, iterations=pad_px)
                
                new_keep = local_mask & ~confirmed_bad
                

                if np.array_equal(new_keep, local_mask):
                    break
                
                local_mask = new_keep
                
    
            axs_sc[1, i].plot(wave_st, flux_st, c='b', lw=0.8)
            axs_sc[1, i].plot(wave_st, smoothed, c='orange', alpha=0.75, lw=0.8)
            axs_sc[1, i].plot(wave_st, fit_syn, c='r', alpha=0.4, lw=0.8)
            axs_sc[1, i].fill_between(wave_st, 0.4, 1.4, where=local_mask,
                                    alpha=.25, color='cyan', edgecolor="none")
            axs_sc[1, i].set_ylim(0.3, 1.4)
            axs_sc[1, i].set_title("%d passes"%(n_1))

            mask_peakproc[mask_st] = local_mask
    
        fig_sc.tight_layout()

        if show_plot == "y":
            fig_sc.savefig(savepath + star_name + 'balmer_peakprocess.png')
        else: 
            plt.ion()
        plt.close(fig_sc)
        return mask_peakproc, wave, flux, star

    ###################################################################################

    mask_peakproc, wave, flux, star = process_balmer_peaks(
        ref_file, orig, rv0, 'bestHF', Balmer_lineindex,
        linewaves, linenames, FW_measure, el_input,
        T_temp, g_temp, fixed_He, fixed_Q, fixed_b,
        fixed_abun, fixed_micro, vsini_star, vmac_star,
        savepath, star_name, spec, find_FW_model,
        broadened_FASTWIND_model,
        nsigma=3.0, R=85000)

    core_half_width = 3

    mask_updated = np.empty(len(linenames), dtype=object)
    mask_withcore = np.empty(len(linenames), dtype=object)  # for Halpha
    wave_st_arr = np.empty(len(linenames), dtype=object)
    flux_st_arr = np.empty(len(linenames), dtype=object)

    mask_ind = np.round(np.linspace(0, len(garr) - 1, 5)).astype(int)

    fluxmin = 1
    for ii in range(len(linenames)):
        xrange = (linewaves[ii] - 30, linewaves[ii] + 30)
        sel = (wave > xrange[0]) & (wave < xrange[1])
        wave_st = wave[sel]
        flux_st = flux[sel]
        mask_st = np.asarray(mask_peakproc[sel], dtype=bool)

        mask_stcore = mask_st.copy()
        mask_stcore[(wave_st > linewaves[ii] - core_half_width) &
                    (wave_st < linewaves[ii] + core_half_width)] = True

        mask_updated[ii] = mask_st
        mask_withcore[ii] = mask_stcore
        wave_st_arr[ii] = wave_st
        flux_st_arr[ii] = flux_st

        if np.nanmin(flux_st) < fluxmin:
            fluxmin = np.nanmin(flux_st) - 0.05        

    # FIRST PLOT: SPECTRUM

    print("Showing the Balmer lines with the FASTWIND models that fulfil the ionization equilibrium...")

    rows = 1
    cols = len(Balmer_lineindex)
    fluxes_syn = []
    fig2, axs2 = plt.subplots(rows, cols, figsize=(7 * cols, 10))
 
    colors_cb = {"black": "#000000", "orange" : "#E69F00", "sky blue": "#56B4E9", "bluish green": "#009E73", 
            "yellow": "#F0E442", "blue": "#0072B2", "vermillion": "#D55E00", "reddish purple": "#CC79A7", 
            "light orange": "#ffd881", "light green" : "#00f6b3"}    

    color_list = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00'] 
 
    handles = [Line2D([0], [0], color="None", label=r"$\uparrow$ Models")]
    legend_onoff = "on"

 
    for i, ind in zip(range(len(Balmer_lineindex)), Balmer_lineindex):
        fluxes_syn.append([])
        wave_st = wave_st_arr[ind]
        flux_st = flux_st_arr[ind]
 
        axs2[i].plot(wave_st, flux_st, label=linenames[ind], c='b')

        count = 0
 
        for gi, Tarri, g_ind in zip(garr, Tarr, range(len(garr))):
            modeli = find_FW_model(FW_measure['Grid_name'][FW_measure["He"] == str(fixed_He)][0], 
                                    el_input,str(int(Tarri)), str(int(gi)),
                                    fixed_He, fixed_Q, fixed_b, fixed_abun)
            wavemodi, fluxmodi = broadened_FASTWIND_model(
                modeli, fixed_micro, linenames[ind], int(vsini_star),
                int(vmac_star), R=85000)
            fit_syn = np.interp(wave_st, wavemodi, fluxmodi)
            fluxes_syn[-1].append(fit_syn)

            if (gi > (np.amax(garr[mask_ind]) - 0.05)) and i > 0:
                mask_window = (fit_syn < 0.9999)
                mask_updated[ind] = mask_updated[ind] & mask_window
                mask_withcore[ind] = mask_withcore[ind] & mask_window

            if g_ind in mask_ind:
                count = count + 1
                axs2[i].plot(wavemodi, fluxmodi, alpha=0.5,
                                label='log g=%.1f, Teff=%i' % (float(gi) / 100, int(Tarri * 100)))
                if linenames[ind] != "HALPHA":

                    if legend_onoff != "off":
                        handles.append(Line2D(
                            [0], [0], color=color_list[count], linewidth=1,
                            linestyle='dashed',
                            label='log g=%.1f dex \nTeff=%i K' % (float(gi) / 100, int(Tarri * 100))))
                        legend_onoff = "saving"
 
        axs2[i].fill_between(wave_st, 0, 2, where=mask_updated[ind],
                              alpha=.25, color='cyan')
        axs2[i].axvline(linewaves[ind], color='red', ls='--', alpha=0.75)

        xrange = (linewaves[ind] - 25, linewaves[ind] + 25)
        axs2[i].set_xlim(*xrange)
        axs2[i].set_ylim(fluxmin, 1.5)
        axs2[i].legend(loc='upper right', fontsize="15", fancybox=True, framealpha=0.7)
        axs2[i].set_xlabel(r'$\lambda$ ($\AA$)')
 
        if legend_onoff == "saving":
            legend_onoff = "off"
 
    axs2[0].set_ylabel('Normalized flux')
    fig2.suptitle(
        '%s (%s) \n %s \n' % (star_name, spc_star, ref_file) +
        "He=%d Q=%d b=%d abSi=%.2f micro=%d" % (
            int(fixed_He), int(fixed_Q), int(fixed_b), float(fixed_abun), int(fixed_micro)),
        fontsize='24')
 
    if save == 'y':
        fig2.savefig(savepath + star_name + 'balmer.png')

    stop = input("Press enter to continue with the analysis.")
 
    plt.tight_layout()
    plt.show(block=False)

    #################################################################
    ################## SECOND PLOT: CHI^2  ##########################
    #################################################################    
    
    colors=['purple','orange','green','navy','red']

    print("Showing the Chi^2 fit of the observed vs model Balmer lines...")
    
    fig, ax = plt.subplots(2,1,figsize=(7,10)) # CHI^2 

    Balmer_lineindex=[index for index, i in enumerate(linenames) if (i[0]=='H') and (i[1:3]!='EI') and ((i[1]!='A'))]

    chi2min=np.empty(len(Balmer_lineindex))
    gcentr_i = np.empty(len(Balmer_lineindex))
    Balmer_lines = []
    chi2=np.empty(len(Balmer_lineindex), dtype=object)

    for i,ind in zip(range(len(Balmer_lineindex)),Balmer_lineindex):
        sigma=1
        wave_st=wave_st_arr[ind]
        flux_st=flux_st_arr[ind]
                            
        chi2i=np.zeros(len(garr),dtype='float64')
                
        for gi,Tarri,g_ind in zip(garr,Tarr,range(len(garr))):
            
            try: 
                fit_syn =fluxes_syn[ind][g_ind]
                
                chi2i[g_ind]=(1/len(wave_st[mask_updated[ind]]))*np.sum((np.array(fit_syn[mask_updated[ind]])-np.array(flux_st[mask_updated[ind]]))**2/sigma**2)  
                
            except:
                print('No model with the following parameters:')
                print(FW_measure['Grid_name'][FW_measure["He"] == fixed_He][0], 
                      el_input, str(int(Tarri)), str(int(gi)),
                      fixed_He, fixed_Q, fixed_b, fixed_abun)    
                chi2i[g_ind]=np.nan # DO NOT TAKE THESE POINTS INTO ACCOUNT IN THE ANALYSIS
        
        #CHI^2
        ax[0].plot(garr/100,chi2i,'.-',label=linenames[ind],color=colors[ind])

        chi2[i] = chi2i
        Balmer_lines.append(linenames[ind])
        gcentr_i[i] = int(garr[np.nanargmin(chi2i)]) 

        chi2min[i]=np.nanmin(chi2i) 

    ax[0].set_ylabel(r'$\chi^2_i$')
    
    g_median = np.nanmedian(gcentr_i)
    g_std = np.nanstd(gcentr_i)
    
    ax[0].axvline(g_median/100, color='grey',ls='--', alpha=0.75)
    ax[0].axvline(g_median/100-g_std/100, color='grey',ls='--', alpha=0.5)
    ax[0].axvline(g_median/100+g_std/100, color='grey',ls='--', alpha=0.5)

    g_lines=[]

    for i,ind in zip(range(len(Balmer_lineindex)),Balmer_lineindex):
        linei = linenames[ind]
        if (gcentr_i[i]<g_median+g_std) & (gcentr_i[i]>g_median-g_std):
            g_lines.append(linei)

    if len(g_lines)<2:
        for i,ind in zip(range(len(Balmer_lineindex)),Balmer_lineindex):
            linei = linenames[ind]
            if (gcentr_i[i]<g_median+11) & (gcentr_i[i]>g_median-11):
                g_lines.append(linei)     

    if len(g_lines)<2: # GO A BIT FURTHER
        for i,ind in zip(range(len(Balmer_lineindex)),Balmer_lineindex):
            linei = linenames[ind]
            if (gcentr_i[i]<g_median+21) & (gcentr_i[i]>g_median-21):
                g_lines.append(linei)                         

    padding_chi2 = 60 # TO SET THE LIMITS OF THE CHI^2

    mask=np.ones(len(garr),dtype=bool)
    for i in range(len(chi2)):
        mask=mask & (~np.isnan(chi2[i]))  # do not take into accound the nans

    for i,ind in zip(range(len(Balmer_lineindex)),Balmer_lineindex):

        if linenames[ind] in g_lines:
            chi2i = chi2[i]
            gmini = garr[np.nanargmin(chi2i)]
            order=2

            mask_chi2 = mask & (garr>gmini-padding_chi2) & (garr<gmini+padding_chi2)

            fit_i=np.polyfit(garr[mask_chi2],
                            chi2i[mask_chi2],order) 
            chi2new_i=np.polyval(fit_i,garr[mask_chi2])

            ax[0].plot(garr[mask_chi2]/100,chi2new_i,'--',
                       color=colors[ind], alpha=0.5, label=r'$\chi^2_i$') #! PINTAR TODO!!!        
            
        else: 
            chi2i = chi2[i] 

    Balmer_mask = np.array([True if Balmer_lines[i] in g_lines else False for i in range(len(Balmer_lines))], dtype=bool)
    chi2 = np.array([chi2[i] for i in range(len(chi2)) if Balmer_mask[i]], dtype=float)   # <-- filter chi2 too
    chi2min = chi2min[Balmer_mask]
    gcentr_i = gcentr_i[Balmer_mask]
    Balmer_lines = np.array(Balmer_lines)[Balmer_mask]

    glims = [np.nanmin(garr) - 5, np.nanmax(garr) + 5]   # from the earlier fix, not padding_chi2

    mask = np.ones(len(garr), dtype=bool)
    for c in chi2:
        mask = mask & (~np.isnan(c))

    garr = garr[mask]

    chi2_masked=np.empty(len(chi2), dtype=object)
    for i in range(len(chi2)):
        chi2_masked[i] = np.zeros(len(garr),dtype='float64')
        chi2_masked[i] = chi2[i][mask]
    
    chi2 = chi2_masked.copy()
    chi2final = np.zeros(len(garr))

    for g_ind in range(len(garr)):
        for i in range(len(chi2)):          # no more Balmer_lineindex needed here
            chi2final[g_ind] += chi2[i][g_ind] / chi2min[i]

    ax[1].set_title(" ".join(["%.1f"%chi2final[i] for i in range(len(chi2final))]), fontsize = "6")

    assert len(chi2) == len(chi2min) == len(Balmer_lines), "chi2 arrays out of sync with Balmer_lines"
    
    chi2max=0
    for i in range(len(chi2)):
        chi2i = chi2[i]
        gmini = garr[np.nanargmin(chi2i)]
        order=2
        mask_chi2 = (garr>gmini-padding_chi2) & (garr<gmini+padding_chi2)
        if np.nanmax(chi2i[mask_chi2])> chi2max:
            chi2max=np.nanmax(chi2i[(garr>glims[0]) & (garr<glims[1])])

    ax[0].legend(loc='upper right',fontsize="15",fancybox=True, framealpha=0.7) 
    ax[0].set_ylim(0,chi2max)
    
    from matplotlib import ticker
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    formatter.set_powerlimits((-1,1)) 
    ax[0].yaxis.set_major_formatter(formatter)   

    gnew=np.linspace(glims[0],glims[1],1000)

    order=2

    #ADDED ON 2026-02-10

    # “Confidence intervals were obtained from the χ² profile by interpolating χ²(θ), numerically locating the minimum, and determining the parameter values for which Δχ² = 1.”

    # FIRST INTERPOLATE CHI² 
    from scipy.interpolate import PchipInterpolator

    mask_nan = ~np.isnan(chi2final)
    y = chi2final.copy()[mask_nan]
    x = garr.copy()[mask_nan]

    f = PchipInterpolator(x, y)

    x_dense = np.linspace(x.min(), x.max(), 1000)
    y_dense = f(x_dense)
    
    # SECOND: Find the minimum of the interpolated χ²
    from scipy.optimize import minimize_scalar
    res = minimize_scalar(f, bounds=(x.min(), x.max()), method='bounded')
    x_min = res.x
    chi2_min = f(x_min) # This gives you a sub-grid minimum without assuming a parabola.
    
    gcentr = round(x_min)

    # THIRD: Solve Δχ² = 1

    def dchi2(x):
        return f(x) - chi2_min - 1.0
    
    ax[1].plot(garr/100,chi2final,'o',color='blue',label=r'Sampled $\chi^2$') 
    ax[1].plot(x_dense/100, y_dense, '-', color="blue", alpha=0.5) #,label='PCHIP interpolation')
    ax[1].axvline(x_min/100, linestyle='--', label='Interpolated minimum')
    ax[1].plot(gnew/100,(chi2_min+1)*np.ones(len(gnew)),'r--')

    handles.append(Line2D([0], [0],color=colors_cb["blue"], linestyle = "solid", 
            label=r'$\chi^2_{Total}$'))
    
    handles.append(Line2D([0], [0],color=colors_cb["vermillion"], linestyle = "dashed", 
            label=r'$1 \sigma$ '))
    
    # FIND ROOTS ON EACH SIDE

    from scipy.optimize import brentq

    edge_tol = 1  # how many grid points from the boundary still counts as "edge"
    argmin_discrete = np.nanargmin(chi2final)

    if argmin_discrete <= edge_tol:
        FLAG = "EDGE"
        x_right = brentq(dchi2, x_min, x.max())
        x_left = x_min - 10
    elif argmin_discrete >= len(garr) - 1 - edge_tol:
        FLAG = "EDGE"
        x_left = brentq(dchi2, x.min(), x_min)
        x_right = x_min + 10
    else:
        FLAG = "OK"
        try:
            x_left = brentq(dchi2, x.min(), x_min)
        except Exception:
            x_left = np.amin(garr)
        try:
            x_right = brentq(dchi2, x_min, x.max())
        except Exception:
            x_right = np.amax(garr)

    print("FLAG:", FLAG)    

    g1 = round(x_left)
    g2 = round(x_right)

    ax[1].legend(loc='upper right',fontsize="15",fancybox=True, framealpha=0.7) 
    ax[1].set_xlim(np.array(glims)/100)
    ax[0].set_xlim(np.array(glims)/100)

    ax[1].set_xlabel('log g') 
    ax[1].set_ylabel(r'$\chi^2_T$') 

    fig.suptitle('%s (%s) \n %s \n' %(star_name, spc_star,ref_file) + "He=%d Q=%d b=%d abSi=%.2f micro=%d"%(int(fixed_He),int(fixed_Q), int(fixed_b), float(fixed_abun), int(fixed_micro)), fontsize='24')

    if save=='y':
        fig.savefig(savepath+star_name+'chi2.png')
    plt.tight_layout() 
    plt.show(block=False)

    stop = input("Press enter to continue with the analysis.")


    # WE INTERPOLATE THE VALUE OF T THAT CORRESPONDS TO THE NEW G 
    garr_copy=list(table_calcTg['g'])
    Teff_copy=list(table_calcTg['Teff']) #Tarr, pero exacto 
    
    nearest_index=[min(range(len(garr_copy)),key=lambda i : abs(int(garr_copy[i])-int(gcentr)))]
    Terr=int(np.around(np.mean(np.array([table_calcTg['ATeff+'][nearest_index],
                    table_calcTg['ATeff-'][nearest_index]],dtype=int)),-2)) #! MEAN INSTEADO OF MAX? 
    if Terr<100: Terr = 100
     
    fitint=np.polyfit(garr_copy,Teff_copy,deg=1)
    
    Tcentr=np.polyval(fitint,gcentr) 

    gcentr = garr_copy[nearest_index[0]]
    Tcentr = Teff_copy[nearest_index[0]]

    try:
        T1=np.polyval(fitint,g1)    
    except:
        T1=np.amin(Tarr)
        print('WARNING: ERROR IN G IS TOO BIG')
    try:
        T2=np.polyval(fitint,g2) 
    except:
        T2=np.amax(Tarr)
        print('WARNING: ERROR IN G IS TOO BIG')
        
    ###################################################### 

    #################################################################
    ############ TERCER PLOT: HALPHA & HBETA VS Q  ##################
    #################################################################
    # 
    # UPDATE 2026 SEPT 23: FIJAMOS BETA (o no)

    Balmer_lineindex=[index for index, i in enumerate(linenames) if (i[0]=='H') and (i[1:3]!='EI')]

    garr_new=sorted(set(FW_measure['g']))

    nearest_g=garr_new[min(range(len(garr_new)),key=lambda i : abs(int(garr_new[i])-int(gcentr)))]
    Tarr_new=sorted(set(FW_measure['T'][FW_measure['g']==str(nearest_g)]))
    nearest_T=Tarr_new[min(range(len(Tarr_new)),key=lambda i : abs(int(Tarr_new[i])-int(Tcentr)/100))]
    
    mask=((FW_measure['T']==str(nearest_T)) & (FW_measure['g']==str(nearest_g)))

    
    fig1=plt.figure(figsize=(7,8),constrained_layout=True)
    gs1=fig1.add_gridspec(1,2)

    print("Showing the Chi^2 fit of Halpha and Hbeta (for determining the wind strength parameter logQ)...")

    count=0
    for bi in sorted(set(FW_measure['b'])):
        for Qi in sorted(set(FW_measure['Q'])):
            try:
                modeli=find_FW_model(FW_measure['Grid_name'][FW_measure["He"] == str(fixed_He)][0], 
                                     el_input, str('%d'%(int(nearest_T))), str(int(nearest_g)),
                              fixed_He, Qi, bi,  fixed_abun)
                
                wavemodi,fluxmodi=broadened_FASTWIND_model(modeli,fixed_micro, linenames[0] ,int(vsini_star),
                                                               int(vmac_star),R=85000)
        
                count = count + 1 
            except:
                print('Error in model with b=%s and Q=%s. Skipping...'%(bi,Qi))


    for c,ind in zip(range(len(Balmer_lineindex[0:2])),Balmer_lineindex[0:2]):        
        ax=fig1.add_subplot(gs1[0,c])
        ax.axvline(linewaves[ind], color='red',ls='--', alpha=0.75)

        ax.plot(wave_st_arr[ind],flux_st_arr[ind],label=linenames[ind],c='navy')
   
        for bi in sorted(set(FW_measure['b'])):  
            for Qi in sorted(set(FW_measure['Q'])):
                try:
                    modeli=find_FW_model(FW_measure['Grid_name'][FW_measure["He"] == str(fixed_He)][0], 
                                    el_input, str('%d'%(int(nearest_T))), str(int(nearest_g)),
                                    fixed_He, Qi, bi,  fixed_abun)
                
                    wavemodi,fluxmodi=broadened_FASTWIND_model(modeli,fixed_micro, linenames[ind] ,int(vsini_star),
                                                               int(vmac_star),R=85000)
                    ax.plot(wavemodi,fluxmodi,alpha=0.75,label='Q=%i,b=%i'%(int(Qi),int(bi)),linewidth=1)

                except:
                    print('Error in model with b=%s and Q=%s. Skipping...'%(bi,Qi))

        ax.set_xlim(np.amin(wave_st_arr[ind]),np.amax(wave_st_arr[ind]))       
        ax.set_ylim(-0.5,1.3)       
        ax.set_xlabel(r'$\lambda$ ($\AA$)')
        ax.legend(loc='lower right',fontsize="15",fancybox=True, framealpha=0.7)
        if c<1:
            ax.set_ylabel('Normalized flux')

    
    fig1.suptitle('%s (%s) \n %s \n' %(star_name, spc_star,ref_file) + "log g=%.1f dex, Teff=%i K, He=%d Q=%d b=%d abSi=%.2f micro=%d"%(int(nearest_g)/100,int(nearest_T)*100,int(fixed_He),int(fixed_Q), int(fixed_b), float(fixed_abun), int(fixed_micro)), fontsize='22') #+ "Q=" + fixed_Q +  "b=" + fixed_b
    
    if save=='y':
        fig1.savefig(savepath+star_name+'_Q.png')   

    stop = input("Press enter to continue with the analysis.")
 
    
    plt.show(block=False) 
    
    if show == "y": print('\n SPECTRAL SYNTHESIS FOR Q ')
    chi2min_Q=np.empty(len(Balmer_lineindex[0:2]))
    chi2_Q=[]
    for c,ind in zip(range(len(Balmer_lineindex[0:2])),Balmer_lineindex[0:2]):        
        sigma=1
        wave_st=wave_st_arr[ind]
        flux_st=flux_st_arr[ind]      
        
        chi2i_Q=np.zeros(len(sorted(set(FW_measure['b'])))*len(sorted(set(FW_measure['Q']))),dtype='float64')
        bb=np.empty(len(sorted(set(FW_measure['b'])))*len(sorted(set(FW_measure['Q']))),dtype='object')
        QQ=np.empty(len(sorted(set(FW_measure['b'])))*len(sorted(set(FW_measure['Q']))),dtype='object')

        count_Q=0
        for bi in sorted(set(FW_measure['b'])):
            for Qi in sorted(set(FW_measure['Q'])):
                try:
                    modeli=find_FW_model(FW_measure['Grid_name'][FW_measure["He"] == str(fixed_He)][0], 
                                        el_input, str('%d'%(int(nearest_T))), str(int(nearest_g)),
                                    fixed_He, Qi, fixed_b,  fixed_abun)
                    wavemodi,fluxmodi=broadened_FASTWIND_model(modeli,fixed_micro, linenames[ind] ,int(vsini_star),
                                                                    int(vmac_star),R=85000)
                    fit_syn=np.interp(wave_st, wavemodi, fluxmodi)

                    chi2i_Q[count_Q]=(1/len(wave_st[mask_withcore[ind]]))*np.sum((np.array(fit_syn[mask_withcore[ind]])-np.array(flux_st[mask_withcore[ind]]))**2/sigma**2)*1e4
                    bb[count_Q]=bi
                    QQ[count_Q]=Qi
                    count_Q=count_Q+1
                except:
                    chi2i_Q[count_Q]=np.nan
                    bb[count_Q]=bi
                    QQ[count_Q]=Qi
                    count_Q=count_Q+1
                    print('No model with b=%s and Q=%s. Skipping...'%(fixed_b,Qi))

        chi2_Q.append(chi2i_Q)
        chi2min_Q[ind]=np.nanmin(chi2i_Q)   
        
        chi2i_argsort=np.argsort(chi2i_Q)
        chi2i_sort=chi2i_Q[chi2i_argsort]

        # if show == "y": 
            # print('\n %s'%(linenames[ind]))
            # print('best chi2:',chi2i_sort[0:5])
            # print('best beta',bb[chi2i_argsort][0:5])
            # print('best Q', QQ[chi2i_argsort][0:5])
            
            # print(r'\n %s: b=%s, Q=%s, $\chi^2_{min}$=%f \n'%
            #     (linenames[ind],bb[np.nanargmin(chi2i_Q)],QQ[np.nanargmin(chi2i_Q)],np.nanmin(chi2i_Q)  ))

        # USAMOS SOLO HALPHA PARA EL VIENTO        
        # if linenames[ind] == "HALPHA":
        #     best_QQ = QQ[chi2i_argsort][0]
        #     best_bb= bb[chi2i_argsort][0]

    # USAMOS HALPHA Y HBETA PARA EL VIENTO 

    chi2final = np.zeros(len(bb))
    for ind, bi, Qi in zip(range(len(bb)), bb, QQ ):
        for i in range(len(chi2min_Q)):
            chi2final[ind] += chi2_Q[i][ind] / chi2min_Q[i]

    chi2i_argsort=np.argsort(chi2final)
    chi2i_sort= np.sort(chi2final)

    # print('\n total chi2:')
    # print('best chi2:',chi2i_sort[0:5])
    # print('best beta',bb[chi2i_argsort][0:5])
    # print('best Q', QQ[chi2i_argsort][0:5])    

    best_QQ = QQ[chi2i_argsort][0]
    best_bb= bb[chi2i_argsort][0]    

    if show == "y": print("beta = %d, log Q = %d" %(int(best_bb), int(best_QQ)))

    #################################################################
    ########## CUARTO PLOT: LINES H HE WITH BEST FIT  ###############
    #################################################################
    
    rows0=2
    cols0=len(He_lineindex)*len(Balmer_lineindex)
        
    fig0=plt.figure(figsize=(20,10),constrained_layout=True)    
    gs0=fig0.add_gridspec(rows0,cols0)
    
    #Primero dibujamos las lineas de Balmer
    
    row0=0
    start=0
    step=int(cols0/len(Balmer_lineindex))
    fluxmin=1
    for c,ind in zip(range(len(Balmer_lineindex)),Balmer_lineindex):
        if np.nanmin(flux_st_arr[ind])<fluxmin:
            fluxmin=np.nanmin(flux_st_arr[ind])-0.05
    
    for c,ind in zip(range(len(Balmer_lineindex)),Balmer_lineindex):        
        ax=fig0.add_subplot(gs0[0,int(start):int(start+step)])
        ax.axvline(linewaves[ind], color='red',ls='--', alpha=0.75)

        ax.plot(wave_st_arr[ind],flux_st_arr[ind],label=linenames[ind],c='navy')
        
        nearest_g=garr_new[min(range(len(garr_new)),key=lambda i : abs(int(garr_new[i])-int(gcentr)))]
        Tarr_new=sorted(set(FW_measure['T'][FW_measure['g']==str(nearest_g)]))
        nearest_T=Tarr_new[min(range(len(Tarr_new)),key=lambda i : abs(int(Tarr_new[i])-int(Tcentr)/100))]
        mask=((FW_measure['T']==nearest_T) & (FW_measure['g']==nearest_g))

        modeli=find_FW_model(FW_measure['Grid_name'][FW_measure["He"] == str(fixed_He)][0], 
                             el_input, str('%d'%(int(nearest_T))), str(int(nearest_g)),
                      fixed_He, fixed_Q, fixed_b,  fixed_abun)
        
        wavemodi,fluxmodi=broadened_FASTWIND_model(modeli,fixed_micro, linenames[ind] ,int(vsini_star),
                                                       int(vmac_star),R=85000) 
                        
        ax.plot(wavemodi,fluxmodi,alpha=0.75,label='Q=%i, b=%i'%(int(fixed_Q),int(fixed_b)),c='green')
        
        ax.set_xlim(np.amin(wave_st_arr[ind]),np.amax(wave_st_arr[ind]))       
        ax.set_ylim(fluxmin,1.3)       
        ax.set_xlabel(r'$\lambda$ ($\AA$)')
        ax.legend(loc='upper right',fontsize="15",fancybox=True, framealpha=0.7)
        if start<1:
            ax.set_ylabel('Normalized flux')
        
        start=start+step
          
    #Ahora dibujamos las lineas de He
    
    row0=1
    start=0
    step=int(cols0/len(He_lineindex))
    fluxmin=1
    for c,ind in zip(range(len(He_lineindex)),He_lineindex):
        if np.nanmin(flux_st_arr[ind])<fluxmin:
            fluxmin=np.nanmin(flux_st_arr[ind])-0.05    
    for c,ind in zip(range(len(He_lineindex)),He_lineindex):
        ax=fig0.add_subplot(gs0[1,int(start):int(start+step)])
        ax.axvline(linewaves[ind], color='red',ls='--', alpha=0.75)

        ax.plot(wave_st_arr[ind],flux_st_arr[ind],label=linenames[ind],c='navy')
        
        modeli=find_FW_model(FW_measure['Grid_name'][FW_measure["He"] == str(fixed_He)][0], 
                             el_input, str('%d'%(int(nearest_T))), str(int(nearest_g)),
                      fixed_He, fixed_Q, fixed_b,  fixed_abun)
        
        wavemodi,fluxmodi=broadened_FASTWIND_model(modeli,fixed_micro, linenames[ind] ,int(vsini_star),
                                                       int(vmac_star),R=85000)
                       
        ax.plot(wavemodi,fluxmodi,alpha=0.75,label='He=%i'%(int(fixed_He)),c='green')
        
        ax.set_xlim(np.amin(wave_st_arr[ind]),np.amax(wave_st_arr[ind]))
        ax.set_ylim(fluxmin,1.3)       
        ax.legend(loc='upper right',fontsize="15",fancybox=True, framealpha=0.7)
        ax.set_xlabel(r'$\lambda$ ($\AA$)')
        ax.legend(loc='upper right',fontsize="15",fancybox=True, framealpha=0.7)
        if start<1:
            ax.set_ylabel('Normalized flux')

        start=start+step    
        # ax.set_title( "He= " + fixed_He + "Q=" + fixed_Q + "b=" + fixed_b + "abSi=" + fixed_abun + "micro" + fixed_micro)  
    # fig0.suptitle('%s %s - log g=%.1f dex, Teff=%i K' %(star_name,spc_star,int(nearest_g)/100,int(nearest_T)*100), fontsize=25)

    # fig0.suptitle('%s (%s) \n %s \n' %(star_name, spc_star,ref_file) + "log g=%.1f dex, Teff=%i K"%(int(nearest_g)/100,int(nearest_T)*100)   + "Q=" + fixed_Q +  "b=" + fixed_b + "He= " + fixed_He  + "abSi=" + fixed_abun + "micro" + fixed_micro, fontsize='22') #+ "Q=" + fixed_Q +  "b=" + fixed_b
    fig0.suptitle('%s (%s) \n %s \n' %(star_name, spc_star,ref_file) + "log g=%.1f dex, Teff=%i K, He=%d Q=%d b=%d abSi=%.2f micro=%d"%(int(nearest_g)/100,int(nearest_T)*100,int(fixed_He),int(fixed_Q), int(fixed_b), float(fixed_abun), int(fixed_micro)), fontsize='22') #+ "Q=" + fixed_Q +  "b=" + fixed_b

    if save=='y':
        fig0.savefig(savepath+star_name+'_HHe.png')    
    
    plt.show(block=False) 

    ###########################################################################
    
    ###########################################################################
    ################### TABLA FINAL CON LOS RESULTADOS  #######################
    ###########################################################################

    tab_Tg=Table([[str(int(Tcentr)),str(int(Tcentr+Terr)),str(int(Tcentr-Terr)),str(int(T1)),str(int(T2))],
                    [str(int(gcentr)),str(int(gcentr)),str(int(gcentr)),str(int(g1)),str(int(g2))],
                    ['10','10','10','10','10'],
                    np.array(int(best_QQ)*np.ones(5),dtype='int').astype('str'),
                    np.array(int(best_bb)*np.ones(5),dtype='int').astype('str'),
                    [fixed_micro,fixed_micro,fixed_micro,fixed_micro,fixed_micro],
                    [fixed_abun,fixed_abun,fixed_abun,fixed_abun,fixed_abun],                    
                    ],
                    names=['T','g','He','Q','b','micro','abun'],
                    dtype=['S10','S10','S10','S10','S10','S10','S10'])
    
    # T_output.write(savepath+star_name+'_Tg_final.ecsv', format='ascii.ecsv', overwrite=True)

    ###########################################################################    
    
    tab_inputparameters=Table([[ref_file], [str(int(Tcentr))],[str(int(Terr))],[str(int(T1))],[str(int(T2))],[str(int(gcentr))],[str(int(g1))],[str(int(g2))],
                        [fixed_He],[best_QQ],[best_bb],[fixed_micro],[fixed_abun]],
                       names=["Ref_file",'T','Terr','T+','T-','g','g+','g-','He','Q','b','micro','abun'],
                       dtype=["S200",'S10','S10','S10','S10','S10','S10','S10','S10','S10','S10','S10','S10'])    
    
    parameterstab2=findtable(star_name+'_input_parameters.ecsv',format='ecsv',path=savepath,delimiter='\s')
    if parameterstab2==None:
        if save=='y':
            tab_inputparameters.write(savepath+star_name+'_input_parameters.ecsv', format='ascii.ecsv', overwrite=True)
    else:
        try:
            tab_inputparameters=vstack([parameterstab2,tab_inputparameters])
        except: 
            parameterstab2.pprint_all()
            tab_inputparameters.pprint_all()
        if save=='y':
            tab_inputparameters.write(savepath+star_name+'_input_parameters.ecsv', format='ascii.ecsv', overwrite=True)
    # print(T_output_new)

    g_lines = "-".join(g_lines)

    
    return tab_Tg,tab_inputparameters, best_QQ, best_bb, g_lines, FLAG 


'''
##############################################################################
######################     CURVE OF GROWTH     ###############################
##############################################################################
'''

def tables_curve_of_growth(ref_file, el_input, spc_star, star_measure, syn_measure, #table,
                           inputparameters,linetable, micro_selection, savepath, micro_center=None,
                           linesvec=['SiIII4552','SiIII4567','SiIII4574'],mult_use="all",
                           do_pdf=True,error=5,ewcut=5, tol=50, point='0',save='y', show="y"):
    
    '''
    Programa que genera las tablas necesarias para el método de la curva de crecimiento.
    
    Parameters:
    ----------  
    star_name: str
        El nombre de la estrella.  
        
    spc_star: str
        Tipo espectral de la estrella.
        
    star_measure: str
        Nombre de la tabla en la que se encuentran las medidas de EW de la estrella.
        
    syn_measure: str
        Nombre de la tabla en la que se encuentran las medidas de EW de los espectros sintéticos. 
        
    table: Table -> EN PROCESO DE ELIMINAR -> PROBAMOS A HACER EL CÁLCULO CON UN ÚNICO T, G
        Tabla con la información de la T y g calculadas con sus respectivos errores, generada por
        la función best_Tg. 
        
    linetable: str
        Enter the wavelenght(s) of the line(s) to show, either in a coma-separated
        string, or in a .txt/.lst file containing the lines.
        
    el_input: str
        Nombre del elemento del que quieres calcular la abundancia.
        Opciones: Si, C, N, O
         
    micro_selection: str
        Tres opciones:
        'up': crea las tablas teniendo en cuenta los valores de microturbulencia más altos del grid.
        'down': crea las tablas teniendo en cuenta los valores de microturbulencia más bajos del grid.
        'mid': crea las tablas teniendo en cuenta los valores de microturbulencia centrales del grid.
            
    savepath: str
        Path donde se guardan las figuras y tablas generadas. 
   
    do_pdf: Bool
        If True, the function generates and saves a pdf with the graphs for the central T, g. 
        Default is True. 
        
    error: int
        Error considered for the EW measurement. 
        Default is 10.
        [NOTE: this would need to be updated to consider the actual uncertainty in the
         measurement of the EW.]
        
    ewcut: int
        EW threshold value for a line to be considered as detected. 
        Default is 10.
    
    tol: int
        Sets the tolerance [km/s] to shifting the spectrum in order to fit the lines 
        after being corrected from radial velocity.
        Default is 50.
    
    point: str
        If '0', does the analysis for the central point, T0,g0
        If '1', does the analysis for the upper uncertainty point, T1, g1
        If '2', does the analysis for the lower uncertainty point, T2, g2

    Returns:
    ----------   
    final_tables: Table
        Tabla con la información de la EW y la abundancia correspondiente a cada T, g, línea, y microturbulencia. 
        Esta tabla también se guarda en los archivos star_name+'_curveofgrowthinput_'+str(ii)+'.ecsv', 
        donde ii son los distintos puntos para los distintos valores de T, g teniendo en cuenta los errores. 
        
    También se genera un pdf mostrando la EW y la abundancia de cada línea y cada microturbulencia para la T, g centrales.
    '''
    
    while el_input not in ['Si','C','N','O']:
        print("Element not supported.")
        el_input = input("Choose element from 'Si','C','N','O':  ")
        
    # syn_measure=findtable(syn_measure)
    star_measure = findtable(star_measure)
    star_measure = star_measure[star_measure['Ref_file']==ref_file][0]
    star_name=star_measure["ID"]
        
    fulllinestab = findtable(linetable,format='csv',path=savepath)
    linestab = fulllinestab ['WL']
    elementstab = fulllinestab['label']
    nn = fulllinestab[ref_file]
    flagtab = fulllinestab['flag']
    
    nn = np.array(nn,dtype=float)
    lines = np.array(linestab,dtype=float)[nn<1] #0 means seen in spectra + good line. 
    elements = np.array(elementstab,dtype=object)[nn<1]
    
    lines_full = np.array(linestab,dtype=float)[nn!=1] # 1 means not seen in spectra, we consider everything else
    elements_full = np.array(elementstab,dtype=object)[nn!=1] 
    flag_full = np.array(flagtab,dtype=object)[nn!=1] 
    nn_full = np.array(nn,dtype=object)[nn!=1] 

    if ref_file in fulllinestab.colnames:
        if "all" in linesvec:
            linesvec = list(fulllinestab['label'][fulllinestab[ref_file]<1])
        else:
            linesvec = [i for i in linesvec if i in list(fulllinestab['label'][fulllinestab[ref_file]<1])]
    
    FW_EW_measuretab=findtable(syn_measure)

    if isinstance(inputparameters,str):
        inputparameters=findtable(inputparameters,format='ecsv',path=savepath,delimiter='\s')
    inputparameters = inputparameters [inputparameters["Ref_file"]==ref_file]

    if point=='0':

        # CALCULATE THE MICRO AND ABUNDANCES FOR A CERTAIN T AND G

        fixed_T=inputparameters['T'][-1]
        fixed_g=inputparameters['g'][-1]
        fixed_He=inputparameters['He'][-1]
        fixed_Q=inputparameters['Q'][-1]
        fixed_b=inputparameters['b'][-1]
        
        garrFW=sorted(set(FW_EW_measuretab['g']))
        fixed_g=garrFW[min(range(len(garrFW)),key=lambda i : abs(int(garrFW[i])-int(fixed_g)))] #NEAREST
        
    elif point=='1':

        # CALCULATE THE ABUNDANCES WITH FIXED MICRO TO ACCOUNT FOR THE UPPER UNCERTAINTY IN T, G

        fixed_T=inputparameters['T+'][-1]
        fixed_g=inputparameters['g+'][-1]
        fixed_He=inputparameters['He'][-1]
        fixed_Q=inputparameters['Q'][-1]
        fixed_b=inputparameters['b'][-1]

        garrFW=sorted(set(FW_EW_measuretab['g']))
        fixed_g=garrFW[min(range(len(garrFW)),key=lambda i : abs(int(garrFW[i])-int(fixed_g)))] #NEAREST


    elif point=='2':

        # CALCULATE THE ABUNDANCES WITH FIXED MICRO TO ACCOUNT FOR THE LOWER UNCERTAINTY IN T, G

        fixed_T=inputparameters['T-'][-1]
        fixed_g=inputparameters['g-'][-1]
        fixed_He=inputparameters['He'][-1]
        fixed_Q=inputparameters['Q'][-1]
        fixed_b=inputparameters['b'][-1]
        
        garrFW=sorted(set(FW_EW_measuretab['g']))
        fixed_g=garrFW[min(range(len(garrFW)),key=lambda i : abs(int(garrFW[i])-int(fixed_g)))] #NEAREST
    
    elif point=='Q': # we take into account the first parameter determination of T, g, He, with low winds!
        fixed_T=inputparameters['T'][-1] 
        fixed_g=inputparameters['g'][-1]
        fixed_He=inputparameters['He'][-1]
        fixed_Q="1400"
        fixed_b="10"
        
        garrFW=sorted(set(FW_EW_measuretab['g']))
        fixed_g=garrFW[min(range(len(garrFW)),key=lambda i : abs(int(garrFW[i])-int(fixed_g)))] #NEAREST
    
    TarrFW=sorted(set(FW_EW_measuretab[FW_EW_measuretab["g"] == "%d"%int(fixed_g)]['T']))
    fixed_T=TarrFW[min(range(len(TarrFW)),key=lambda i : abs(int(TarrFW[i])-int(fixed_T)/100))] #NEAREST

    mask=((FW_EW_measuretab['T']=="%d"%(int(fixed_T)) ) & (FW_EW_measuretab['g']=="%d"%int(fixed_g)) &
                      (FW_EW_measuretab['He']=="%d"%int(fixed_He)) & (FW_EW_measuretab['Q']=="%d"%int(fixed_Q)) & 
                      (FW_EW_measuretab['b']=="%d"%int(fixed_b))
                      )
    
    synth_measure=FW_EW_measuretab[mask]

    micro=sorted(set(synth_measure['micro'])) #! micro = from FW_EW_measure -> ok just micro

    #PROBAMOS CON TODAS LAS MICROS.
    K=5
    if micro_selection =='down':
        micro=micro[0:K]
    elif micro_selection =='up':
        micro=micro[-K:]
    elif micro_selection =='mid':
        middle = len(micro) // 2
        micro=[micro[i] for i in range(middle - K//2, middle + (K//2) +1)]
    elif micro_selection=='full':
        micro=micro
    elif micro_selection=='center':

        ind_nearest=min(range(len(micro)),key=lambda i : abs(int(micro[i])-float(micro_center)))
        if ind_nearest<2:
            micro_new=np.array(micro,dtype=str)[0:5]
        elif ind_nearest>(len(micro)-3):
            micro_new=np.array(micro,dtype=str)[-5:]
        else:
            micro_new=np.array(micro,dtype=str)[ind_nearest-2:ind_nearest+3]  

        micro=list(sorted(micro_new))

    # MICRO = ARRAY OF FIVE NUMBERS E.G. [ '009', '011', '015', '020', '025']

    #################################################################################
    # MONTE CARLO:
    linesvec_final = []
    EWvec=[]
    abvec=[[],[],[],[],[]] #one for each micro
    #################################################################################
    
    final_table=Table(names=['Lambda','Line','Flag','EW','AEW']+ 
                  [j+i for i in [str(int(k)) for k in micro] for j in ['ab'+el_input+'m','Aab'+el_input+'m_']],
                  dtype=['float64','S10','int64','float64','float64']+['float64','float64']*len(micro))
    
    #Primero definimos la tabla con las medidas de EW
    measurements=Table(names=['Lambda','Line','Flag','EW0','AEW0'],dtype=['float64','S10','int64','float64','float64']) 

    for i,line,element in zip(range(len(lines_full)),lines_full,elements_full):
        flag=False #if True, no tenemos en cuenta esta línea. 
        if "RV0" in star_measure.colnames:
            if ((np.ma.getmask(star_measure['EW_'+element])) or (star_measure['EW_'+element]<ewcut)
                or (np.abs(star_measure['RV0']-star_measure['RV_'+element])>0.5*tol)):
                if show=="y": print(element + " does not fit star quality cuts.")
                flag=True

                continue
            else:
                EW0=star_measure['EW_'+element]
                #AEW0=star_measure['EW_'+element]*error
                AEW0=error
                flagline=flag_full[i]
                nnline = nn_full[i]
                row=Table([[line],[element],[nnline],[EW0],[AEW0]],names=['Lambda','Line','Flag','EW0','AEW0']) 
                measurements=vstack([measurements,row])             
        else:
            if ((np.ma.getmask(star_measure['EW_'+element])) or (star_measure['EW_'+element]<ewcut)):
                if show=="y": print(element + " does not fit star quality cuts.")
                flag=True
                continue         
            else:
                EW0=star_measure['EW_'+element]
                #AEW0=star_measure['EW_'+element]*error
                AEW0=error
                flagline=flag_full[i]
                nnline = nn_full[i]
                row=Table([[line],[element],[nnline],[EW0],[AEW0]],names=['Lambda','Line','Flag','EW0','AEW0']) 
                measurements=vstack([measurements,row]) 

    models_list=[]

    for i,line,element in zip(range(len(measurements)),measurements['Lambda'],measurements['Line']):
        flag=False #if True, no tenemos en cuenta esta línea. 
        models=Table(names=['Lambda','Line','m','AO','EW'],
                 dtype=['float64','S10','int64','float64','float64'])
        for jj in micro:
            Si=sorted(set(synth_measure['abun_'+el_input][synth_measure['micro']==str(jj)]))
            for kk in Si: 

                EW = synth_measure['EW_'+element][(synth_measure['micro']==str(jj)) &
                                                (synth_measure['abun_'+el_input]==str(kk))]
                if EW<=0:
                    if show=="y": print(element + " does not fit synth quality cuts HERE !!!.")
                    continue
                row = Table([[line],[element],[int(jj)],[float(kk)],EW], names=['Lambda','Line','m','AO','EW'])
                models = vstack([models,row])

        AO=np.empty(len(micro*2),dtype=float) 
        EW0=measurements['EW0'][measurements['Line']==element][0]
        AEW0=measurements['AEW0'][measurements['Line']==element][0]
        nnline=measurements['Flag'][measurements['Line']==element][0]
        n=-2

        ok = True

        #################################################################################
        # MONTE CARLO TEST:    
        montecarlosize=50000
        montecarlosize=50 #! so it's shorter, for tests... 
        if element in linesvec:
            EWrand=np.random.normal(EW0,AEW0,montecarlosize)
            EWvec.append(EWrand)
            linesvec_final.append(element)
        #################################################################################


        for l,micr in zip(range(len(micro)),micro):
            n=n+2
            mask=(np.array(models['m'],dtype=str)==str(int(micr)))

            labeli=r'$\xi {}$'.format(int(micr))

            if len(models["EW"][mask])>0:

                fit=np.polyfit( models['EW'][mask],models['AO'][mask], 1)
                AO_EW=np.polyval(fit,EW0) 
                AO_EW_inf=np.polyval(fit,EW0-AEW0) 
                AO_EW_sup=np.polyval(fit,EW0+AEW0)
                AAO_EW_sup=np.abs(AO_EW_sup-AO_EW)
                AAO_EW_inf=np.abs(AO_EW_inf-AO_EW)
                #el error en algún momento tendría que pasar a ser asimétrico
                AO_EW_err=np.abs(max([AO_EW_sup,AO_EW_inf])-AO_EW)
                AO[n]=AO_EW
                AO[n+1]=AO_EW_err

                ################elementif oi#################################################################
                if element in linesvec:
                    abvecelem=np.polyval(fit,EWrand)
                    abvec[l].append(abvecelem)
                ################################################################################
            else:
                ok = False

        if not ok or np.any(np.isnan(AO)):
            print("can't fit without info for all micros -> we won't use this line: "+element)
            continue
        else:
            row_line=Table([[line],[element],[nnline],[EW0],[AEW0]]+[[j] for j in AO],
                           names=['Lambda','Line','Flag','EW','AEW']+ 
                  [j+i for i in [str(int(k)) for k in micro] for j in ['ab'+el_input+'m','Aab'+el_input+'m_']],
                  dtype=['float64','S10','int64','float64','float64']+['float64','float64']*len(micro))
            final_table=vstack([final_table,row_line])

            models_list.append(models)
    
    #final formatting:
    final_table['EW'].info.format='7.2f'
    final_table['AEW'].info.format='7.2f'

    ii=int(point)     
    if ii<1:
    #FIGURE 01     
        if do_pdf==True:
            from matplotlib.backends.backend_pdf import PdfPages

            if not os.path.exists(savepath):
                os.makedirs(savepath)
    
            pdf_fitting = PdfPages(savepath + '/abundances_%s_%s.pdf' % (star_measure['ID'],el_input))
                    
            nrows, ncols = even_plot(len(final_table['Line']))

            greenlines = final_table['Line'][final_table['Flag'] <1 ]

            count = 0
            # Avoid showing the spectrum for each line when is already saved in pdf
            show_plot = False
            ppp=25 #plots per page
            for i,line,element in zip(range(len(final_table)),final_table['Lambda'],final_table['Line']):
                models=models_list[i]

                if (i % ppp) == 0:
                    if len(final_table['Lambda']) <= ppp:
                        fig,ax = plt.subplots(nrows, ncols, figsize=(16,16))
                    else:
                        fig,ax = plt.subplots(5, 5, figsize=(16,16))
                    axs =  ax.flatten()

                if len(models)<1:
                    continue
                else:
                    EW0=measurements['EW0'][measurements['Line']==element]
                    AEW0=measurements['AEW0'][measurements['Line']==element]
   
                    if element in greenlines:
                        
                        axs[i % ppp].set_title('%s' % element, c='black')

                        for l,micr in zip(range(len(micro)),micro):
                            labeli=r'$\xi {} $'.format(int(micr))
                            mask=(np.array(models['m'],dtype=str)==str(int(micr)))

                        count = count + 1
                    else:
                        axs[i % ppp].set_title('%s' % element, c='grey')

                    axs[i % ppp].set_ylabel(r'$EW$ [m$\AA$]') 
                    axs[i % ppp].set_xlabel(r'$\epsilon$ (Si)') 
                    for l,micr in zip(range(len(micro)),micro):
                        labeli=r'$\xi {} $'.format(int(micr))
                        mask=(np.array(models['m'],dtype=str)==str(int(micr)))

                        axs[i % ppp].plot(models['AO'][mask],models['EW'][mask], '.-')
                        axs[i % ppp].text(max(models['AO'][mask]),max(models['EW'][mask]),labeli)

                    axs[i % ppp].plot(np.linspace(min(models['AO']),max(models['AO']), 100), EW0*np.ones(100), '--k')
                    axs[i % ppp].plot(np.linspace(min(models['AO']),max(models['AO']), 100), EW0+AEW0*np.ones(100), color='gray', linestyle='dashed')
                    axs[i % ppp].plot(np.linspace(min(models['AO']),max(models['AO']), 100), EW0-AEW0*np.ones(100), color='gray', linestyle='dashed')
                    axs[i % ppp].set_xlim(min(models['AO'])-0.15,max(models['AO'])+0.15)  
                    axs[i % ppp].tick_params(direction='in', top='on', right='on')

                if i % ppp == ppp-1 or (i == len(final_table['Lambda'])-1):
                    fig.suptitle('%s %s %s Abundance Determination' %(star_measure['ID'],spc_star,el_input), fontsize=20)
                    
                    # Delete the unused axes
                    [fig.delaxes(axs[i]) for i in np.arange(i % ppp + 1, len(axs))]
                        
                    fig.tight_layout()
                    fig.subplots_adjust(wspace=0.6, hspace=0.4)
                    if el_input in ["Si"]:
                        plt.savefig(savepath + '/abundances_%s_%s.png' % (star_name, el_input))
                    plt.show(block=False)
                    
                    pdf_fitting.savefig(fig)
                    
            pdf_fitting.close()

    # final_table.write(maindir + 'tables/tests_curvadecrecimiento/measure_input_test_file.fits', format='fits', overwrite=True)         
    final_table.write(savepath + star_name+'_curveofgrowthinput_'+el_input+'_'+str(ii)+'.ecsv',
                        format='ascii.ecsv',overwrite=True)
    if show == "y":
        print('\nTable saved in '+ savepath + star_name+'_curveofgrowthinput_'+el_input+'_'+point+'.ecsv')
        print('Figure saved in '+ savepath + 'abundances_%s.pdf' % star_name)
        print('Line-list in '+ savepath+ linetable)

    if len(linesvec_final)<3:
        if el_input == 'Si':
            print('CANNOT MEASURE ALL SI III TRIPLET: CANNOT DO CURVE OF GROWTH')
            print(linesvec_final)
            print(EWvec)
            print(abvec)
        else:
            print('Cannot measure at least three lines: cannot do curve of growth. ')
        raise Exception("Cannot measure Si III triplet: cannot do curve of growth") 
    return (final_table, linesvec_final, EWvec, abvec)

#%%
def abundances_growthcurve(star_name, ref_file, spc_star, input_table, inputparameters,el_input,linetable,FW_EW_measuretab,parameterstable,
                           savepath, microreduced,linesvec,EWvec,abvec,mult_use="all", mult_highlight=None,
                           do_pdf=True,point='0',save='y', show="y"): #! usar linesvec para marcar con qué líneas calculamos la micro!

    '''
    Programa que calcula la abundancia y microturubulencia de la estrella usando
    el método de la curva de crecimiento.
    
    Parameters:
    ---------- 
    
    star_name: str
        El nombre de la estrella.
        
    spc_star: str
        Tipo espectral de la estrella.
        
     input_table: Table
        Tabla creada por la función abundances_input_file, que contiene,
        para cada línea con su medida de anchura equivalente de una estrella dada,
        la abundancia y el error de un elemento que corresponde a cada valor de
        microturbulencia.    
        
    table_Tg: Table
        Tabla con la información de la T y g calculadas con sus respectivos errores, 
        generada por la función best_Tg.    
        
    el_input: str
        Nombre del elemento del que quieres calcular la abundancia.
        Opciones: Si, (C, N, O)
        Default is Si
        
    m: int 
        constante para empezar a medir las columnas de abundancias según la microturbulencia.
        (Para los ficheros usados actualmente la constante debe ser m = 4,
         ya que antes de llegar a las columnas de abundancias hay cuatro columnas.)
        
        [!]: igual se puede eliminar y poner debajo --> mirar a ver si aparece en más funciones
        
    micro_fixed: int, optional. 
        micro_fixed = None (Default): 
            el programa busca la mejor microturbulencia y calcula la abundancia a partir de este valor
        micro_fixed = int
            el programa calcula la abundancia a partir de un valor de microturbulencia fijo
            
    microreduced: str, optional
        Can be 'yes', 'montecaro' or 'fixed'
            
    Tsel: str
        Opciones: '0','1','2','3','4'. 
        Representa el índice de la tabla creada por abundances_input_file. Cada índice
        tiene información de un valor de T, g diferente, y es necesario para hacer el 
        cálculo de las incertidumbres.
    
    savepath: str
        Path donde se guardan las figuras y tablas generadas.               
    

    Returns:
    ----------   
        
        final_table: Table
            Tabla con la información de los parámetros espectroscópicos calculados con sus 
            respectivas incertidumbres: temperatura, gravedad, microturbulencia, y abundancias. 

    '''
    m=5 #?    

    colors = {"black": "#000000", "orange" : "#E69F00", "sky blue": "#56B4E9", "bluish green": "#009E73", 
              "yellow": "#F0E442", "blue": "#0072B2", "vermillion": "#D55E00", "reddish purple": "#CC79A7", 
              "light orange": "#ffd881", "light green" : "#00f6b3"}    

    if isinstance(input_table, str):
        input_table=findtable(input_table)

    input_table_full=input_table
    input_table=input_table[input_table['Flag'] < 1] 
    # input_table=input_table[input_table['Flag']!='red'] 
    
    if isinstance( inputparameters,str):
        inputparameters=findtable(inputparameters,format='ecsv',path=savepath,delimiter='\s')
        inputparameters = inputparameters[inputparameters["Ref_file"] == ref_file]
    if point=='0':
        # CALCULATE THE MICRO AND ABUNDANCES FOR A CERTAIN T AND G

        fixed_T=inputparameters['T'][-1]
        fixed_g_exact=inputparameters['g'][-1]
        fixed_He=inputparameters['He'][-1]
        fixed_Q=inputparameters['Q'][-1]
        fixed_b=inputparameters['b'][-1]
        fixed_Terr=inputparameters['Terr'][-1]
        FW_EW_measuretab=findtable(FW_EW_measuretab)
        garrFW=sorted(set(FW_EW_measuretab['g']))
        fixed_g=garrFW[min(range(len(garrFW)),key=lambda i : abs(int(garrFW[i])-int(fixed_g_exact)))] #NEAREST
        grid = FW_EW_measuretab["Grid_name"][FW_EW_measuretab["He"] == str(fixed_He)][0]
        
    elif point=='1':
        # CALCULATE THE ABUNDANCES WITH FIXED MICRO TO ACCOUNT FOR THE UPPER UNCERTAINTY IN T, G

        fixed_T=inputparameters['T+'][-1]
        fixed_g_exact=inputparameters['g+'][-1] 
        fixed_He=inputparameters['He'][-1]
        fixed_Q=inputparameters['Q'][-1]
        fixed_b=inputparameters['b'][-1]
        fixed_Terr=inputparameters['Terr'][-1]
        FW_EW_measuretab=findtable(FW_EW_measuretab)
        garrFW=sorted(set(FW_EW_measuretab['g']))
        fixed_g=garrFW[min(range(len(garrFW)),key=lambda i : abs(int(garrFW[i])-int(fixed_g_exact)))] #NEAREST

    elif point=='2':
        # CALCULATE THE ABUNDANCES WITH FIXED MICRO TO ACCOUNT FOR THE LOWER UNCERTAINTY IN T, G

        fixed_T=inputparameters['T-'][-1]
        fixed_g_exact=inputparameters['g-'][-1]
        fixed_He=inputparameters['He'][-1]
        fixed_Q=inputparameters['Q'][-1]
        fixed_b=inputparameters['b'][-1]
        fixed_Terr=inputparameters['Terr'][-1]
        FW_EW_measuretab=findtable(FW_EW_measuretab)
        garrFW=sorted(set(FW_EW_measuretab['g']))
        fixed_g=garrFW[min(range(len(garrFW)),key=lambda i : abs(int(garrFW[i])-int(fixed_g_exact)))] #NEAREST   
    elif point=='Q':
        fixed_T=inputparameters['T'][-1]
        fixed_g_exact=inputparameters['g'][-1]
        fixed_He=inputparameters['He'][-1]
        fixed_Q=inputparameters['Q'][-1]
        fixed_b=inputparameters['b'][-1]
        fixed_Terr=inputparameters['Terr'][-1]
        FW_EW_measuretab=findtable(FW_EW_measuretab)
        garrFW=sorted(set(FW_EW_measuretab['g']))
        fixed_g=garrFW[min(range(len(garrFW)),key=lambda i : abs(int(garrFW[i])-int(fixed_g_exact)))] #NEAREST
        
    if microreduced=='fixed':
        micro_fixed=inputparameters['micro'][-1]
    else:
        micro_fixed=None

    #nmicro=int((len(input_table.colnames)-4)/2)
    nmicro=5
    nlines=len(input_table)
    slope=np.zeros(nmicro)
    Aslope=np.zeros(nmicro)
    absc=np.zeros(nmicro)    
    micro=np.zeros(nmicro,dtype=int)
    aO=np.zeros(nmicro)
    AaO=np.zeros(nmicro)
    n=m

    fulllinestab = findtable(linetable,format='csv',path=savepath)
    # mult_highlight = ("%d"%(int(mult_highlight))) if len(str(mult_highlight))>1 else ("0%d"%(int(mult_highlight)))
 
    #OJO: MULT_HIGHLIGHT NO ESTÁ PREPARADO PARA TENER MÁS DE 99 MULTIPLETES! TIENEN QUE SER DOS CIFRAS, SI NO CAMBIAR TODO

    if el_input=='Si':
        mask_col1=((np.char.startswith(np.array(input_table['Line'],dtype=str),'SiII')) & 
                ~(np.char.startswith(np.array(input_table['Line'],dtype=str),'SiIII')))          
        mask_col2=((np.char.startswith(np.array(input_table['Line'],dtype=str),'SiIII')))
        mask_col3=((np.char.startswith(np.array(input_table['Line'],dtype=str),'SiIV')))
        mask_col4=((input_table['Line']=='SiII6347'))
        mask_col5=((input_table['Line']=='SiIII4552') | (input_table['Line']=='SiIII4567') |
                    (input_table['Line']=='SiIII4574') )
        mask_col6=((input_table['Line']=='SiIV4116'))

        linecolor_col1=colors["sky blue"]
        linecolor_col2=colors["light orange"]
        linecolor_col3=colors["light green"]
        linecolor_col4=colors["blue"]
        linecolor_col5=colors["orange"]
        linecolor_col6=colors["bluish green"]         

        label_col1='SiII'
        label_col2='SiIII'
        label_col3='SiIV'
        label_col4='SiII6347'
        label_col5='SiIII4552-67-74'
        label_col6='SiIV4116'

        mask_col1_full=((np.char.startswith(np.array(input_table_full['Line'],dtype=str),'SiII')) & 
                ~(np.char.startswith(np.array(input_table_full['Line'],dtype=str),'SiIII')))          
        mask_col2_full=((np.char.startswith(np.array(input_table_full['Line'],dtype=str),'SiIII')))
        mask_col3_full=((np.char.startswith(np.array(input_table_full['Line'],dtype=str),'SiIV')))

        if mult_highlight is not None:
            mask_highlight = np.array([True if mult_highlight in fulllinestab["mult_FW"][fulllinestab["label"] == (input_table['Line'][i])][0].split("+")    
                                    else False for i in range(len(input_table['Line']))], 
                                    dtype=bool)
            linecolor_highlihgt = colors["reddish purple"]
            label_highlight = '-'.join(list(input_table['Line'][mask_highlight]) )

            mask_all_cols=[mask_col1,mask_col2,mask_col3,mask_col4,mask_col5,mask_col6, mask_highlight]    
            linecolor_all_cols=[linecolor_col1,linecolor_col2,linecolor_col3,linecolor_col4,linecolor_col5,linecolor_col6, linecolor_highlihgt]
            label_all_cols=[label_col1,label_col2,label_col3,label_col4,label_col5,label_col6, label_highlight]

            mask_highlight_full = np.array([True if mult_highlight in fulllinestab["mult_FW"][fulllinestab["label"] == (input_table_full['Line'][i])][0].split("+")    
                                    else False for i in range(len(input_table_full['Line']))], 
                                    dtype=bool) 
            mask_all_cols_full=[mask_col1_full,mask_col2_full,mask_col3_full,mask_highlight_full]
            linecolor_all_cols_full = ["grey", "grey", "grey", colors["reddish purple"]]

        else:                
            mask_all_cols=[mask_col1,mask_col2,mask_col3,mask_col4,mask_col5,mask_col6]    
            linecolor_all_cols=[linecolor_col1,linecolor_col2,linecolor_col3,linecolor_col4,linecolor_col5,linecolor_col6]
            label_all_cols=[label_col1,label_col2,label_col3,label_col4,label_col5,label_col6]
            
            mask_all_cols_full=[mask_col1_full,mask_col2_full,mask_col3_full]
            linecolor_all_cols_full = ["grey", "grey", "grey"]
        
        label_counter_col=np.zeros(len(mask_all_cols))

    elif el_input=='C':
        mask_col1=((np.char.startswith(np.array(input_table['Line'],dtype=str),'CII')) & 
                ~(np.char.startswith(np.array(input_table['Line'],dtype=str),'CIII')))
        mask_col2=((np.char.startswith(np.array(input_table['Line'],dtype=str),'CIII')))
        mask_col3=((np.char.startswith(np.array(input_table['Line'],dtype=str),'CIV')))

        linecolor_col1=colors["sky blue"]
        linecolor_col2=colors["orange"]
        linecolor_col3=colors["bluish green"]
        
        label_col1='CII'
        label_col2='CIII'
        label_col3='CIV' 
                    
        mask_col1_full=((np.char.startswith(np.array(input_table_full['Line'],dtype=str),'CII')) & 
                ~(np.char.startswith(np.array(input_table_full['Line'],dtype=str),'CIII')))
        mask_col2_full=((np.char.startswith(np.array(input_table_full['Line'],dtype=str),'CIII')))
        mask_col3_full=((np.char.startswith(np.array(input_table_full['Line'],dtype=str),'CIV')))
        
        if mult_highlight is not None:
            mask_highlight = np.array([True if mult_highlight in fulllinestab["mult_FW"][fulllinestab["label"] == (input_table['Line'][i])][0].split("+")    
                                    else False for i in range(len(input_table['Line']))], 
                                    dtype=bool)
            linecolor_highlihgt = colors["reddish purple"]
            label_highlight = '-'.join(list(input_table['Line'][mask_highlight]) )

            mask_all_cols=[mask_col1,mask_col2,mask_col3,mask_highlight]    
            linecolor_all_cols=[linecolor_col1,linecolor_col2,linecolor_col3,linecolor_highlihgt]
            label_all_cols=[label_col1,label_col2,label_col3,label_highlight]

            mask_highlight_full = np.array([True if mult_highlight in fulllinestab["mult_FW"][fulllinestab["label"] == (input_table_full['Line'][i])][0].split("+")    
                                    else False for i in range(len(input_table_full['Line']))], 
                                    dtype=bool) 
            mask_all_cols_full=[mask_col1_full,mask_col2_full,mask_col3_full,mask_highlight_full]
            linecolor_all_cols_full = ["grey", "grey", "grey", colors["reddish purple"]]                

        else:                
            mask_all_cols=[mask_col1,mask_col2,mask_col3]    
            linecolor_all_cols=[linecolor_col1,linecolor_col2,linecolor_col3]
            label_all_cols=[label_col1,label_col2,label_col3]
            mask_all_cols_full=[mask_col1_full,mask_col2_full,mask_col3_full]
            linecolor_all_cols_full = ["grey", "grey", "grey"]                  


        label_counter_col=np.zeros(len(mask_all_cols))                     
        
    elif el_input=='N':
        mask_col1=((np.char.startswith(np.array(input_table['Line'],dtype=str),'NII')) & 
                ~(np.char.startswith(np.array(input_table['Line'],dtype=str),'NIII')))
        mask_col2=((np.char.startswith(np.array(input_table['Line'],dtype=str),'NIII')))
        mask_col3=((np.char.startswith(np.array(input_table['Line'],dtype=str),'NIV')))

        linecolor_col1=colors["sky blue"]
        linecolor_col2=colors["orange"]
        linecolor_col3=colors["bluish green"]            
                    
        label_col1='NII'
        label_col2='NIII'
        label_col3='NIV' 
                
        mask_col1_full=((np.char.startswith(np.array(input_table_full['Line'],dtype=str),'NII')) & 
                ~(np.char.startswith(np.array(input_table_full['Line'],dtype=str),'NIII')))
        mask_col2_full=((np.char.startswith(np.array(input_table_full['Line'],dtype=str),'NIII')))
        mask_col3_full=((np.char.startswith(np.array(input_table_full['Line'],dtype=str),'NIV')))
        

        if mult_highlight is not None:
            mask_highlight = np.array([True if mult_highlight in fulllinestab["mult_FW"][fulllinestab["label"] == (input_table['Line'][i])][0].split("+")    
                                    else False for i in range(len(input_table['Line']))], 
                                    dtype=bool)
            linecolor_highlihgt = colors["reddish purple"]
            label_highlight = '-'.join(list(input_table['Line'][mask_highlight]) )

            mask_all_cols=[mask_col1,mask_col2,mask_col3,mask_highlight]    
            linecolor_all_cols=[linecolor_col1,linecolor_col2,linecolor_col3,linecolor_highlihgt]
            label_all_cols=[label_col1,label_col2,label_col3,label_highlight]

            mask_highlight_full = np.array([True if mult_highlight in fulllinestab["mult_FW"][fulllinestab["label"] == (input_table_full['Line'][i])][0].split("+")    
                                    else False for i in range(len(input_table_full['Line']))], 
                                    dtype=bool) 
            mask_all_cols_full=[mask_col1_full,mask_col2_full,mask_col3_full,mask_highlight_full]
            linecolor_all_cols_full = ["grey", "grey", "grey", colors["reddish purple"]]                


        else:                
            mask_all_cols=[mask_col1,mask_col2,mask_col3]    
            linecolor_all_cols=[linecolor_col1,linecolor_col2,linecolor_col3]
            label_all_cols=[label_col1,label_col2,label_col3]     
            mask_all_cols_full=[mask_col1_full,mask_col2_full,mask_col3_full]
            linecolor_all_cols_full = ["grey", "grey", "grey"]                               
        
        label_counter_col=np.zeros(len(mask_all_cols))                     
                    
    elif el_input=='O':
        mask_col1=((np.char.startswith(np.array(input_table['Line'],dtype=str),'OII')) & 
                ~(np.char.startswith(np.array(input_table['Line'],dtype=str),'OIII')))            
        mask_col2=((np.char.startswith(np.array(input_table['Line'],dtype=str),'OIII')))
        mask_col3=((np.char.startswith(np.array(input_table['Line'],dtype=str),'OIV')))

        linecolor_col1=colors["sky blue"]
        linecolor_col2=colors["orange"]
        linecolor_col3=colors["bluish green"]   
        
        label_col1='OII'
        label_col2='OIII'
        label_col3='OIV' 
    
        mask_col1_full=((np.char.startswith(np.array(input_table_full['Line'],dtype=str),'OII')) & 
                ~(np.char.startswith(np.array(input_table_full['Line'],dtype=str),'OIII')))            
        mask_col2_full=((np.char.startswith(np.array(input_table_full['Line'],dtype=str),'OIII')))
        mask_col3_full=((np.char.startswith(np.array(input_table_full['Line'],dtype=str),'OIV')))
        
        if mult_highlight is not None:
            mask_highlight = np.array([True if mult_highlight in fulllinestab["mult_FW"][fulllinestab["label"] == (input_table['Line'][i])][0].split("+")    
                                    else False for i in range(len(input_table['Line']))], 
                                    dtype=bool)
            linecolor_highlihgt = colors["reddish purple"]
            label_highlight = '-'.join(list(input_table['Line'][mask_highlight]) )

            mask_all_cols=[mask_col1,mask_col2,mask_col3,mask_highlight]    
            linecolor_all_cols=[linecolor_col1,linecolor_col2,linecolor_col3,linecolor_highlihgt]
            label_all_cols=[label_col1,label_col2,label_col3,label_highlight]
            
            mask_highlight_full = np.array([True if mult_highlight in fulllinestab["mult_FW"][fulllinestab["label"] == (input_table_full['Line'][i])][0].split("+")    
                                    else False for i in range(len(input_table_full['Line']))], 
                                    dtype=bool) 
            mask_all_cols_full=[mask_col1_full,mask_col2_full,mask_col3_full,mask_highlight_full]
            linecolor_all_cols_full = ["grey", "grey", "grey", colors["reddish purple"]]                

        else:                
            mask_all_cols=[mask_col1,mask_col2,mask_col3]    
            linecolor_all_cols=[linecolor_col1,linecolor_col2,linecolor_col3]
            label_all_cols=[label_col1,label_col2,label_col3]

            mask_all_cols_full=[mask_col1_full,mask_col2_full,mask_col3_full]
            linecolor_all_cols_full = ["grey", "grey", "grey"]                   
        
        label_counter_col=np.zeros(len(mask_all_cols))   


    nlines_masks =np.array( [np.sum(mask_all_cols[i]) for i in range(0,3)])
    label_ionmaxlines = label_all_cols[np.argmax(nlines_masks)]
    labels_notionmaxlines = np.delete(label_all_cols[0:3], np.argmax(nlines_masks),-1)


    if ref_file in fulllinestab.colnames:
        if "all" in linesvec:
            linesvec = list(fulllinestab['label'][fulllinestab[ref_file]<1])
        else:
            linesvec = [i for i in linesvec if i in list(fulllinestab['label'][fulllinestab[ref_file]<1])]

    mask_linesvec = np.array([True if i in linesvec else False for i in input_table['Line']])

    print("Calculating the microturbulence by looking for a constant abundance in the following lines:")
    print([str(x) for x in linesvec])

    n=m

    range_eSi=np.array([100,-100],dtype=float)
    for i in range(nmicro):
        
        if np.amax(input_table[input_table.colnames[n]])>range_eSi[1]:
            range_eSi[1]=np.amax(input_table[input_table.colnames[n]]+input_table[input_table.colnames[n+1]])+0.2
        if np.amin(input_table[input_table.colnames[n]])<range_eSi[0]:
            range_eSi[0]=np.amin(input_table[input_table.colnames[n]]-input_table[input_table.colnames[n+1]])-0.2
        n=n+2    
    n=m       

    flag_full = np.array(["--------"]*len(input_table_full["Line"]))
    
    for i in range(len(input_table_full["Line"])):
        ab = [input_table_full[coli][i] for coli in input_table_full.colnames if coli.startswith("ab") ]
        ab = np.array(ab, dtype=float)
        ab_FW = np.array(sorted(set(FW_EW_measuretab["abun_%s"%el_input])), dtype=float)
        maskab = ((ab>np.amin(ab_FW)) & (ab<np.amax(ab_FW)))
        flag_full[i] = "OK" if np.sum(maskab)>0 else "EXTRAPOL"
    
    if micro_fixed==None:    
        if do_pdf==True:           
            rows=4
            cols=nmicro
            
            fig=plt.figure(figsize=(18,12),constrained_layout=True)
            gs = fig.add_gridspec(rows, cols)
        
        #cada vez que añada una figura: ax=fig.add_subplot(gs[rowstart:rowend,colstart:colend])
        
        #PLOT 0: E(O) VS LOG EW SEGÚN MICROTURBULENCIA
        
        for i in range(nmicro):
            #para cada micro, plotear EW vs A(0) según la línea                        
            ew=input_table['EW']
            Aew=input_table['AEW']
            abSim1=input_table[input_table.colnames[n]]
            AabSim1=input_table[input_table.colnames[n+1]]
            
            logew=np.log10(ew)
            logAew=Aew/ew

            ew_full=input_table_full['EW']
            Aew_full=input_table_full['AEW']
            abSim1_full=input_table_full[input_table_full.colnames[n]]
            AabSim1_full=input_table_full[input_table_full.colnames[n+1]]
            
            logew_full=np.log10(ew_full)
            logAew_full=Aew_full/ew_full
            
            if microreduced=='yes':
                #abundancia media para cada valor de micro
                aO[i]=np.average(abSim1[mask_linesvec]) #,weights=1/AabSim1[mask_linesvec]**2 #pesada por errores 
                AaO[i]=np.std(abSim1[mask_linesvec]) #solo desviación típica 
                x_fit = np.linspace(0, np.amax(ew[mask_linesvec])+1, 100)
                
                p,pcov  = np.polyfit(ew[mask_linesvec],abSim1[mask_linesvec],1,cov='unscaled',w=1./AabSim1[mask_linesvec]) # ,w=1./AabSim1[mask_linesvec]
                perr = np.sqrt(np.diag(pcov))
                slope[i]=p[0]
                absc[i]=p[1]
                Aslope[i]=perr[0]
                
            if microreduced=='montecarlo': #(el más usado)
                n_mc=len(EWvec[0])
                slope_mc=np.zeros(n_mc,dtype=float)
                absc_mc=np.zeros(n_mc,dtype=float)

                #abundancia media para cada valor de micro
                for imc in range(n_mc):
                    EWvec_mc = [EWvec[jj][imc] for jj in range(len(linesvec))]
                    abvec_mc = [abvec[i][jj][imc] for jj in range(len(linesvec))]

                    # EWvec_mc=[EWvec[0][imc],EWvec[1][imc],EWvec[2][imc]]
                    # abvec_mc=[abvec[i][0][imc],abvec[i][1][imc],abvec[i][2][imc]]

                    p=np.polyfit(EWvec_mc,abvec_mc,1)
                    slope_mc[imc]=p[0]
                    absc_mc[imc]=p[1]
                
                slope[i]=np.mean(slope_mc)
                absc[i]=np.mean(absc_mc)               
                Aslope[i]=np.std(slope_mc)  
                
                aO[i]=np.average(abSim1[mask_linesvec]) 
                AaO[i]=np.std(abSim1[mask_linesvec])
                
            else: #usamos todas las líneas (no solemos...)
                
                #abundancia media para cada valor de micro
                aO[i]=np.average(abSim1) #,weights=1/AabSim1**2 #pesada por errores 
                AaO[i]=np.std(abSim1) #solo desviación típica 
                x_fit = np.linspace(0, np.amax(logew)+1, 100)
                
                p,pcov  = np.polyfit(ew,abSim1,1,cov='unscaled') #,w=1./AabSim1 #cov=1/sigma: reliable estimate of the uncertainty 
                perr = np.sqrt(np.diag(pcov))
                slope[i]=p[0]
                Aslope[i]=perr[0]
                
            #micro[i]=int(input_table.colnames[n].split('abSim')[1])
            micro[i]=int(input_table.colnames[n].split('ab'+el_input+'m')[1])
            
            if do_pdf==True:
                ax=fig.add_subplot(gs[0,i])

                for i_cols in range(len(mask_all_cols_full)):
                    ax.errorbar(ew_full[mask_all_cols_full[i_cols]],abSim1_full[mask_all_cols_full[i_cols]], 
                                xerr=Aew_full[mask_all_cols_full[i_cols]], fmt='o', 
                                color=linecolor_all_cols_full[i_cols],ecolor='black', ms=10,
                                capsize=2,elinewidth=1, alpha=0.5,zorder=1) 
                                        
                    ax.scatter(ew_full[((mask_all_cols_full[i_cols]) & (flag_full == "EXTRAPOL"))],
                                abSim1_full[((mask_all_cols_full[i_cols]) & (flag_full == "EXTRAPOL"))], 
                                marker="o", s=100,color="none",edgecolors="red",linewidths=1, 
                                zorder=500)


                for i_cols in range(len(mask_all_cols)):
                    ax.errorbar(ew[mask_all_cols[i_cols]],abSim1[mask_all_cols[i_cols]], 
                                xerr=Aew[mask_all_cols[i_cols]], fmt='o',
                                color=linecolor_all_cols[i_cols],ecolor='black', ms=10,
                                capsize=2,elinewidth=1)  

                    #, yerr=AabSim1[mask_all_cols[i_cols]] 
                    
                    ###### CHECK DIAGONAL UNCERTAINTIES IN EPS(SI) VS EW ##########
                    
                    for i_point in range(len(AabSim1[mask_all_cols[i_cols]])):
                        ax.plot([ (ew[mask_all_cols[i_cols]])[i_point]-np.abs(Aew[mask_all_cols[i_cols]])[i_point],
                                    (ew[mask_all_cols[i_cols]])[i_point]+np.abs(Aew[mask_all_cols[i_cols]])[i_point]  ],
                                   [ (abSim1[mask_all_cols[i_cols]])[i_point]-np.abs(AabSim1[mask_all_cols[i_cols]])[i_point],
                                    (abSim1[mask_all_cols[i_cols]])[i_point]+np.abs(AabSim1[mask_all_cols[i_cols]])[i_point] ],
                                   color='k',marker='_' ,lw=0.75,ms=6 )
                
                xx=np.linspace(0,np.amax(ew)+50)
                #ax.plot(xx,np.polyval(p,xx),'r--')
                ax.plot(xx,slope[i]*xx+absc[i],'--', c = colors["vermillion"]) 
                ax.plot(xx,aO[i]*np.ones(50),'k--',alpha=0.5)
                ax.plot(xx,(aO[i]+AaO[i])*np.ones(50),'k--',alpha=0.25)
                ax.plot(xx,(aO[i]-AaO[i])*np.ones(50),'k--',alpha=0.25)
    
                ax.set_xlabel(r'$EW$ [m$\AA$]')
                if i<1:
                    ax.set_ylabel(r'$\epsilon$ (%s)'%(el_input))
                else:
                    ax.set_yticks([])
                ax.set_ylim(range_eSi)
                ewrange = [0,np.amax(ew)+10 ]
                ax.set_xlim(ewrange)
                #ax.legend(fancybox=True, framealpha=0.7)
                ax.set_title(r'$\xi$={:d} km/s'.format(micro[i]))
            n=n+2

         #PLOT 1: SLOPE VS MICRO -> CÁLCULO DE MICRO PARA SLOPE = 0 
        order=2

        fit=np.polyfit(micro,slope,order)
        fit_sup=np.polyfit(micro,slope+Aslope,order)
        fit_inf=np.polyfit(micro,slope-Aslope,order)
        
        micro_val=np.roots(fit)[(np.abs(np.roots(fit)-np.median(micro))).argmin()]
        micro_sup=np.roots(fit_sup)[(np.abs(np.roots(fit_sup)-np.median(micro))).argmin()]
        micro_inf=np.roots(fit_inf)[(np.abs(np.roots(fit_inf)-np.median(micro))).argmin()]

        if micro_val<1:
            micro_val=1
        if micro_val-micro_inf<0:
            micro_inf=0
            
        Amicro=[[np.abs(micro_inf-micro_val)],[np.abs(micro_sup-micro_val)]]
        Amicro_sup = np.abs(micro_sup-micro_val)
        Amicro_inf = np.abs(micro_inf-micro_val)
        Amicro_val=np.mean(Amicro)
               
        xx=np.linspace(np.amin(micro)-1,np.amax(micro)+1)
               
        if do_pdf==True:
            ax=fig.add_subplot(gs[1,0:int(np.floor(cols/2))])
            ax.errorbar(micro,slope,yerr=Aslope,fmt='o',color='k',ecolor='k',capsize=2,elinewidth=1)
    
            ax.plot(xx,np.polyval(fit,xx),'--', c = colors["vermillion"])
            ax.plot(xx,np.polyval(fit_sup,xx),'--',c = colors["vermillion"], alpha=0.3)
            ax.plot(xx,np.polyval(fit_inf,xx),'--',c = colors["vermillion"], alpha=0.3)        
            #ax.fill_between(xx, np.polyval(fit_sup,xx), np.polyval(fit_inf,xx), alpha=.25,color='yellow')            
            
            ax.plot(np.linspace(np.amin(micro)-1,np.amax(micro)+1),0*np.ones(50),'k--',alpha=0.3)
            ax.errorbar(micro_val,0,xerr=Amicro,fmt='none',ecolor=colors["blue"],capsize=2,zorder=99)
            ax.scatter(micro_val,0,color = colors["blue"],marker="s",s=100,zorder=100,
                       label=r'$\xi$=%.2f$\pm$ %.2f km/s'%(micro_val,Amicro_val))
            ax.set_xlim(np.amin(micro)-1,np.amax(micro)+1)
            ax.set_ylim(np.amin(slope-Aslope),np.amax(slope+Aslope))
            ax.set_xlabel(r'$\xi$ [km/s]')
            ax.set_ylabel('Slope')
            ax.legend(fancybox=True, framealpha=0.7)
        
        #PLOT 2: AB VS MICRO -> CÁLCULO DE AB PARA MICRO CALCULADA
        
        order=2
        
        fit=np.polyfit(micro, aO, order) 
        fit_sup=np.polyfit(micro,aO+AaO,order) 
        fit_inf=np.polyfit(micro,aO-AaO,order) 
        
        aO_val=np.polyval(fit,micro_val)
        # aO_sup=np.polyval(fit,micro_val-Amicro_val)
        # aO_inf=np.polyval(fit,micro_val+Amicro_val)
        aO_sup=np.polyval(fit,micro_val-Amicro_inf)
        aO_inf=np.polyval(fit,micro_val+Amicro_sup)
        AaO_v=[[np.abs(aO_inf-aO_val)],[np.abs(aO_sup-aO_val)]]
        AaO_val=np.mean(AaO_v)
        AaO_m_inf = np.abs(aO_inf-aO_val)
        AaO_m_sup = np.abs(aO_sup-aO_val)

        xx=np.linspace(np.amin(micro)-1,np.amax(micro)+1)
        
        if do_pdf==True:
            ax=fig.add_subplot(gs[2,0:int(np.floor(cols/2))])
            ax.errorbar(micro,aO,yerr=AaO,fmt='o',color='k',ecolor='black',capsize=2,elinewidth=1)
    
            ax.plot(xx,np.polyval(fit,xx),'--', c = colors["vermillion"])
            ax.plot(xx,np.polyval(fit_sup,xx),'--',c = colors["vermillion"], alpha=0.3)
            ax.plot(xx,np.polyval(fit_inf,xx),'--',c = colors["vermillion"], alpha=0.3)
            ax.axvline(micro_val, color='black',ls='--', alpha=0.5)
            ax.axvline(micro_val+Amicro_val, color='black',ls='--', alpha=0.3)
            ax.axvline(micro_val-Amicro_val, color='black',ls='--', alpha=0.3)
            ax.plot(np.linspace(np.amin(micro)-1,np.amax(micro)+1),aO_val*np.ones(50),'k--',alpha=0.5)
            ax.plot(np.linspace(np.amin(micro)-1,np.amax(micro)+1),aO_sup*np.ones(50),'k--',alpha=0.3)
            ax.plot(np.linspace(np.amin(micro)-1,np.amax(micro)+1),aO_inf*np.ones(50),'k--',alpha=0.3)
                        
            ax.errorbar(micro_val,aO_val,yerr=AaO_v,fmt='none',ecolor=colors["blue"],capsize=2,zorder=99)
            ax.scatter(micro_val,aO_val,color=colors["blue"],marker="s",s=100,zorder=100,
                       label=r'log(%s/H) = %.2f $\pm$ %.2f dex'%(el_input,aO_val,AaO_val))
            ax.set_xlim(np.amin(micro)-1,np.amax(micro)+1)
            ax.set_xlabel(r'$\xi$ [km/s]')
            ax.set_ylabel(r'$\epsilon$ (%s)'%(el_input))
            ax.legend(fancybox=True, framealpha=0.7,loc='upper right')
        
            #PLOT 3: E(0) VS MICRO PARA e [km/s]'microreduced

            #Para la micro calculada, para cada línea calcular la abundancia que sale de aquí. 
          
            # ax=fig.add_subplot(gs[1:3,int(np.floor(cols/2)):])   
            ax=fig.add_subplot(gs[3,0:int(np.floor(cols/2))])
        
        abline=np.empty(len(input_table))
        Aabline=np.empty(len(input_table))

        for j,line in zip(range(len(input_table['Line'])),input_table['Line']):
            ab=np.empty(len(micro))
            Aab=np.empty(len(micro))
            n=m
            ########################### COLORES ##############################
            for i_cols in range(len(mask_all_cols)):
                if mask_all_cols[-1-i_cols][j]:
                    linecolori=linecolor_all_cols[-1-i_cols]
                    if label_counter_col[-1-i_cols]<1:
                        label=label_all_cols[-1-i_cols]
                        label_counter_col[-1-i_cols]=label_counter_col[-1-i_cols]+1
                    else:
                        label='no'
                    break
                else:
                    pass
           
            for l in range(len(micro)):
                ab[l]=input_table[input_table.colnames[n]][input_table['Line']==line][0]
                Aab[l]=input_table[input_table.colnames[n+1]][input_table['Line']==line][0]
                n=n+2 
            
            #COLORES AQUÍ (+ PONER EL LABEL UNA VEZ?)
            if do_pdf==True:
                if label=='no':
                    ax.plot(micro,ab,'.-',color=linecolori)
                else:
                    ax.plot(micro,ab,'.-',color=linecolori,label=label)
                    
                ax.errorbar(micro,ab,yerr=Aab,fmt='o',color=linecolori,ecolor='black',capsize=2,elinewidth=1)
                ax.scatter(micro,ab,color='k')
            fit=np.polyfit(micro, ab, order)
            fit_sup=np.polyfit(micro,ab+Aab,order)
            fit_inf=np.polyfit(micro,ab-Aab,order)
            abline[j]=np.polyval(fit,micro_val)
            abline_sup=np.polyval(fit_sup,micro_val)
            abline_inf=np.polyval(fit_inf,micro_val)
            Aabline[j]=np.mean([np.abs(abline[j]-abline_sup),np.abs(abline[j]-abline_inf)])

        #! + añadir líneas que NO están en input_table("Line")
         
        if do_pdf==True:
            ax.axvline(micro_val, color='black',ls='--', alpha=0.5,
                       label=r'$\xi$=%.2f$\pm$ %.2f km/s'%(micro_val,Amicro_val))
            ax.axvline(micro_val+Amicro_val, color='black',ls='--', alpha=0.3)
            ax.axvline(micro_val-Amicro_val, color='black',ls='--', alpha=0.3)

            
            ax.set_xlim(np.amin(micro)-1, np.amax(micro)+1)
            ax.set_xlabel(r'$\xi$ [km/s]')
            ax.set_ylabel(r'$\epsilon$ (%s)'%(el_input)) 
    
        #######################################################################
        
        #Calculamos la abundancia de cada línea usando la micro calculada:
        collines=['ab'+el_input+'_'+linei for linei in input_table_full['Line']]
        ablinefull=np.empty(len(input_table_full['Line']))
        ewfull=input_table_full['EW']
        logewfull=np.log10(input_table_full['EW'])
        for j,line in zip(range(len(input_table_full['Line'])),input_table_full['Line']):
            ab=np.empty(len(micro))
            Aab=np.empty(len(micro))
            n=m   
           
            for l in range(len(micro)):                
                ab[l]=input_table_full[input_table_full.colnames[n]][input_table_full['Line']==line][0]
                Aab[l]=input_table_full[input_table_full.colnames[n+1]][input_table_full['Line']==line][0]
                n=n+2 
            fit=np.polyfit(micro, ab, order)
            ablinefull[j]=np.polyval(fit,micro_val)   

        #######################################################################

        #PLOT 4: 
        ab_final=np.median(abline) #weights=1/Aabline**2
        Aab_final=np.std(abline) 

        if do_pdf==True:
            ax.plot(np.linspace(np.amin(micro)-1,np.amax(micro)+1),ab_final*np.ones(50),'k--',alpha=0.5,
                    label=r'log(%s/H) = %.2f $\pm$ %.2f dex'%(el_input,ab_final,Aab_final))
            ax.plot(np.linspace(np.amin(micro)-1,np.amax(micro)+1),(ab_final+np.std(abline))*np.ones(50),'k--',alpha=0.3)
            ax.plot(np.linspace(np.amin(micro)-1,np.amax(micro)+1),(ab_final-np.std(abline))*np.ones(50),'k--',alpha=0.3)            
            ax.legend(fancybox=True, framealpha=0.7, loc='upper right')

            x_fit = np.linspace(0, np.amax(ewfull)+10, 100)
            ax=fig.add_subplot(gs[1:3,int(np.floor(cols/2)):])
            
            ax.plot(x_fit,ab_final*np.ones(100),'k--',alpha=0.5,
                    label=r'log(%s/H) = %.2f $\pm$ %.2f dex'%(el_input,ab_final,Aab_final))
            ax.plot(x_fit,ab_final+np.std(abline)*np.ones(100),'k--',alpha=0.3)
            ax.plot(x_fit,ab_final-np.std(abline)*np.ones(100),'k--',alpha=0.3)

            #mask the outliers!

            mask_full = ( (ablinefull < ab_final + 5*Aab_final) & 
                                   (ablinefull > ab_final - 5*Aab_final)  &
                                   (ewfull < ewrange[1]) )
                                #    & (ewfull < np.median(ewfull) + 5*np.std(ewfull)) & 
                                #    (ewfull > np.median(ewfull) - 5*np.std(ewfull)) )

            ewrange = [0,np.amax(ewfull[mask_full])+10]
            abrange = [np.nanmin(ablinefull[mask_full])-0.25, np.nanmax([np.nanmax(abline+Aabline),np.nanmax(ablinefull[mask_full])])+0.25]

            ew_outliers = np.array([np.amax(ewfull[mask_full]) + 6
                            if ewfull[ii] > np.amax(ewfull[mask_full]) + 10 else ewfull[ii]
                            for ii in range(len(ewfull))])  
            
            ab_outliers = np.array([np.amax(ablinefull[mask_full]) + 0.2
                                    if ablinefull[ii]>abrange[1]
                                    else np.amin(ablinefull[mask_full]) - 0.2
                                    if ablinefull[ii]<abrange[0]
                                    else ablinefull[ii] 
                                    for ii in range(len(ablinefull))])
            
           
            marker_outliers = np.array(['^' if ablinefull[ii]>ab_final
                                        else 'v' for ii in range(len(ablinefull))]) 

            for i_cols in range(len(mask_all_cols_full)):

                ax.scatter(ew_outliers[mask_full & mask_all_cols_full[i_cols]],
                           ab_outliers[mask_full &mask_all_cols_full[i_cols]],
                           color='none',
                           edgecolors=linecolor_all_cols[i_cols], #linecolor_all_cols_full[i_cols]
                           marker='o',alpha=0.9,linewidths=3, zorder = 100, s = 100)
                
                # COMMENTING ON 24.03.2026            
                ax.scatter(ew_outliers[(mask_full & (mask_all_cols_full[i_cols]) & (flag_full == "EXTRAPOL"))],
                            ab_outliers[(mask_full & (mask_all_cols_full[i_cols]) & (flag_full == "EXTRAPOL"))], 
                            marker="o", s=150,color="none",edgecolors="red",linewidths=1, 
                            zorder=500)                
                
                #! DIBUJAR LOS OUTLIERS!
                if len(ewfull[~mask_full & mask_all_cols_full[i_cols]])>0:

                    for ii in range(len(ewfull)):
                        if ((mask_full[ii] == False) & (mask_all_cols_full[i_cols] [ii]==True)):
                            ax.scatter(ew_outliers[ii],ab_outliers[ii],color='none', edgecolors=linecolor_all_cols[i_cols], #linecolor_all_cols_full[i_cols]
                                    marker=marker_outliers[ii], alpha=0.9,linewidths=3, zorder = 100, s = 100) 

            for i_cols in range(len(mask_all_cols)):

                ax.errorbar(ew[mask_all_cols[i_cols]],abline[mask_all_cols[i_cols]],
                            fmt='o',color= linecolor_all_cols[i_cols],ecolor='k',capsize=5,elinewidth=1,
                            zorder = 500, ms = np.sqrt(100) )  
                
                #,yerr=Aabline[mask_all_cols[i_cols]]
                
                ######## CHECK DIAGONAL UNCERTAINTIES IN EPS(SI) VS EW ############
                
                for i_point in range(len(AabSim1[mask_all_cols[i_cols]])):

                        ax.plot([ (ew[mask_all_cols[i_cols]])[i_point]-np.abs(Aew[mask_all_cols[i_cols]])[i_point],
                                    (ew[mask_all_cols[i_cols]])[i_point]+np.abs(Aew[mask_all_cols[i_cols]])[i_point]  ],
                                   [ (abline[mask_all_cols[i_cols]])[i_point]-np.abs(AabSim1[mask_all_cols[i_cols]])[i_point],
                                    (abline[mask_all_cols[i_cols]])[i_point]+np.abs(AabSim1[mask_all_cols[i_cols]])[i_point] ],
                                   color='k',marker='_' ,lw=0.75,ms=15) #,s=3,linewidths=5
            ew_ordered,line_ordered, abline_ordered, flag_ordered =zip(*sorted(zip(ew_outliers,input_table_full['Line'], 
                                                                                  ab_outliers, flag_full)))
            
            for i in range(len(input_table_full)):
                # COMMENTING ON 24.03.2026  
                # coli = "red" if flag_ordered[i]=="EXTRAPOL" else "black"
                coli = "black"

                if abline_ordered[i]<=ab_final:
                    if line_ordered[i] in input_table['Line']:
                        ax.text(ew_ordered[i],abline_ordered[i],"  [%s]"%i, size=14,
                                color=coli, va="top", ha="left", alpha=1)
                    else:
                        ax.text(ew_ordered[i],abline_ordered[i],"  [%s]"%i, size=14,
                                color=coli, va="top", ha="left", alpha=0.5)                     

                else:
                    if line_ordered[i] in input_table['Line']:
                        ax.text(ew_ordered[i], abline_ordered[i],"  [%s]"%i, 
                                size=14,color=coli, va="bottom", ha="right", alpha=1)
                    else:
                        ax.text(ew_ordered[i], abline_ordered[i],"  [%s]"%i, 
                                size=14,color=coli, va="bottom", ha="right", alpha=0.5)
                        
            ax.set_ylim(abrange)
            ax.set_xlim(ewrange)      

            ax.set_ylabel(r'$\epsilon$ (%s)'%(el_input))
            ax.set_xlabel(r'$EW$ [m$\AA$]')
            ax.legend(fancybox=True, framealpha=0.4,loc='upper right')
            
            # ax=fig.add_subplot(gs[3,-1])
            ax=fig.add_subplot(gs[3,int(np.floor(cols/2)):])
            ax.set_ylim(-len(input_table),1)

            i1 = ceil(len(line_ordered)/4)
            
            for i in range(len(input_table_full)):   
                # COMMENTING ON 24.03.2026  
                # coli = "red" if flag_ordered[i]=="EXTRAPOL" else "black" 
                coli = "black"
                if i<i1:
                    if line_ordered[i] in input_table['Line']:
                        ax.text(0,-i,"[%i] %s"%(i,line_ordered[i]),fontsize='14',color=coli, alpha=1)
                    else:
                        ax.text(0,-i,"[%i] %s"%(i,line_ordered[i]),fontsize='14',color=coli, alpha=0.5)
                elif i<(i1*2):
                    if line_ordered[i] in input_table['Line']:
                        ax.text(5,-i+i1,"[%i] %s"%(i,line_ordered[i]),fontsize='14',color=coli, alpha=1)
                    else:
                        ax.text(5,-i+i1,"[%i] %s"%(i,line_ordered[i]),fontsize='14', color=coli, alpha=0.5)
                elif i<(i1*3):
                    if line_ordered[i] in input_table['Line']:
                        ax.text(10,-i+(i1*2),"[%i] %s"%(i,line_ordered[i]),fontsize='14',color=coli, alpha=1)
                    else:
                        ax.text(10,-i+(i1*2),"[%i] %s"%(i,line_ordered[i]),fontsize='14',color=coli, alpha=0.5)
                else:
                    if line_ordered[i] in input_table['Line']:
                        ax.text(15,-i+(i1*3),"[%i] %s"%(i,line_ordered[i]),fontsize='14',color=coli, alpha=1)
                    else:
                        ax.text(15,-i+(i1*3),"[%i] %s"%(i,line_ordered[i]),fontsize='14',color=coli, alpha=0.5)                        

            ax.set_xlim(-1,18)
            ax.set_ylim(-5.5,1)

            ax.axis('off')
            # ax.set_xticks([])
            # ax.set_yticks([])
            
            #print(star_name,spc_star,(fixed_T),(fixed_g),float(fixed_Q),float(fixed_b/10))
            
            fig.suptitle('%s (%s) - T=%i K, g=%.1f dex, logQ=%.1f, b=%.1f, He = %.2f \n %s'
                         %(star_name,spc_star,int(fixed_T),float(fixed_g)/100,-float(fixed_Q)/100,float(fixed_b)/10, float(fixed_He)/100, ref_file),
                         fontsize='20') 
            
            if mult_highlight is not None:
                filename=star_name+'_'+el_input+'_abundances_%s'%mult_highlight
            else:
                filename=star_name+'_'+el_input+'_abundances'
            #fig.suptitle('%s (%s) - T=%i, g=%.1f'%(star_name,spc_star,table_Tg['T'][int(Tsel)],table_Tg['g'][int(Tsel)]/100),fontsize='25')
            #filename=star_name+'_'+el_input+'_abundances'+Tsel
            
            if save=='y':
                if point=='0':
                    # plt.savefig(savepath+filename+'_MAIN_%s.png'% time.strftime('%Y%m%d_%H%M%S'))
                    fig.savefig(savepath+filename+'_MAIN.png')
                elif point=='1':
                    # plt.savefig(savepath+filename+'_1_%s.png'% time.strftime('%Y%m%d_%H%M%S'))
                    fig.savefig(savepath+filename+'_1.png')
                  
                elif point=='2':
                    # plt.savefig(savepath+filename+'_2_%s.png'% time.strftime('%Y%m%d_%H%M%S'))                
                    fig.savefig(savepath+filename+'_2.png')
                
            #plt.savefig(savepath+filename+'.png')
            plt.show(block=False)

        
        ######################################################################

        # if point==0 do the table and overwrite. If point==1 or 2, add a new row
        # (if there are already enough rows overwrite the last one)

        #if point==0 and el_input !=Si, we might need to read the old table + add columns... -> because Si is how we get T, g. 
        #if point =1,2 and el_input != Si, then we just fill the values. 
        
        if point=='0':
            final_table=Table([[star_name],[ref_file],[grid],["WIP"],[micro_val],[Amicro_val],[Amicro_inf],[Amicro_sup],[ab_final],[Aab_final],[AaO_val],[AaO_m_inf],[AaO_m_sup],[0.0],[0.0],
                               [fixed_T],[int(fixed_Terr)],[-999999.9],[-999.9],[fixed_g_exact],[-999999.9],[-999.9],[fixed_He],[fixed_Q],[fixed_b]]+
                              [[j] for j in ablinefull],
                names=(['ID','Ref_file','Grid','FLAG','micro','Amicro','Amicro_sup','Amicro_inf','ab'+el_input,'Aab'+el_input+'_s','Aab'+el_input+'_m','Aab'+el_input+'_m_inf','Aab'+el_input+'_sup','ab'+el_input+'+','ab'+el_input+'-',
                        'T','Terr','T+','T-','g','g+','g-','He','Q','b']+collines),
                dtype=(['S10',"S200","S50","S10",'float64','float64','float64','float64','float64','float64','float64','float64','float64','float64','float64','float64',
                        'int64','float64','float64','int64','float64','float64','int64','int64','int64'] + ['float64']*len(ablinefull)
                       ))
            
            # GUARDAR ESTO EN LA TABLA DE ANALISIS DE TODAS LAS ESTRELLAS CON UN FLAG DE UNFINISHED

            if save=='y':  
                    format_i = "csv" if ".csv" in parameterstable else "ecsv" if ".ecsv" in parameterstable else "ascii"
                    delimiter_i = "," if ".csv" in parameterstable else "\s" if ".ecsv" in parameterstable else "\s"    

                    parameterstable_tab = findtable(parameterstable,format=format_i,path=savepath,delimiter=delimiter_i)
                    if parameterstable_tab is not None:
                            if el_input == "Si":            
                                parameterstable_tab = parameterstable_tab[parameterstable_tab["Ref_file"]!=ref_file]
                                parameterstable_tab = vstack([parameterstable_tab,final_table])
                                parameterstable_tab.write(savepath + parameterstable,format='ascii.%s'%format_i,overwrite=True)
                            else:
                                if microreduced == "fixed":
                                    final_table_el = final_table[["Ref_file",'ab'+el_input,'Aab'+el_input+'_s','Aab'+el_input+'_m','Aab'+el_input+'_m_inf','Aab'+el_input+'_sup','ab'+el_input+'+','ab'+el_input+'-'] + collines ]
                                    colnames_new = ['ab'+el_input,'Aab'+el_input+'_s','Aab'+el_input+'_m','Aab'+el_input+'_m_inf','Aab'+el_input+'_sup','ab'+el_input+'+','ab'+el_input+'-' ] + collines
                                else:
                                    final_table.rename_columns(['micro','Amicro','Amicro_sup','Amicro_inf'],
                                                               ['micro'+el_input,'Amicro'+el_input,'Amicro_sup'+el_input,'Amicro_inf'+el_input])

                                    final_table_el = final_table[["Ref_file",'micro'+el_input,'Amicro'+el_input,'Amicro_sup'+el_input,'Amicro_inf'+el_input,
                                                                'ab'+el_input,'Aab'+el_input+'_s','Aab'+el_input+'_m','Aab'+el_input+'_m_inf','Aab'+el_input+'_sup','ab'+el_input+'+','ab'+el_input+'-']+collines]
                                    colnames_new = ['micro'+el_input,'Amicro'+el_input,'Amicro_sup'+el_input,'Amicro_inf'+el_input,'ab'+el_input,'Aab'+el_input+'_s','Aab'+el_input+'_m','Aab'+el_input+'_m_inf','Aab'+el_input+'_sup','ab'+el_input+'+','ab'+el_input+'-']+ collines

                                parameterstable_tab_onestar = parameterstable_tab[parameterstable_tab["Ref_file"]==ref_file]
                                colnames_old = parameterstable_tab_onestar.colnames 
                                colnames_old = [coli for coli in colnames_old if coli not in colnames_new]
                                parameterstable_tab_onestar = parameterstable_tab_onestar[colnames_old]

                                final_table = join(parameterstable_tab_onestar,final_table_el,keys = "Ref_file", join_type="outer")

                                parameterstable_tab = parameterstable_tab[parameterstable_tab["Ref_file"]!=ref_file]
                                parameterstable_tab = vstack([parameterstable_tab,final_table])
                                parameterstable_tab.write(savepath + parameterstable,format='ascii.%s'%format_i,overwrite=True)
                    else:
                        final_table.write(savepath + parameterstable,format='ascii.%s'%format_i,overwrite=True)

        elif point=='1':
            format_i = "csv" if ".csv" in parameterstable else "ecsv" if ".ecsv" in parameterstable else "ascii"
            delimiter_i = "," if ".csv" in parameterstable else "\s" if ".ecsv" in parameterstable else "\s"    

            parameterstable_tab = findtable(parameterstable,format=format_i,path=savepath,delimiter=delimiter_i)
            if len(parameterstable_tab[parameterstable_tab["Ref_file"] == ref_file])>0:

                parameterstable_tab['ab'+el_input+'+'][parameterstable_tab["Ref_file"] == ref_file] = ab_final
                parameterstable_tab['T+'][parameterstable_tab["Ref_file"] == ref_file] = fixed_T
                # final_table['g-']=parameterstab['g+'][-1] #??
                parameterstable_tab['g+'][parameterstable_tab["Ref_file"] == ref_file]=fixed_g_exact #? que es el que he usado para la abun 
                # GUARDAR ESTO EN LA TABLA DE ANALISIS DE TODAS LAS ESTRELLAS CON UN FLAG DE UNFINISHED
                if save=='y':
                    parameterstable_tab.write(savepath + parameterstable,format='ascii.%s'%format_i,overwrite=True)
            else:
                print("ERROR: CAN'T FIND INFO ON MAIN PARAMETERS SO CAN'T WRITE THE UNCERTAINTY")
        elif point=='2':
            format_i = "csv" if ".csv" in parameterstable else "ecsv" if ".ecsv" in parameterstable else "ascii"
            delimiter_i = "," if ".csv" in parameterstable else "\s" if ".ecsv" in parameterstable else "\s"    

            parameterstable_tab = findtable(parameterstable,format=format_i,path=savepath,delimiter=delimiter_i)
            if len(parameterstable_tab[parameterstable_tab["Ref_file"] == ref_file])>0:

                parameterstable_tab['ab'+el_input+'-'][parameterstable_tab["Ref_file"] == ref_file] = ab_final
                parameterstable_tab['T-'][parameterstable_tab["Ref_file"] == ref_file] = fixed_T
                # final_table['g-']=parameterstab['g-'][-1] #?? 
                parameterstable_tab['g-'][parameterstable_tab["Ref_file"] == ref_file]=fixed_g_exact #? que es el que he usado para la abun 
                # GUARDAR ESTO EN LA TABLA DE ANALISIS DE TODAS LAS ESTRELLAS CON UN FLAG DE UNFINISHED
                if save=='y':
                    parameterstable_tab.write(savepath + parameterstable,format='ascii.%s'%format_i,overwrite=True)
            else:
                print("ERROR: CAN'T FIND INFO ON MAIN PARAMETERS SO CAN'T WRITE THE UNCERTAINTY")            

        # borramos la opción de point = "Q"
                
        # elif point=='Q':
        #     final_table=findtable(star_name+'_'+el_input+'_parameters.ecsv',
        #                           path=savepath,format='ecsv',delimiter='\s')
            
        #     final_table['T_Q1400']=fixed_T
        #     final_table['g_Q1400']=fixed_g_exact
        #     final_table['micro_Q1400']=micro_val
        #     final_table['ab_Q1400']=ab_final

        #      # GUARDAR ESTO EN LA TABLA DE ANALISIS DE TODAS LAS ESTRELLAS CON UN FLAG DE UNFINISHED

        #     if save=='y':    
        #         final_table.write(savepath + star_name+'_'+el_input+'_parameters.ecsv',
        #               format='ascii.ecsv',overwrite=True) 

    ######################################################################################
    ############################## MICRO IS FIXED ########################################
    ######################################################################################

    else: # micro_fixed != None: micro is fixed!
    
        abline=np.empty(len(input_table))
        Aabline=np.empty(len(input_table))
        micro_fixed=float(micro_fixed)
        # Aabline=np.empty(len(input_table))   
        order=2
        n=m

        for i in range(nmicro):
            ###################################################
            # (this is just for the plot)
            ew=input_table['EW']
            Aew=input_table['AEW']
            abSim1=input_table[input_table.colnames[n]]
            AabSim1=input_table[input_table.colnames[n+1]]
            
            logew=np.log10(ew)
            logAew=Aew/ew

            # esta parte la quitamos porque microreduced == "fixed" en este caso 

            # if microreduced=='yes':
            #     #abundancia media para cada valor de micro
            #     aO[i]=np.average(abSim1[mask_linesvec]) #,weights=1/AabSim1[mask_linesvec]**2 #pesada por errores 
            #     AaO[i]=np.std(abSim1[mask_linesvec]) #solo desviación típica 
            #     x_fit = np.linspace(0, np.amax(ew[mask_linesvec])+1, 100)
                
            #     p,pcov  = np.polyfit(ew[mask_linesvec],abSim1[mask_linesvec],1,cov='unscaled',w=1./AabSim1[mask_linesvec]) # ,w=1./AabSim1[mask_linesvec]
            #     perr = np.sqrt(np.diag(pcov))
            #     slope[i]=p[0]
            #     absc[i]=p[1]
            #     Aslope[i]=perr[0]
                
            # if microreduced=='montecarlo': #(el más usado)
            #     n_mc=len(EWvec[0])
            #     slope_mc=np.zeros(n_mc,dtype=float)
            #     absc_mc=np.zeros(n_mc,dtype=float)

            #     #abundancia media para cada valor de micro
            #     for imc in range(n_mc):
            #         EWvec_mc = [EWvec[jj][imc] for jj in range(len(linesvec))]
            #         abvec_mc = [abvec[i][jj][imc] for jj in range(len(linesvec))]

            #         # EWvec_mc=[EWvec[0][imc],EWvec[1][imc],EWvec[2][imc]]
            #         # abvec_mc=[abvec[i][0][imc],abvec[i][1][imc],abvec[i][2][imc]]

            #         p=np.polyfit(EWvec_mc,abvec_mc,1)
            #         slope_mc[imc]=p[0]
            #         absc_mc[imc]=p[1]
                
            #     slope[i]=np.mean(slope_mc)
            #     absc[i]=np.mean(absc_mc)               
            #     Aslope[i]=np.std(slope_mc)  
                
            #     aO[i]=np.average(abSim1[mask_linesvec]) 
            #     AaO[i]=np.std(abSim1[mask_linesvec])
            
            #usamos todas las líneas:

            #abundancia media para cada valor de micro
            aO[i]=np.average(abSim1) #,weights=1/AabSim1**2 #pesada por errores 
            AaO[i]=np.std(abSim1) #solo desviación típica 
            x_fit = np.linspace(0, np.amax(logew)+1, 100)
            
            p,pcov  = np.polyfit(ew[mask_linesvec],abSim1[mask_linesvec],1,cov='unscaled') #,w=1./AabSim1 #cov=1/sigma: reliable estimate of the uncertainty 
            perr = np.sqrt(np.diag(pcov))
            slope[i]=p[0]
            absc[i]=p[1]
            Aslope[i]=perr[0]
            ###################################################

            micro[i]=int(input_table.colnames[n].split('ab'+el_input+'m')[1])
            n=n+2

        for j,line in zip(range(len(input_table['Line'])),input_table['Line']):
            ab=np.empty(len(micro))
            Aab=np.empty(len(micro))
            n=m
            for l in range(len(micro)):
                ab[l]=input_table[input_table.colnames[n]][input_table['Line']==line][0]
                Aab[l]=input_table[input_table.colnames[n+1]][input_table['Line']==line][0]
                n=n+2  

            fit=np.polyfit(micro, ab, order)
            fit_sup=np.polyfit(micro,ab+Aab,order)
            fit_inf=np.polyfit(micro,ab-Aab,order)
            abline[j]=np.polyval(fit,micro_fixed)
            # print(abline)
            abline_sup=np.polyval(fit_sup,micro_fixed)
            abline_inf=np.polyval(fit_inf,micro_fixed)
            Aabline[j]=np.mean([np.abs(abline[j]-abline_sup),np.abs(abline[j]-abline_inf)])
            # input('check') # ! 
        ab_final=np.median(abline) #weights=1/Aabline**2
        Aab_final=np.std(abline) 

        collines=['ab'+el_input+'_'+linei for linei in input_table_full['Line']]
        ablinefull=np.empty(len(input_table_full['Line']))
        ewfull=input_table_full['EW']
        logewfull=np.log10(input_table_full['EW'])
        for j,line in zip(range(len(input_table_full['Line'])),input_table_full['Line']):
            ab=np.empty(len(micro))
            Aab=np.empty(len(micro))
            n=m   
           
            for l in range(len(micro)):                
                ab[l]=input_table_full[input_table_full.colnames[n]][input_table_full['Line']==line][0]
                Aab[l]=input_table_full[input_table_full.colnames[n+1]][input_table_full['Line']==line][0]
                n=n+2 
            fit=np.polyfit(micro, ab, order)
            ablinefull[j]=np.polyval(fit,micro_fixed)  

        #PLOT
        if do_pdf==True:           
            rows=4
            cols=nmicro
            
            plt.rcParams.update({
                    'xtick.labelsize' : 14,
                    'ytick.labelsize' : 14,
                    'axes.titlesize' : 14,
                    'axes.titlepad' : 6,
                    "font.size": 12,
                    "axes.titlesize": 14,
                    "axes.labelsize": 14,
                    "legend.fontsize": 12,
                    "figure.titlesize": 14,
                    })
            
            fig=plt.figure(figsize=(18,13),constrained_layout=True)
            gs = fig.add_gridspec(rows, cols)
        n=m 

        #PLOT 0: E(O) VS LOG EW SEGÚN MICROTURBULENCIA
        for i in range(nmicro):

            ew=input_table['EW']
            Aew=input_table['AEW']
            abSim1=input_table[input_table.colnames[n]]
            AabSim1=input_table[input_table.colnames[n+1]]

            logew=np.log10(ew)
            logAew=Aew/ew

            if do_pdf==True:
                ax=fig.add_subplot(gs[0,i])

                for i_cols in range(len(mask_all_cols)):
                    ax.errorbar(ew[mask_all_cols[i_cols]],abSim1[mask_all_cols[i_cols]], xerr=Aew[mask_all_cols[i_cols]], 
                            fmt='o', color=linecolor_all_cols[i_cols],ecolor='black', capsize=2,elinewidth=1)  
                    
                    #, yerr=AabSim1[mask_all_cols[i_cols]] 
                    
                    ###### CHECK DIAGONAL UNCERTAINTIES IN EPS(SI) VS EW ##########
                    
                    for i_point in range(len(AabSim1[mask_all_cols[i_cols]])):
                        ax.plot([ (ew[mask_all_cols[i_cols]])[i_point]-np.abs(Aew[mask_all_cols[i_cols]])[i_point],
                                    (ew[mask_all_cols[i_cols]])[i_point]+np.abs(Aew[mask_all_cols[i_cols]])[i_point]  ],
                                   [ (abSim1[mask_all_cols[i_cols]])[i_point]-np.abs(AabSim1[mask_all_cols[i_cols]])[i_point],
                                    (abSim1[mask_all_cols[i_cols]])[i_point]+np.abs(AabSim1[mask_all_cols[i_cols]])[i_point] ],
                                   color='k',marker='_' ,lw=0.75,ms=6 )
                
                xx=np.linspace(0,np.amax(ew)+1,50)
                #ax.plot(xx,np.polyval(p,xx),'r--')
                ax.plot(xx,slope[i]*xx+absc[i],'--', c = colors["vermillion"]) 
                ax.plot(xx,aO[i]*np.ones(50),'k--',alpha=0.5)
                ax.plot(xx,(aO[i]+AaO[i])*np.ones(50),'k--',alpha=0.25)
                ax.plot(xx,(aO[i]-AaO[i])*np.ones(50),'k--',alpha=0.25)
    
                ax.set_xlabel(r'$EW$ [m$\AA$]')
                if i<1:
                    ax.set_ylabel(r'$\epsilon$ (%s)'%(el_input))
                else:
                    ax.set_yticks([])
                ax.set_ylim(range_eSi)
                #ax.legend(fancybox=True, framealpha=0.7)
                ax.set_title(r'$\xi$={:d} km/s'.format(micro[i]))
            n=n+2

        order=2

        fit=np.polyfit(micro,slope,order)
        fit_sup=np.polyfit(micro,slope+Aslope,order)
        fit_inf=np.polyfit(micro,slope-Aslope,order)
        
        micro_val=np.roots(fit)[(np.abs(np.roots(fit)-np.median(micro))).argmin()]
        micro_sup=np.roots(fit_sup)[(np.abs(np.roots(fit_sup)-np.median(micro))).argmin()]
        micro_inf=np.roots(fit_inf)[(np.abs(np.roots(fit_inf)-np.median(micro))).argmin()]

        if micro_val<1:
            micro_val=1
        if micro_val-micro_inf<0:
            micro_inf=0
            
        Amicro=[[np.abs(micro_inf-micro_val)],[np.abs(micro_sup-micro_val)]]
        Amicro_sup = np.abs(micro_sup-micro_val)
        Amicro_inf = np.abs(micro_inf-micro_val)
        Amicro_val=np.mean(Amicro)

        xx=np.linspace(np.amin(micro)-1,np.amax(micro)+1)

        if do_pdf==True: 
            # PLOT 1: SLOPE VS MICRO -> CÁLCULO DE MICRO PARA SLOPE = 0 
            ax=fig.add_subplot(gs[1,0:int(np.floor(cols/2))]) 
            ax.errorbar(micro,slope,yerr=Aslope,fmt='o',color='k',ecolor='k',capsize=2,elinewidth=1) 

            ax.plot(xx,np.polyval(fit,xx),'--', c = colors["vermillion"])
            ax.plot(xx,np.polyval(fit_sup,xx),'--',c = colors["vermillion"], alpha=0.3)
            ax.plot(xx,np.polyval(fit_inf,xx),'--',c = colors["vermillion"], alpha=0.3)        
            #ax.fill_between(xx, np.polyval(fit_sup,xx), np.polyval(fit_inf,xx), alpha=.25,color='yellow')            
            
            ax.plot(np.linspace(np.amin(micro)-1,np.amax(micro)+1),0*np.ones(50),'k--',alpha=0.2)
            ax.axvline(micro_fixed, color='black',ls='--', alpha=0.9,label=r'$\xi=%.2f$ km/s'%(micro_fixed ))

            # ax.errorbar(micro_val,0,xerr=Amicro,fmt='none',ecolor=colors["blue"],capsize=2,zorder=99)
            # ax.scatter(micro_val,0,color = colors["blue"],marker="s",s=100,zorder=100,
            #            label=r'$\xi$=%.2f$\pm$ %.2f km/s'%(micro_val,Amicro_val))
            ax.set_xlim(np.amin(micro)-1,np.amax(micro)+1)
            ax.set_ylim(np.amin(slope-Aslope),np.amax(slope+Aslope))
            ax.set_xlabel(r'$\xi$ [km/s]')
            ax.set_ylabel('Slope')
            ax.legend(fancybox=True, framealpha=0.7,loc='upper right')

            order=2
        
            fit=np.polyfit(micro, aO, order) 
            fit_sup=np.polyfit(micro,aO+AaO,order) 
            fit_inf=np.polyfit(micro,aO-AaO,order) 
            
            aO_val=np.polyval(fit,micro_val)
            # aO_sup=np.polyval(fit,micro_val-Amicro_val)
            # aO_inf=np.polyval(fit,micro_val+Amicro_val)
            aO_sup=np.polyval(fit,micro_val-Amicro_inf)
            aO_inf=np.polyval(fit,micro_val+Amicro_sup)
            AaO_v=[[np.abs(aO_inf-aO_val)],[np.abs(aO_sup-aO_val)]]
            AaO_val=np.mean(AaO_v)
            AaO_m_inf = np.abs(aO_inf-aO_val)
            AaO_m_sup = np.abs(aO_sup-aO_val)

            #PLOT 2: AB VS MICRO -> CÁLCULO DE AB PARA MICRO CALCULADA
            ax=fig.add_subplot(gs[2,0:int(np.floor(cols/2))]) 
            ax.errorbar(micro,aO,yerr=AaO,fmt='o',color='k',ecolor='black',capsize=2,elinewidth=1)
    
            ax.plot(xx,np.polyval(fit,xx),'--', c = colors["vermillion"])
            ax.plot(xx,np.polyval(fit_sup,xx),'--',c = colors["vermillion"], alpha=0.3)
            ax.plot(xx,np.polyval(fit_inf,xx),'--',c = colors["vermillion"], alpha=0.3)
            # ax.axvline(micro_fixed, color='black',ls='--', alpha=0.5)
            # ax.axvline(micro_val+Amicro_val, color='black',ls='--', alpha=0.3)
            # ax.axvline(micro_val-Amicro_val, color='black',ls='--', alpha=0.3)
            ax.plot(np.linspace(np.amin(micro)-1,np.amax(micro)+1),ab_final*np.ones(50),'k--',alpha=0.5)
            ax.plot(np.linspace(np.amin(micro)-1,np.amax(micro)+1),(ab_final+Aab_final)*np.ones(50),'k--',alpha=0.3)
            ax.plot(np.linspace(np.amin(micro)-1,np.amax(micro)+1),(ab_final-Aab_final)*np.ones(50),'k--',alpha=0.3)
            
            ax.axvline(micro_fixed, color='black',ls='--', alpha=0.9,label=r'$\xi=%.2f$ km/s'%(micro_fixed ))

            # ax.errorbar(micro_val,aO_val,yerr=AaO_v,fmt='none',ecolor=colors["blue"],capsize=2,zorder=99)
            # ax.scatter(micro_val,aO_val,color=colors["blue"],marker="s",s=100,zorder=100,
            #            label=r'log (%s/H) = %.2f $\pm$ %.2f dex'%(el_input,aO_val,AaO_val))

            ax.set_xlim(np.amin(micro)-1,np.amax(micro)+1)
            ax.set_xlabel(r'$\xi$ [km/s]')
            ax.set_ylabel(r'$\epsilon$ (%s)'%(el_input))
            ax.legend(fancybox=True, framealpha=0.7,loc='upper right')    
            
            #PLOT 3: AB VS MICRO -> CÁLCULO DE AB PARA CADA LÍNEA
            ax=fig.add_subplot(gs[3,0:int(np.floor(cols/2))])

            for j,line in zip(range(len(input_table['Line'])),input_table['Line']):
                ab=np.empty(len(micro))
                Aab=np.empty(len(micro))                
                n=m
                ########################### COLORES ##############################

                for i_cols in range(len(mask_all_cols)):
                    if mask_all_cols[-1-i_cols][j]:
                        linecolori=linecolor_all_cols[-1-i_cols]
                        if label_counter_col[-1-i_cols]<1:
                            label=label_all_cols[-1-i_cols]
                            label_counter_col[-1-i_cols]=label_counter_col[-1-i_cols]+1
                        else:
                            label='no'
                        break
                    else:
                        pass

                for l in range(len(micro)):
                    ab[l]=input_table[input_table.colnames[n]][input_table['Line']==line][0]
                    Aab[l]=input_table[input_table.colnames[n+1]][input_table['Line']==line][0]
                    n=n+2 
                            
                #COLORES AQUÍ (+ PONER EL LABEL UNA VEZ?)
                
                if label=='no':
                    ax.plot(micro,ab,'.-',color=linecolori)
                else:
                    ax.plot(micro,ab,'.-',color=linecolori,label=label)
                    
                #ax.text(6,max(ab),line)
                ax.errorbar(micro,ab,yerr=Aab,fmt='o',color=linecolori,ecolor='black',capsize=2,elinewidth=1)
                ax.scatter(micro,ab,color='k')

            ax.axvline(micro_fixed, color='black',ls='--', alpha=0.9,label=r'$\xi$=%.2f km/s'%(micro_fixed ))

            ax.plot(np.linspace(np.amin(micro)-1,np.amax(micro)+1),ab_final*np.ones(50),'k--',alpha=0.5)
            ax.plot(np.linspace(np.amin(micro)-1,np.amax(micro)+1),(ab_final+Aab_final)*np.ones(50),'k--',alpha=0.3)
            ax.plot(np.linspace(np.amin(micro)-1,np.amax(micro)+1),(ab_final-Aab_final)*np.ones(50),'k--',alpha=0.3)

            ax.set_xlim(np.amin(micro)-1, np.amax(micro)+1)
            ax.set_xlabel(r'$\xi$ [km/s]')
            ax.set_ylabel(r'$\epsilon$ (%s)'%(el_input))  
            
            ax.legend(fancybox=True, framealpha=0.7,loc='upper right')    

            #PLOT 4: AB VS EW PARA CADA LÍNEA
            ax=fig.add_subplot(gs[1:3,int(np.floor(cols/2)):])

            xx=np.linspace(0,np.amax(ew)+1,50)

            ax.plot(xx,ab_final*np.ones(50),'k--',alpha=0.5,
                    label=r'log(%s/H) = %.2f $\pm$ %.2f dex'%(el_input,ab_final,Aab_final))
            ax.plot(xx,ab_final+Aab_final*np.ones(50),'k--',alpha=0.3)
            ax.plot(xx,ab_final-Aab_final*np.ones(50),'k--',alpha=0.3)

            #mask the outliers!

            mask_full = ( (ablinefull < ab_final + 5*Aab_final) & 
                                   (ewfull > ab_final - 5*Aab_final)  )
                                #    & (ewfull < np.median(ewfull) + 5*np.std(ewfull)) & 
                                #    (ewfull > np.median(ewfull) - 5*np.std(ewfull)) )

            for i_cols in range(len(mask_all_cols_full)): 
                
                ax.scatter(ewfull[mask_full & mask_all_cols_full[i_cols]],
                           ablinefull[mask_full &mask_all_cols_full[i_cols]],
                           color='none',
                           edgecolors=linecolor_all_cols[i_cols], #linecolor_all_cols_full[i_cols]
                           marker='o',alpha=0.9,linewidths=3, zorder = 100, s = 100)
                
                #! DIBUJAR LOS OUTLIERS!
                if len(ewfull[~mask_full & mask_all_cols_full[i_cols]])>0:
                    ab_outliers = np.array([np.amax(ablinefull[mask_full]) + 0.2
                                   if ablinefull[~mask_full &mask_all_cols_full[i_cols]][ii]>ab_final
                                   else np.amin(ablinefull[mask_full]) - 0.2
                                   for ii in range(len(ewfull[~mask_full & mask_all_cols_full[i_cols]]))])
                    marker_outliers = np.array(['^' if ablinefull[~mask_full &mask_all_cols_full[i_cols]][ii]>ab_final
                                   else 'v' for ii in range(len(ewfull[~mask_full & mask_all_cols_full[i_cols]]))]) 
                    for ii in range(len(ewfull[~mask_full & mask_all_cols_full[i_cols]])):
                        ax.scatter(ewfull[~mask_full & mask_all_cols_full[i_cols]][ii],ab_outliers[ii],
                                    color='none', edgecolors=linecolor_all_cols[i_cols], #linecolor_all_cols_full[i_cols]
                                    marker=marker_outliers[ii], alpha=0.9,linewidths=3, zorder = 100, s = 100) 

            for i_cols in range(len(mask_all_cols)):

                ax.errorbar(ew[mask_all_cols[i_cols]],abline[mask_all_cols[i_cols]],
                            fmt='o',color= linecolor_all_cols[i_cols],ecolor='k',capsize=5,elinewidth=1,
                            zorder = 500, ms = np.sqrt(100) )  
                
                #,yerr=Aabline[mask_all_cols[i_cols]]
                
                ######## CHECK DIAGONAL UNCERTAINTIES IN EPS(SI) VS EW ############
                
                for i_point in range(len(AabSim1[mask_all_cols[i_cols]])):

                        ax.plot([ (ew[mask_all_cols[i_cols]])[i_point]-np.abs(Aew[mask_all_cols[i_cols]])[i_point],
                                    (ew[mask_all_cols[i_cols]])[i_point]+np.abs(Aew[mask_all_cols[i_cols]])[i_point]  ],
                                   [ (abline[mask_all_cols[i_cols]])[i_point]-np.abs(AabSim1[mask_all_cols[i_cols]])[i_point],
                                    (abline[mask_all_cols[i_cols]])[i_point]+np.abs(AabSim1[mask_all_cols[i_cols]])[i_point] ],
                                   color='k',marker='_' ,lw=0.75,ms=15) #,s=3,linewidths=5
            
            ew_ordered,line_ordered, abline_ordered=zip(*sorted(zip(input_table_full['EW'],input_table_full['Line'], ablinefull)))
            
            for i in range(len(input_table_full)):

                if abline_ordered[i]<=ab_final:
                    if (abline_ordered[i] >= (ab_final - 5*Aab_final)):
                        if line_ordered[i] in input_table['Line']:
                            ax.text(ew_ordered[i],abline_ordered[i],"  [%s]"%i, size=14,color='black', va="top", ha="left")
                        else:
                            ax.text(ew_ordered[i],abline_ordered[i],"  [%s]"%i, size=14,color='grey', va="top", ha="left")
                    else:
                        if line_ordered[i] in input_table['Line']:
                            ax.text(ew_ordered[i], np.amin(ablinefull[mask_full]) - 0.2,
                                    "  [%s]"%i, size=14,color='black', va="top", ha="left")
                        else:
                            ax.text(ew_ordered[i], np.amin(ablinefull[mask_full]) - 0.2,
                                    "  [%s]"%i, size=14,color='grey', va="top", ha="left")                            

                else:
                    if (abline_ordered[i] <= (ab_final + 5*Aab_final)):
                        if line_ordered[i] in input_table['Line']:
                            ax.text(ew_ordered[i],abline_ordered[i],"[%s]  "%i, 
                                    size=14,color='black', va="bottom", ha="right")
                        else:
                            ax.text(ew_ordered[i],abline_ordered[i],"[%s]  "%i, 
                                    size=14,color='grey', va="bottom", ha="right")
                    else:
                        if line_ordered[i] in input_table['Line']:
                            ax.text(ew_ordered[i], np.amax(ablinefull[mask_full]) + 0.2,
                                    "  [%s]"%i, size=14,color='black', va="bottom", ha="right")
                        else:
                            ax.text(ew_ordered[i], np.amax(ablinefull[mask_full]) + 0.2,
                                    "  [%s]"%i, size=14,color='grey', va="bottom", ha="right")                             
            
            # ax.axhspan(np.amin(ab_FW), np.amax(ab_FW), facecolor = "yellow", edgecolor = None, alpha = 0.1, zorder = -999)               
            
            ax.set_ylim(np.nanmin(ablinefull[mask_full])-0.4,
                        np.nanmax([np.nanmax(abline+Aabline),np.nanmax(ablinefull[mask_full])])+0.4)        
            ax.set_xlim(1,
                        np.amax(ewfull[mask_full])+5)           
                
            ax.set_ylabel(r'$\epsilon$ (%s)'%(el_input))
            ax.set_xlabel(r'$EW$ [m$\AA$]')
            ax.legend(fancybox=True, framealpha=0.4,loc='upper right')
            
            # ax=fig.add_subplot(gs[3,-1])
            ax=fig.add_subplot(gs[3,int(np.floor(cols/2)):])
            ax.set_ylim(-len(input_table),1)

            i1 = ceil(len(line_ordered)/4)
            
            for i in range(len(input_table_full)):            
                if i<i1:
                    if line_ordered[i] in input_table['Line']:
                        ax.text(0,-i,"[%i] %s"%(i,line_ordered[i]),fontsize='14',color='black')
                    else:
                        ax.text(0,-i,"[%i] %s"%(i,line_ordered[i]),fontsize='14',color='grey')
                elif i<(i1*2):
                    if line_ordered[i] in input_table['Line']:
                        ax.text(5,-i+i1,"[%i] %s"%(i,line_ordered[i]),fontsize='14', color='black')
                    else:
                        ax.text(5,-i+i1,"[%i] %s"%(i,line_ordered[i]),fontsize='14', color='grey')
                elif i<(i1*3):
                    if line_ordered[i] in input_table['Line']:
                        ax.text(10,-i+(i1*2),"[%i] %s"%(i,line_ordered[i]),fontsize='14', color='black')
                    else:
                        ax.text(10,-i+(i1*2),"[%i] %s"%(i,line_ordered[i]),fontsize='14', color='grey')
                else:
                    if line_ordered[i] in input_table['Line']:
                        ax.text(15,-i+(i1*3),"[%i] %s"%(i,line_ordered[i]),fontsize='14', color='black')
                    else:
                        ax.text(15,-i+(i1*3),"[%i] %s"%(i,line_ordered[i]),fontsize='14', color='grey')                        
            
            ax.set_xlim(-1,18)
            ax.set_ylim(-5.5,1)

            ax.axis('off')
            # ax.set_xticks([])
            # ax.set_yticks([])
            
            #print(star_name,spc_star,(fixed_T),(fixed_g),float(fixed_Q),float(fixed_b/10))
            
            fig.suptitle('%s (%s) - T=%i K, g=%.1f dex, logQ=%.1f, b=%.1f'
                         %(star_name,spc_star,int(fixed_T),float(fixed_g)/100,-float(fixed_Q)/100,float(fixed_b)/10),
                         fontsize='20') 
            if mult_highlight is not None:
                filename=star_name+'_'+el_input+'_abundances_%s'%mult_highlight
            else:
                filename=star_name+'_'+el_input+'_abundances'
            #fig.suptitle('%s (%s) - T=%i, g=%.1f'%(star_name,spc_star,table_Tg['T'][int(Tsel)],table_Tg['g'][int(Tsel)]/100),fontsize='25')
            #filename=star_name+'_'+el_input+'_abundances'+Tsel
            if save=='y':
                if point=='0':
                    # plt.savefig(savepath+filename+'_MAIN_%s.png'% time.strftime('%Y%m%d_%H%M%S'))
                    plt.savefig(savepath+filename+'_MAIN_fixed.png')
                elif point=='1':
                    # plt.savefig(savepath+filename+'_1_%s.png'% time.strftime('%Y%m%d_%H%M%S'))
                    plt.savefig(savepath+filename+'_1_fixed.png')
                  
                elif point=='2':
                    # plt.savefig(savepath+filename+'_2_%s.png'% time.strftime('%Y%m%d_%H%M%S'))                
                    plt.savefig(savepath+filename+'_2_fixed.png')
                
            #plt.savefig(savepath+filename+'.png')
            plt.show(block=False)

        if point=='0':
            final_table=Table([[star_name],[ref_file],[grid],["Unfinished"],[micro_val],[Amicro_val],[Amicro_inf],[Amicro_sup],[ab_final],[Aab_final],[AaO_val],[AaO_m_inf],[AaO_m_sup],[0.0],[0.0],
                               [fixed_T],[int(fixed_Terr)],[-999999.9],[-999.9],[fixed_g_exact],[-999999.9],[-999.9],[fixed_He],[fixed_Q],[fixed_b]]+
                              [[j] for j in ablinefull],
                names=(['ID','Ref_file','Grid','FLAG','micro','Amicro','Amicro_sup','Amicro_inf','ab'+el_input,'Aab'+el_input+'_s','Aab'+el_input+'_m','Aab'+el_input+'_m_inf','Aab'+el_input+'_sup','ab'+el_input+'+','ab'+el_input+'-',
                        'T','Terr','T+','T-','g','g+','g-','He','Q','b']+collines),
                dtype=(['S10',"S200","S50","S10",'float64','float64','float64','float64','float64','float64','float64','float64','float64','float64','float64','float64',
                        'int64','float64','float64','int64','float64','float64','int64','int64','int64'] + ['float64']*len(ablinefull) 
                        ))
            # GUARDAR ESTO EN LA TABLA DE ANALISIS DE TODAS LAS ESTRELLAS CON UN FLAG DE UNFINISHED
            if save=='y':  
                format_i = "csv" if ".csv" in parameterstable else "ecsv" if ".ecsv" in parameterstable else "ascii"
                delimiter_i = "," if ".csv" in parameterstable else "\s" if ".ecsv" in parameterstable else "\s"    

                parameterstable_tab = findtable(parameterstable,format=format_i,path=savepath,delimiter=delimiter_i)
                if parameterstable_tab is not None:
                        if el_input == "Si":            
                            parameterstable_tab = parameterstable_tab[parameterstable_tab["Ref_file"]!=ref_file]
                            parameterstable_tab = vstack([parameterstable_tab,final_table])
                            parameterstable_tab.write(savepath + parameterstable,format='ascii.%s'%format_i,overwrite=True)
                        else:
                            if microreduced == "fixed":

                                final_table_el = final_table[["Ref_file",'ab'+el_input,'Aab'+el_input+'_s','Aab'+el_input+'_m','Aab'+el_input+'_m_inf','Aab'+el_input+'_sup','ab'+el_input+'+','ab'+el_input+'-' ]+collines]
                                colnames_new = ['ab'+el_input,'Aab'+el_input+'_s','Aab'+el_input+'_m','Aab'+el_input+'_m_inf','Aab'+el_input+'_sup','ab'+el_input+'+','ab'+el_input+'-']+collines
                            else:
                                final_table.rename_columns(['micro','Amicro','Amicro_sup','Amicro_inf'],
                                                            ['micro'+el_input,'Amicro'+el_input,'Amicro_sup'+el_input,'Amicro_inf'+el_input])

                                final_table_el = final_table[["Ref_file",'micro'+el_input,'Amicro'+el_input,'Amicro_sup'+el_input,'Amicro_inf'+el_input,
                                                            'ab'+el_input,'Aab'+el_input+'_s','Aab'+el_input+'_m','Aab'+el_input+'_m_inf','Aab'+el_input+'_sup','ab'+el_input+'+','ab'+el_input+'-' ]+collines]
                                colnames_new = ['micro'+el_input,'Amicro'+el_input,'Amicro_sup'+el_input,'Amicro_inf'+el_input,'ab'+el_input,'Aab'+el_input+'_s','Aab'+el_input+'_m','Aab'+el_input+'_m_inf','Aab'+el_input+'_sup','ab'+el_input+'+','ab'+el_input+'-' ]+collines

                            parameterstable_tab_onestar = parameterstable_tab[parameterstable_tab["Ref_file"]==ref_file]
                            colnames_old = parameterstable_tab_onestar.colnames 
                            colnames_old = [coli for coli in colnames_old if coli not in colnames_new]
                            parameterstable_tab_onestar = parameterstable_tab_onestar[colnames_old]

                            final_table = join(parameterstable_tab_onestar,final_table_el,keys = "Ref_file", join_type="outer")

                            parameterstable_tab = parameterstable_tab[parameterstable_tab["Ref_file"]!=ref_file]
                            parameterstable_tab = vstack([parameterstable_tab,final_table])
                            parameterstable_tab.write(savepath + parameterstable,format='ascii.%s'%format_i,overwrite=True)

                else:
                    final_table.write(savepath + parameterstable,format='ascii.%s'%format_i,overwrite=True)

        elif point=='1':
            format_i = "csv" if ".csv" in parameterstable else "ecsv" if ".ecsv" in parameterstable else "ascii"
            delimiter_i = "," if ".csv" in parameterstable else "\s" if ".ecsv" in parameterstable else "\s"    

            parameterstable_tab = findtable(parameterstable,format=format_i,path=savepath,delimiter=delimiter_i)
            if len(parameterstable_tab[parameterstable_tab["Ref_file"] == ref_file])>0:

                parameterstable_tab['ab'+el_input+'+'][parameterstable_tab["Ref_file"] == ref_file] = ab_final
                parameterstable_tab['T+'][parameterstable_tab["Ref_file"] == ref_file] = fixed_T
                # final_table['g-']=parameterstab['g+'][-1] #??
                parameterstable_tab['g+'][parameterstable_tab["Ref_file"] == ref_file]=fixed_g_exact #? que es el que he usado para la abun 
                # GUARDAR ESTO EN LA TABLA DE ANALISIS DE TODAS LAS ESTRELLAS CON UN FLAG DE UNFINISHED
                final_table = parameterstable_tab[parameterstable_tab["Ref_file"] == ref_file]
                if save=='y':
                    parameterstable_tab.write(savepath + parameterstable,format='ascii.%s'%format_i,overwrite=True)
            else:
                print("ERROR: CAN'T FIND INFO ON MAIN PARAMETERS SO CAN'T WRITE THE UNCERTAINTY")
        elif point=='2':
            format_i = "csv" if ".csv" in parameterstable else "ecsv" if ".ecsv" in parameterstable else "ascii"
            delimiter_i = "," if ".csv" in parameterstable else "\s" if ".ecsv" in parameterstable else "\s"    

            parameterstable_tab = findtable(parameterstable,format=format_i,path=savepath,delimiter=delimiter_i)
            if len(parameterstable_tab[parameterstable_tab["Ref_file"] == ref_file])>0:

                parameterstable_tab['ab'+el_input+'-'][parameterstable_tab["Ref_file"] == ref_file] = ab_final
                parameterstable_tab['T-'][parameterstable_tab["Ref_file"] == ref_file] = fixed_T
                # final_table['g-']=parameterstab['g-'][-1] #?? 
                parameterstable_tab['g-'][parameterstable_tab["Ref_file"] == ref_file]=fixed_g_exact #? que es el que he usado para la abun 
                # GUARDAR ESTO EN LA TABLA DE ANALISIS DE TODAS LAS ESTRELLAS CON UN FLAG DE UNFINISHED
                final_table = parameterstable_tab[parameterstable_tab["Ref_file"] == ref_file]
                if save=='y':
                    parameterstable_tab.write(savepath + parameterstable,format='ascii.%s'%format_i,overwrite=True)
            else:
                print("ERROR: CAN'T FIND INFO ON MAIN PARAMETERS SO CAN'T WRITE THE UNCERTAINTY")      

    return final_table 
    
#%%

def abundances(ref_file,spc_star,el_input,star_db,star_measure,FW_EW_measure,
               inputparameters,parameterstable,linetable,savepath,
               linesvec = ['SiIII4552','SiIII4567','SiIII4574'],
               mult_use="all", mult_highlight=None,interpolate=True,
               ewcut=5, redo_count=0,tol=50, microreduced="montecarlo", 
               auto="n", show="y", unc="y"):

    star_measuretab=findtable(star_measure)
    FW_EW_measuretab=findtable(FW_EW_measure)

    star_name = star_measuretab["ID"][star_measuretab["Ref_file"] == ref_file][0]
    linetab=findtable(linetable,format='csv',path=savepath)
    linelabels=linetab['label']

    if isinstance(mult_use, (list, tuple, np.ndarray)):
        mult_use = np.array(["%s"%(i) for i in mult_use])
        mu = 'y'
    elif mult_use == "all":
        mu = 'n'
    else:
        mu = 'n'

    if ref_file not in linetab.colnames:

        ind = (linetab.colnames).index("flag")
        linetab.add_column(np.zeros(len(linelabels),dtype=int),name=ref_file,index = ind+1)

        for lineii in linetab['label']:
            if "RV0" in star_measuretab.colnames:
                if (star_measuretab['EW_'+lineii][star_measuretab['Ref_file']==ref_file][0]<ewcut 
                    or (np.ma.getmask(star_measuretab['EW_'+lineii][star_measuretab['Ref_file']==ref_file][0]) )
                    or (np.abs(star_measuretab['RV0'][star_measuretab['Ref_file']==ref_file][0]-
                            star_measuretab['RV_'+lineii][star_measuretab['Ref_file']==ref_file][0])>0.5*tol)):
                    linetab[ref_file][linetab['label']==lineii]=1
            else:
                if (star_measuretab['EW_'+lineii][star_measuretab['Ref_file']==ref_file][0]<ewcut 
                    or np.ma.getmask(star_measuretab['EW_'+lineii][star_measuretab['Ref_file'][0]==ref_file]) ):
                    linetab[ref_file][linetab['label']==lineii][0]=1
                    
            if linetab['flag'][linetab['label']==lineii] not in ['green']:
                linetab[ref_file][linetab['label']==lineii]=2
                
            if mu == 'y': 
                use = np.any(np.array([True if i in mult_use else False for i in linetab['mult_FW'][linetab['label']==lineii][0].split("+")]))
                if use == False:
                    linetab[ref_file][linetab['label']==lineii]=3

        linetab.write(savepath+linetable, format='ascii.csv', overwrite=True, delimiter=(','))     

    # stop=input('\nPlease check the line fitting for '+star_name+' and mark any bad fit as 1 in the line table, ' +linetable)
    
    parameterstab2=findtable(inputparameters,format='ecsv',path=savepath,delimiter='\s')

    parameterstab2 = parameterstab2[parameterstab2["Ref_file"] == ref_file]
    
    if "micro"+el_input in parameterstab2.colnames:
        micro_center = parameterstab2['micro'+el_input][-1]
    else:
        micro_center = parameterstab2['micro'][-1]

    micro_center=parameterstab2['micro'][-1]
    micro_selection='center'
    microlist=sorted(set(FW_EW_measuretab['micro']))
    micro_old=microlist[min(range(len(microlist)),key=lambda i : abs(int(microlist[i])-float(micro_center)))]

    stop='y'
    microstop='y'

    while microstop=='y':
        
        while stop=='y':
            
            parameterstab2=findtable(inputparameters,format='ecsv',path=savepath,delimiter='\s')
            parameterstab2 = parameterstab2[parameterstab2["Ref_file"] == ref_file]
            #GET THE NEAREST ONES (T AND G PARTICULARLY)

            fixed_T=parameterstab2['T'][-1]
            fixed_g=parameterstab2['g'][-1]
            fixed_He=parameterstab2['He'][-1]
            fixed_Q=parameterstab2['Q'][-1]
            fixed_b=parameterstab2['b'][-1]

            print("Teff = %d K; log g = %.2f dex"%(int(fixed_T), int(fixed_g)/100))

            garrFW=sorted(set(FW_EW_measuretab['g']))
            fixed_g=garrFW[min(range(len(garrFW)),key=lambda i : abs(int(garrFW[i])-int(fixed_g)))] #NEAREST
            
            mask1=((FW_EW_measuretab['g']=="%d"%int(fixed_g)) &
                              (FW_EW_measuretab['He']=="%d"%int(fixed_He)) & (FW_EW_measuretab['Q']=="%d"%int(fixed_Q)) & 
                              (FW_EW_measuretab['b']=="%d"%int(fixed_b)) 
                              )
            
            if np.sum(mask1)<1:
                print(fixed_g, fixed_He, fixed_Q, fixed_b) 
                print(sorted(set(FW_EW_measuretab['g'])))
                print(sorted(set(FW_EW_measuretab['He'])))
                print(sorted(set(FW_EW_measuretab['Q'])))
                print(sorted(set(FW_EW_measuretab['b'])))

                print("Please review the parameters g, He, Q, b, and choose a combination that exists in the FASTWIND"
                      +"grid that is currently being used.")
                break 

            if interpolate == True: 

                fixed_T=str('%d'%(int(fixed_T)/100))
                maskT=((FW_EW_measuretab['T']==(fixed_T)) & (mask1 ))
                        
                if np.sum(maskT)<1:
                
                    if show == "y": print("No information for T=%d in the FASTWIND tables. Interpolating... \n"%(int(fixed_T)))
                    
                    FASTWIND_tab_interpolate(fixed_T,FW_EW_measure,fixed_g,fixed_He,fixed_Q,fixed_b,el_input,n=4)
                    
                    FW_EW_measuretab=findtable(FW_EW_measure)

            else: 
                
                fixed_T =  float(fixed_T)/100
                TarrFW=sorted(set(FW_EW_measuretab[mask1]['T']))
                fixed_T_interpolate=TarrFW[min(range(len(TarrFW)),key=lambda i : abs(int(TarrFW[i])-int(fixed_T)))] #NEAREST

                if np.abs(int(fixed_T_interpolate)-int(fixed_T))>15:
                    print("ERROR: GRID IS NOT BIG ENOUGH TO USE T=%d K, logg=%.2f dex" %(int(fixed_T)*100, int(fixed_g)/100))
                    return None
                
                fixed_T = fixed_T_interpolate
                

                parameterstab2.add_row(parameterstab2[-1])
                parameterstab2["Ref_file"][-1] = ref_file
                parameterstab2['T'][-1] = int(fixed_T)*100
                
                parameterstab2.write(savepath+inputparameters,  #! parameterstab
                                 format='ascii.ecsv',overwrite=True) 
            
            # WE CALCULATE THE CENTRAL ABUNDANCE + MICROTURBULENCE

            print("Calculating the curve of growth plots for each line...")
            final_tables,linesvec,EWvec,abvec=tables_curve_of_growth(ref_file,el_input,spc_star,
                                            star_measure,FW_EW_measure,
                                            inputparameters,linetable,
                                            micro_selection,savepath,
                                            ewcut=ewcut,tol=tol,micro_center=micro_center,linesvec=linesvec,show=show)
            if auto == "n":
                stop = input("If you want to change the lines used in file %s, press 'y'. \n"%linetable +
                         "Press ENTER when you're happy with the line selection." )
            
            stop='n'

            final_table = abundances_growthcurve(star_name,ref_file, spc_star, final_tables, inputparameters, 
                                                 el_input, linetable,FW_EW_measure,parameterstable,savepath, microreduced,
                                                 linesvec,EWvec,abvec, mult_highlight=mult_highlight,show=show)

            if "micro%s"%el_input in final_table.colnames:
                micro_center=final_table['micro%s'%el_input][-1]
            else: 
                micro_center=final_table['micro'][-1]

            if auto=="n":
            
                microstop=input('Do you want to work with a different range of microturbulence? [y/n] - Press enter to automatically center the microturbulence: ') 
                if microstop=='y':
                    stop='y'
                    print('\nCurrent microturbulence range: ', micro_selection)
                    micro_selection=input('Choose the range of microturbulence 0-30 km/h from [down/mid/up/center]: ')
                    while micro_selection not in ['down','mid','up','center']:
                        print('Please choose an available option:')
                        micro_selection=input('Choose the range of microturbulence 0-30 km/h from [down/mid/up/center]: ')

                # Tstop=input('Do you want to work with a different temperature? [y/n]: ') 
                # if Tstop=='y':
                #     stop='y'
                    
                #     parameters=findtable(star_name+'_input_parameters.ecsv', #! parameterstab
                #                         path=savepath,format='ecsv',delimiter='\s')

                #     Told_selection=parameters['T'][-1]
                #     print('\nCurrent Temperature : ', Told_selection)
                #     Tnew_selection=input('Choose how much you want to change the Temperature (e.g. : +500 / -500 -> ')
                    
                #     if '+' in Tnew_selection:
                #         Tnew_selection = Tnew_selection.replace('+','')
                #         Tnew_selection=int(Tnew_selection)

                #         parameters.add_row(parameters[-1])
                #         parameters["Ref_file"][-1] = ref_file
                #         parameters['T'][-1]='%d'%(int(parameters['T'][-2])+Tnew_selection)
                #         parameters['T+'][-1]='%d'%(int(parameters['T+'][-2])+Tnew_selection)
                #         parameters['T-'][-1]='%d'%(int(parameters['T-'][-2])+Tnew_selection)
                        
                #         parameters.write(savepath+star_name+'_input_parameters.ecsv',  #! parameterstab
                #                     format='ascii.ecsv',overwrite=True)                
                        
                #     elif '-' in Tnew_selection:
                #         print("RECONGNIZING THE -")
                        
                #         Tnew_selection = Tnew_selection.replace('-','')
                #         Tnew_selection=int(Tnew_selection) 

                #         parameters.add_row(parameters[-1])
                #         parameters["Ref_file"][-1] = ref_file
                #         print("--------------------------")
                #         print(int(parameters['T'][-2]))
                #         print(Tnew_selection)
                #         print(int(parameters['T'][-2]) - Tnew_selection)
                #         print('%d'%(int(parameters['T'][-2]) - Tnew_selection))
                #         print("--------------------------")
                #         parameters['T'][-1]='%d'%(int(parameters['T'][-2]) - Tnew_selection)
                #         parameters['T+'][-1]='%d'%(int(parameters['T+'][-2]) - Tnew_selection)
                #         parameters['T-'][-1]='%d'%(int(parameters['T-'][-2]) - Tnew_selection)

                #         parameters.write(savepath+star_name+'_input_parameters.ecsv',
                #                         format='ascii.ecsv',overwrite=True)
                #     print("New temperature: ",parameters['T'][-1] )
                fullstop=input('If you want to change any other parameters, change the file '+
                            star_name+'_'+el_input+'_parameters.ecsv and press [y]. Press any key to continue. ' )
                
                if fullstop=='y':
                    stop='y'
            
            elif auto=="y":
                micro_selection="center"

            else:
                print("'auto' keyword not 'y' or 'n': assuming 'y' and continuing in auto mode.")
                
        micro_new=microlist[min(range(len(microlist)),key=lambda i : abs(int(microlist[i])-float(micro_center)))]
        # if show == "y":
            # print('Old micro %s, New micro %s'%(micro_old , micro_new))
            # print('Currently on iteration number %d \n'%(redo_count))
        
        if ((int(micro_old)==int(micro_new)) or (redo_count>3)):
            print("Stopping the iterations\n")
            microstop = 'n'
            stop = 'n'
            redo_count=redo_count+1

        else:
            print("Iterating again")
            microstop = 'y'
            stop = 'y'
            micro_old = micro_new 

        ######################################################################
        ########### CALCULATING UNCERTAINTIES FOR T,G,Q ######################
        ######################################################################

        if ((microstop == 'n') & (stop=='n') & (unc == "y")):
            print('Repeating the analysis for different T, g values to get the uncertainty...')

            # UPPER T, G UNCERTAINTY 
  
            point='1'
            fixed_T=parameterstab2['T+'][-1]
            fixed_g=parameterstab2['g+'][-1]
            fixed_g=garrFW[min(range(len(garrFW)),key=lambda i : abs(int(garrFW[i])-int(fixed_g)))] #NEAREST
            # garrFW.remove(fixed_g)    
            
            mask1=((FW_EW_measuretab['g']=="%d"%int(fixed_g)) &
                            (FW_EW_measuretab['He']=="%d"%int(fixed_He)) & (FW_EW_measuretab['Q']=="%d"%int(fixed_Q)) & 
                            (FW_EW_measuretab['b']=="%d"%int(fixed_b) )
                            )
            
            if np.sum(mask1)<1:
                input("Please review the parameters g, He, Q, b, and choose a combination that exists in the FASTWIND"+
                    "grid that is currently being used.")
            
            fixed_T=str('%d'%(int(fixed_T)/100))
            maskT=((FW_EW_measuretab['T']==(fixed_T)) & (mask1 ))

            if interpolate == True: 
 
                if np.sum(maskT)<1:
                
                    if show == "y": print("No information for T=%d in the FASTWIND tables. Interpolating... \n"%(int(fixed_T)))
                    
                    FASTWIND_tab_interpolate(fixed_T,FW_EW_measure,fixed_g,fixed_He,fixed_Q,fixed_b,el_input,n=4)
                    
                    FW_EW_measuretab=findtable(FW_EW_measure)
                    mask1=((FW_EW_measuretab['g']=="%d"%int(fixed_g)) &
                            (FW_EW_measuretab['He']=="%d"%int(fixed_He)) & (FW_EW_measuretab['Q']=="%d"%int(fixed_Q)) & 
                            (FW_EW_measuretab['b']=="%d"%int(fixed_b) )   )                  
                    maskT=((FW_EW_measuretab['T']==(fixed_T)) & (mask1 ))                        
            
            else: 
                
                TarrFW=sorted(set(FW_EW_measuretab[mask1]['T']))
                fixed_T=TarrFW[min(range(len(TarrFW)),key=lambda i : abs(int(TarrFW[i])-int(fixed_T)))] #NEAREST

                maskT=((FW_EW_measuretab['T']==(fixed_T)) & (mask1 ))
                       
            final_tables,linesvec,EWvec,abvec=tables_curve_of_growth(ref_file,el_input,spc_star,star_measure,
                                            FW_EW_measure,inputparameters,linetable,
                                            micro_selection='center',savepath=savepath,ewcut=ewcut,
                                            micro_center=micro_new,linesvec=linesvec, do_pdf=False,
                                            point=point, show="n")
                                            
            final_table = abundances_growthcurve(star_name,ref_file, spc_star, final_tables, inputparameters, 
                                                  el_input, linetable,FW_EW_measure,parameterstable,
                                                savepath, 'fixed',linesvec,EWvec,abvec,mult_highlight=mult_highlight, 
                                                do_pdf=True,point=point, show="n")  
            
            # LOWER T, G UNCERTAINTY 

            point='2'
            fixed_T=parameterstab2['T-'][-1]
            fixed_g=parameterstab2['g-'][-1]
            fixed_g=garrFW[min(range(len(garrFW)),key=lambda i : abs(int(garrFW[i])-int(fixed_g)))] #NEAREST
            # garrFW.remove(fixed_g)            
            mask1=((FW_EW_measuretab['g']=="%d"%int(fixed_g)) &
                            (FW_EW_measuretab['He']=="%d"%int(fixed_He)) & (FW_EW_measuretab['Q']=="%d"%int(fixed_Q)) & 
                            (FW_EW_measuretab['b']=="%d"%int(fixed_b) )
                            )
            
            if np.sum(mask1)<1:
                input("Please review the parameters g, He, Q, b, and choose a combination that exists in the FASTWIND"+
                    "grid that is currently being used.")
            
            fixed_T=str('%d'%(int(fixed_T)/100))
            maskT=((FW_EW_measuretab['T']==(fixed_T)) & (mask1 ))

            if interpolate == True: 
 
                if np.sum(maskT)<1:
                
                    if show == "y": print("No information for T=%d in the FASTWIND tables. Interpolating... \n"%(int(fixed_T)))
                    
                    FASTWIND_tab_interpolate(fixed_T,FW_EW_measure,fixed_g,fixed_He,fixed_Q,fixed_b,el_input,n=4)
                    FW_EW_measuretab=findtable(FW_EW_measure)
                    fixed_T=str('%d'%(int(fixed_T)/100))
                    maskT=((FW_EW_measuretab['T']==(fixed_T)) & (FW_EW_measuretab['g']=="%d"%int(fixed_g)) &
                            (FW_EW_measuretab['He']=="%d"%int(fixed_He)) & (FW_EW_measuretab['Q']=="%d"%int(fixed_Q)) & 
                            (FW_EW_measuretab['b']=="%d"%int(fixed_b)) )
            
            else: 
                
                TarrFW=sorted(set(FW_EW_measuretab[mask1]['T']))
                fixed_T=TarrFW[min(range(len(TarrFW)),key=lambda i : abs(int(TarrFW[i])-int(fixed_T)))] #NEAREST

                maskT=((FW_EW_measuretab['T']==(fixed_T)) & (FW_EW_measuretab['g']=="%d"%int(fixed_g)) &
                            (FW_EW_measuretab['He']=="%d"%int(fixed_He)) & (FW_EW_measuretab['Q']=="%d"%int(fixed_Q)) & 
                            (FW_EW_measuretab['b']=="%d"%int(fixed_b)))
                       
            final_tables,linesvec,EWvec,abvec=tables_curve_of_growth(ref_file,el_input,spc_star,star_measure,
                                            FW_EW_measure,inputparameters,linetable,
                                            micro_selection='center',savepath=savepath,ewcut=ewcut,
                                            micro_center=micro_new,linesvec=linesvec, do_pdf=False,
                                            point=point, show="n")
                                            
            final_table = abundances_growthcurve(star_name,ref_file, spc_star, final_tables, inputparameters, 
                                                  el_input, linetable,FW_EW_measure,parameterstable,
                                                savepath, 'fixed',linesvec,EWvec,abvec,mult_highlight=mult_highlight, 
                                                do_pdf=True,point=point, show="n")  
    return final_table 


def parameters(ref_file,el_input,star_db,star_measure,finallinelist,linestab,lines,
               FW_EW_measure,parameterstable,savepath, auto = "n", show="y", unc="y"):
    
    '''
    star_name: str

    el_input: str

    star_db: str

    star_measure: str

    finallinelist: str

    linestab: str

    lines: str

    FW_EW_measure: str

    grid_name: str

    parameterstable: str

    savepath: str

    auto:str 
        "y"/"n"

    show:str
        if show="n", does not show any outputs, just saves everything. 

    '''

    print("\n-----------------\nSTARTING THE ANALYSIS\n-----------------\n")

    # FIRST, WE DEFINE ALL THE RELEVANT PARAMETERS FOR THE ANALYSIS. 


    star_name=findtable(star_db)['ID'][findtable(star_db)['Ref_file']==ref_file][0]  

    if "Class." in   findtable(star_db).colnames:
        spc_star=findtable(star_db)['Class.'][findtable(star_db)['Ref_file']==ref_file][0]
    elif "SpType" in findtable(star_db).colnames:
        spc_star=findtable(star_db)['SpType'][findtable(star_db)['Ref_file']==ref_file][0]
    else: 
        spc_star = "-"

    print("%s (%s)"%(star_name, spc_star))

    inputparameters=star_name+'_input_parameters.ecsv'
    
    FW_EW_measure_tab=findtable(FW_EW_measure)
    grid_name = FW_EW_measure_tab["Grid_name"][0]

    inputparameters=star_name+'_input_parameters.ecsv'
    FW_EW_measure_tab=findtable(FW_EW_measure)

    input_parameters=findtable(inputparameters,path=savepath,format='ecsv',delimiter='\s')
    if input_parameters is not None:
        nn_0 = len(input_parameters)
    else:
        nn_0 = 0

    # INITIAL PARAMETERS OF THE STAR (FROM PREVIOUS ANALYSIS OR STANDARD PARAMETERS)


    format_i = "csv" if ".csv" in parameterstable else "ecsv" if ".ecsv" in parameterstable else "ascii"
    delimiter_i = "," if ".csv" in parameterstable else "\s" if ".ecsv" in parameterstable else "\s"    
    final_parameters=findtable(parameterstable,path=savepath,format=format_i,delimiter=delimiter_i)
    
    if final_parameters is not None:
        if ref_file in final_parameters["Ref_file"]:
            final_parameters = final_parameters[final_parameters["Ref_file"] == ref_file]
            final_Q=final_parameters['Q'][-1]
            final_He=final_parameters['He'][-1]
            final_b = final_parameters['b'][-1]
            final_micro=final_parameters['micro'][-1]
            # final_abun = final_parameters['ab%s'%el_input][-1]
            # ablist=sorted(set(FW_EW_measure_tab['abun_%s'%el_input]))
            # final_abun = ablist[min(range(len(ablist)),key=lambda i : abs(float(ablist[i])-float(final_abun)))]
            final_abun='centr'
            microlist=sorted(set(FW_EW_measure_tab['micro']))
            final_micro=microlist[min(range(len(microlist)),key=lambda i : abs(int(microlist[i])-float(final_micro)))]
                    
        else:
            final_Q='1400'
            final_He='10'
            final_b="10"
            final_abun='centr'   
            if 'bsg' in FW_EW_measure:
                final_micro='020'
            elif 'bms' in FW_EW_measure:
                final_micro='005'
            else:
                final_micro='001' 
    else:
        final_Q='1400'
        final_He='10'
        final_b="10"
        final_abun='centr'   
        if 'bsg' in FW_EW_measure:
            final_micro='020'
        elif 'bms' in FW_EW_measure:
            final_micro='005'
        else:
            final_micro='001'
    
    n = 0 

    if show=="y": print('Starting with initial micro of %s km/s'%final_micro)
    
    if show=="y": print(ref_file)
    
    if '.fits' in ref_file:
        orig='IACOB'
    else:
        orig='txt'    

    linetab=findtable(linestab,format='csv',path=savepath)
    if linetab is None: 
        raise FileNotFoundError("Could not find the file %s in path %s"%(linestab, savepath))    

    # DEFINIMOS LAS LÍNEAS CON LAS QUE CALCULAREMOS LA MICRO! 

    if el_input == 'Si':
        linesvec = ['SiIII4552','SiIII4567','SiIII4574'] 
    else:
        #! Usar todas las líneas ! 
        linesvec = ["all"] # USAMOS TODAS LAS LÍNEAS PARA CALCULAR LA MICRO

    print("\n-----------------\nPLOTTING SOME RELEVANT SPECTRAL LINES\n-----------------\n")        
    
    plot_lines(ref_file,spc_star,orig,star_measure,'lines_T_g.lst',savepath)

    if auto == "n": stop = input("Press return to continue with the analysis.")

    plt.close("all")

    redo='y'
    repeat = "y"
    redo_count=0

    # FIRST APPROXIMATION OF T, G, Q, BETA USING BETA = 1, LOGQ = -14


    while n < 1 : 

        # first approximation -> beta = 1 and log Q = -14

        fixed_He=final_He

        beta_fixed = True
        fixed_b="10"
        fixed_Q="1400"
        fixed_micro=final_micro
        fixed_abun=final_abun

        print("\n-----------------\nIONIZATION EQUILIBRIUM\n-----------------\n")
        table_calcTg, flag, T_lines = calc_Tg(ref_file,spc_star,orig,el_input,star_measure,FW_EW_measure,lines,star_db,savepath,
                  fixed_Q,fixed_He,fixed_b,fixed_micro,fixed_abun,show=show)
        
        if auto == "n": stop = input("Press return to continue with the analysis.")

        plt.close("all")

        print("\n-----------------\nBALMER LINES\n-----------------\n")
        tab_Tgfinal, tab_inputparameters, best_QQ, best_bb, g_lines, FLAG_EDGE = best_Tg(ref_file,spc_star,orig,table_calcTg,star_measure,FW_EW_measure,
                                                      star_db,parameterstable,savepath,fixed_Q,fixed_He,fixed_b,fixed_micro,fixed_abun,
                                                      beta_fixed = beta_fixed,show=show)
        
        if auto == "n": stop = input("Press return to continue with the analysis.")

        plt.close("all")
        
        # this gives us an approximation on T, g, Q and beta.

        if int(best_QQ) >= 1350 :  
            beta_fixed = True
            final_b = "10"
            final_Q = best_QQ

        else: 
            beta_fixed = True

            # WE FIX THE BETA PARAMETER FOLLOWING CASTRO, URBANEJA ET AL. (2012)
            T = int(tab_Tgfinal["T"][0])
            if T<10000: 
                best_bb = 3.6
            elif T>30000:
                best_bb = 0.7
            else: 
                best_bb = 5 - 1.4 * T*1e-4 

            betalist=sorted(set(FW_EW_measure_tab['b']))
            final_b=betalist[min(range(len(betalist)),key=lambda i : abs(int(betalist[i])-float(best_bb)*10 ))]

        if (fixed_Q == final_Q) & (fixed_b == final_b): repeat ="n"

        n = n + 1

    if show == "y": 
        print("First approximation of T, g, Q and b:")
        print("Teff = %d K , logg =%.1f dex, log Q = %d dex, beta = %d "%(int(tab_Tgfinal["T"][0]), float(tab_Tgfinal["g"][0]), 
                                                                          int(final_Q), int(final_b) ))
        
    plt.close("all")

    while redo == "y": 
        fixed_Q=final_Q
        fixed_He=final_He
        fixed_b=final_b
        fixed_micro=final_micro
        fixed_abun=final_abun        

        count = 0

        while repeat == "y": 
            # we don't need to do it again if the wind parameters didn't change. 
            count=count+1


            print("\n-----------------\nIONIZATION EQUILIBRIUM -- TEFF CALCULATION\n-----------------\n")
            table_calcTg, flag, T_lines = calc_Tg(ref_file,spc_star,orig,el_input,star_measure,FW_EW_measure,lines,star_db,savepath,
                    fixed_Q,fixed_He,fixed_b,fixed_micro,fixed_abun,show=show)
            
            if auto == "n": stop = input("Press return to continue with the analysis.")

            plt.close("all")

            print("\n-----------------\nBALMER LINES -- LOG G CALCULATION\n-----------------\n")
            tab_Tgfinal, tab_inputparameters, best_QQ, best_bb, g_lines, FLAG_EDGE = best_Tg(ref_file,spc_star,orig,table_calcTg,star_measure,FW_EW_measure,
                                                        star_db,parameterstable,savepath,fixed_Q,fixed_He,fixed_b,fixed_micro,fixed_abun,
                                                        beta_fixed = beta_fixed,show=show)
            
            if auto == "n": stop = input("Press return to continue with the analysis.")

            plt.close("all")

            
            if final_b == "10" :  
                beta_fixed = True
                final_Q = best_QQ

            else: 
                beta_fixed = True

                # FIJAMOS BETA SEGÚN CASTRO, URBANEJA ET AL. (2012)
                T = int(tab_Tgfinal["T"][0])
                if T<10000: 
                    best_bb = 3.6
                elif T>30000:
                    best_bb = 0.7
                else: 
                    best_bb = 5 - 1.4 * T*1e-4    

                betalist=sorted(set(FW_EW_measure_tab['b']))
                final_b=betalist[min(range(len(betalist)),key=lambda i : abs(int(betalist[i])-float(best_bb)*10 ))]                             

            if (fixed_Q == final_Q) & (fixed_b == final_b): 
                repeat ="n"
            if count>=3: 
                repeat = "n"
                
        if FLAG_EDGE == "EDGE":
            flag = FLAG_EDGE     
        
        redo='n'
        print("\n-----------------\nABUNDANCES\n-----------------\n")
        abundances(ref_file,spc_star,el_input,star_db,star_measure,FW_EW_measure,
                        inputparameters,parameterstable,linestab,savepath,
                        linesvec=linesvec,ewcut=5,redo_count=redo_count, 
                        auto=auto, show=show, unc=unc)
        
        plt.close("all")

        format_i = "csv" if ".csv" in parameterstable else "ecsv" if ".ecsv" in parameterstable else "ascii"
        delimiter_i = "," if ".csv" in parameterstable else "\s" if ".ecsv" in parameterstable else "\s"    
        final_parameters=findtable(parameterstable,path=savepath,format=format_i,delimiter=delimiter_i)
        final_parameters = final_parameters[final_parameters["Ref_file"] == ref_file]

        #!!! añadir referencia al ref_file
        final_Q=str(int(final_parameters['Q'][-1]))
        final_He=str(int(final_parameters['He'][-1]))
        final_b=str(int(final_parameters['b'][-1]))
        final_micro=final_parameters['micro'][-1]
        microlist=sorted(set(FW_EW_measure_tab['micro']))
        nearest=microlist[min(range(len(microlist)),key=lambda i : abs(int(microlist[i])-float(final_micro)))]
        final_micro=nearest
        final_abun=final_parameters['ab'+el_input][-1]
        abunlist=sorted(set(FW_EW_measure_tab['abun_'+el_input]))
        nearest=abunlist[min(range(len(abunlist)),key=lambda i : abs(float(abunlist[i])-float(final_abun)))]
        final_abun=nearest
        final_g=final_parameters['g'][-1]
        final_T=final_parameters['T'][-1]

        if final_micro!=fixed_micro: # 

            input_parameters=findtable(inputparameters,path=savepath,format='ecsv',delimiter='\s')
            nn_i = len(input_parameters)

            if np.abs(nn_i - nn_0)>10:
                redo = "n"
            else:
                if show=="y":
                    print("Repeating the analysis since the begining with the new parameters...")
                    print("New parameters:", final_Q, final_He, final_b, final_micro)
                    print("Old parameters", fixed_Q, fixed_He, fixed_b, fixed_micro)

                #CHANGE INPUT PARAMETERS:
            
                input_parameters.add_row(input_parameters[-1])
                #!!! añadir referencia al ref_file
                input_parameters['Q'][-1]=final_Q
                input_parameters['b'][-1]=final_b
                input_parameters['micro'][-1]=final_micro
                input_parameters['abun'][-1]=final_abun

                input_parameters.write(savepath+inputparameters, format='ascii.ecsv',
                                    overwrite=True)
                redo='y'
                repeat = "y"
                redo_count=redo_count+1

        else:
            redo='n'
    
    # COMMENTED: FUNCTION BEST_HE IS NOT READY YET

    # final_He = best_He(ref_file, spc_star, orig, star_measure, FW_EW_measure, star_db, parameterstable, savepath)


    # FORMATTING THE TABLE WITH THE RESULTS: 

    format_i = "csv" if ".csv" in parameterstable else "ecsv" if ".ecsv" in parameterstable else "ascii"
    delimiter_i = "," if ".csv" in parameterstable else "\s" if ".ecsv" in parameterstable else "\s"    
    final_parameters=findtable(parameterstable,path=savepath,format=format_i,delimiter=delimiter_i)
    # final_parameters = final_parameters[final_parameters["Ref_file"] == ref_file][-1]

    if "T_lines" not in final_parameters.colnames:
        final_parameters["T_lines"]=["-"*30]
    if "g_lines" not in final_parameters.colnames:
        final_parameters["g_lines"]=["-"*30]
    if "FLAG" not in final_parameters.colnames:
        final_parameters["FLAG"]=["-"*10]

    final_parameters['T_lines']=final_parameters['T_lines'].astype('30str')
    final_parameters['g_lines']=final_parameters['g_lines'].astype('30str')
    final_parameters['FLAG']=final_parameters['FLAG'].astype('10str')

    final_parameters["T_lines"][final_parameters["Ref_file"] == ref_file] = T_lines
    final_parameters["g_lines"][final_parameters["Ref_file"] == ref_file] = str(g_lines)
    final_parameters["FLAG"][final_parameters["Ref_file"] == ref_file] = flag

    final_parameters.write(savepath+parameterstable,format='ascii.%s'%format_i, overwrite=True)

    final_parameters = final_parameters[final_parameters["Ref_file"] == ref_file][-1]
    print("Final results: T=%s, g=%s, Q=%s, b=%s, He=%s, micro=%s, abun=%s"%
          (final_parameters['T'],final_parameters['g'],final_parameters['Q'],
           final_parameters['b'],final_parameters['He'],final_parameters['micro'], 
           final_parameters['abSi']))
    
    clean="y"
    if clean=='y':
        if show=="y": print('Cleaning the FW table, removing the interpolation data...')
        FW_EW_measure_tab=FW_EW_measure_tab[FW_EW_measure_tab['Grid_name']!= "symetric-interp"]
        FW_EW_measure_tab=FW_EW_measure_tab[FW_EW_measure_tab['Grid_name']!= "asymetric-interp"]
        FW_EW_measure_tab.write(maindir+'tables/'+FW_EW_measure,overwrite=True)

    finallinelist='lines_finalplot.lst'
    print("\n-----------------\nFINAL PLOT\n-----------------\n")
    plot_best_model(ref_file,el_input,star_db,star_measure,
                    finallinelist,FW_EW_measure,parameterstable,savepath)   

    if auto == "n": stop = input("Press return to finish with the analysis.")

    plt.close("all")

    print("\n-----------------\nDONE!!\n-----------------\n")
    print('DONE! FLAG = %s \n\n'%flag)

    return final_parameters

def make_summary_pdf(ref_file,savepath):

    print("Making a pdf summary of all the plots in the analysis...")

    star_name = ref_file.split("_")[0]
    
    import os
    import glob 

    print()
    
    if os.path.isfile(savepath+ref_file.split(".")[0]+'_bestfit.png'):

        from matplotlib.backends.backend_pdf import PdfPages
        import matplotlib.image as mpimg

        output=savepath+ref_file.split(".")[0]+'_summaryplots.pdf'

        with PdfPages(output) as pdf:

            bar = pb.ProgressBar(maxval=8,widgets=[pb.Bar('=','[',']'),' ',pb.Percentage()])
            bar.start()                  

            fig, ax = plt.subplots(figsize=(8.27,11.69))
            img =  mpimg.imread(savepath + ref_file.split(".")[0]+'_bestfit.png')
            h, w = img.shape[:2]
            aspect = w / h
            img_h = 10 / 2.54  # 10cm in inches
            img_w = img_h * aspect
            ax.imshow(img)
            ax.axis('off')
            fig.subplots_adjust(0, 0, 1, 1)
            pdf.savefig(fig)
            plt.close(fig)      
            bar.update(0)


            fig, ax = plt.subplots(figsize=(8.27,11.69))
            img =  mpimg.imread(savepath + star_name+'EW2EW1_vs_Teff.png')
            h, w = img.shape[:2]
            aspect = w / h
            img_h = 10 / 2.54  # 10cm in inches
            img_w = img_h * aspect
            ax.imshow(img)
            ax.axis('off')
            fig.subplots_adjust(0, 0, 1, 1)
            pdf.savefig(fig)
            plt.close(fig)      
            bar.update(1)    

            fig, ax = plt.subplots(figsize=(8.27,11.69))
            img =  mpimg.imread(savepath +star_name+'balmer.png')
            h, w = img.shape[:2]
            aspect = w / h
            img_h = 10 / 2.54  # 10cm in inches
            img_w = img_h * aspect
            ax.imshow(img)
            ax.axis('off')
            fig.subplots_adjust(0, 0, 1, 1)
            pdf.savefig(fig)
            plt.close(fig)      
            bar.update(2)

            fig, ax = plt.subplots(figsize=(8.27,11.69))
            img =  mpimg.imread(savepath +star_name+'chi2.png')
            h, w = img.shape[:2]
            aspect = w / h
            img_h = 10 / 2.54  # 10cm in inches
            img_w = img_h * aspect
            ax.imshow(img)
            ax.axis('off')
            fig.subplots_adjust(0, 0, 1, 1)
            pdf.savefig(fig)
            plt.close(fig)      
            bar.update(3) 


            fig, ax = plt.subplots(figsize=(8.27,11.69))
            img =  mpimg.imread(savepath + star_name+'_Q.png')
            h, w = img.shape[:2]
            aspect = w / h
            img_h = 10 / 2.54  # 10cm in inches
            img_w = img_h * aspect
            ax.imshow(img)
            ax.axis('off')
            fig.subplots_adjust(0, 0, 1, 1)
            pdf.savefig(fig)
            plt.close(fig)      
            bar.update(4)     

            fig, ax = plt.subplots(figsize=(8.27,11.69))
            img =  mpimg.imread(savepath +'abundances_'+star_name+'_Si.png')
            h, w = img.shape[:2]
            aspect = w / h
            img_h = 10 / 2.54  # 10cm in inches
            img_w = img_h * aspect
            ax.imshow(img)
            ax.axis('off')
            fig.subplots_adjust(0, 0, 1, 1)
            pdf.savefig(fig)
            plt.close(fig)      
            bar.update(5)           

            fig, ax = plt.subplots(figsize=(8.27,11.69))
            img =  mpimg.imread(savepath +star_name+'_Si_abundances_MAIN.png')
            h, w = img.shape[:2]
            aspect = w / h
            img_h = 10 / 2.54  # 10cm in inches
            img_w = img_h * aspect
            ax.imshow(img)
            ax.axis('off')
            fig.subplots_adjust(0, 0, 1, 1)
            pdf.savefig(fig)
            plt.close(fig)      
            bar.update(6)                                        
         
            fig, ax = plt.subplots(figsize=(8.27,11.69))
            img =  mpimg.imread(savepath + star_name+'_Si_abundances_1_fixed.png')
            h, w = img.shape[:2]
            aspect = w / h
            img_h = 10 / 2.54  # 10cm in inches
            img_w = img_h * aspect
            ax.imshow(img)
            ax.axis('off')
            fig.subplots_adjust(0, 0, 1, 1)
            pdf.savefig(fig)
            plt.close(fig)      
            bar.update(7)   

            fig, ax = plt.subplots(figsize=(8.27,11.69))
            img =  mpimg.imread(savepath + star_name+'_Si_abundances_2_fixed.png')
            h, w = img.shape[:2]
            aspect = w / h
            img_h = 10 / 2.54  # 10cm in inches
            img_w = img_h * aspect
            ax.imshow(img)
            ax.axis('off')
            fig.subplots_adjust(0, 0, 1, 1)
            pdf.savefig(fig)
            plt.close(fig)      
            bar.update(8)   
        
            
        bar.finish()

        files_rm = glob.glob(savepath+"*%s*"%star_name)
        for filei in files_rm:
            if  "summaryplots.pdf" not in filei:
                if "csv" not in filei: 
                    print("Deleting ", filei)
                    os.unlink(filei)

    else:
        print("Skipping %s because plot files don't exist... \n"%(star_name))
    
    return None

def abundances_element(ref_file,el_input,star_db,star_measure,linestab,
                       FW_EW_measure,parameterstable,
                       savepath,mult_use="all",mult_highlight=None,interpolate=False,
                       microreduced="montecarlo", unc = "y", auto="y"):
    
    star_name=findtable(star_db)['ID'][findtable(star_db)['Ref_file']==ref_file][0]    
    if "Class." in findtable(star_db).colnames:
        spc_star=findtable(star_db)['Class.'][findtable(star_db)['Ref_file']==ref_file][0]
    elif "SpType" in findtable(star_db).colnames:
        spc_star=findtable(star_db)['SpType'][findtable(star_db)['Ref_file']==ref_file][0]        
    else:
        spc_star = "-"

    inputparameters=star_name+'_input_parameters.ecsv'
    parameterstab2=findtable(inputparameters,format='ecsv',path=savepath,delimiter='\s')
    
    FW_EW_measure_tab = findtable(FW_EW_measure)

    format_i = "csv" if ".csv" in parameterstable else "ecsv" if ".ecsv" in parameterstable else "ascii"
    delimiter_i = "," if ".csv" in parameterstable else "\s" if ".ecsv" in parameterstable else "\s"    
    final_parameters=findtable(parameterstable,path=savepath,format=format_i,delimiter=delimiter_i)
    parameterstab = final_parameters[final_parameters["Ref_file"]==ref_file]

    if el_input=='O':
        abun_centr='8.75'
    elif el_input=='N':
        abun_centr='7.80'
    elif el_input=='Mg':
        abun_centr='7.55'
    elif el_input=='H_He':
        abun_centr='7.80'
    elif el_input=='C':
        abun_centr='8.35'
    elif el_input=='Si':
        abun_centr='7.50'

    if parameterstab2 is None:
        
        print(ref_file,parameterstab["T"][0],parameterstab["Terr"][0],parameterstab["T+"][0], parameterstab["T-"][0],
                                parameterstab["g"][0],parameterstab["g+"][0], parameterstab["g-"][0], parameterstab["He"][0],parameterstab["Q"][0],
                                parameterstab["b"][0], parameterstab["micro"][0], [abun_centr])
        # SI INPUT_PARAMETERS NO EXISTE, CREARLO! 
        parameterstab2 = Table([[ref_file],parameterstab["T"],parameterstab["Terr"],parameterstab["T+"], parameterstab["T-"],
                                parameterstab["g"],parameterstab["g+"], parameterstab["g-"], parameterstab["He"],parameterstab["Q"],
                                parameterstab["b"], parameterstab["micro"], [abun_centr]],
                               names = ["Ref_file", "T", "Terr", "T+", "T-", "g", "g+", "g-", "He", "Q", "b", "micro", "abun"])
        parameterstab2.write(savepath + inputparameters,format='ascii.ecsv',overwrite=True)
    else:
        parameterstab2.add_row(parameterstab2[-1])
        parameterstab2["abun"][-1] = float(abun_centr)
        parameterstab2.write(savepath + inputparameters,format='ascii.ecsv',overwrite=True)

    if microreduced != "fixed":
        if "micro"+el_input in parameterstab.colnames:
            final_micro = parameterstab['micro'+el_input][0]
        else:
            final_micro = parameterstab['micro'][0]
        microlist = sorted(set(FW_EW_measure_tab['micro']))
        final_micro = microlist[min(range(len(microlist)),key=lambda i : abs(int(microlist[i])-float(final_micro)))]
    else: final_micro = parameterstab['micro'][0]

    redo = 'y'
    redo_count=0

    linetab=findtable(linestab,format='csv',path=savepath)

    # DEFINIMOS LAS LÍNEAS CON LAS QUE CALCULAREMOS LA MICRO! 

    if el_input == 'Si':
        linesvec = ['SiIII4552','SiIII4567','SiIII4574']
    elif el_input == "O":
        linesvec = linetab["label"][linetab["flag"] == "green"]
        linesvec = [linei for linei in linesvec if "II" in linei and "III" not in linei]
    elif el_input == "C":
        linesvec = linetab["label"][linetab["flag"] == "green"]
        linesvec = [linei for linei in linesvec if "II" in linei and "III" not in linei]
    elif el_input == "N":
        linesvec = linetab["label"][linetab["flag"] == "green"]
        # linesvec = [linei for linei in linesvec if "II" in linei and "III" not in linei]                
    else:
        linesvec = linetab["label"][linetab["flag"] == "green"]

    while redo == 'y':

        fixed_micro = final_micro

        abundances(ref_file, spc_star, el_input, star_db, star_measure, FW_EW_measure,
                   inputparameters,parameterstable,linestab,savepath,linesvec = linesvec,  
                   mult_use=mult_use, mult_highlight=mult_highlight, interpolate= interpolate, 
                   ewcut=5,redo_count=redo_count,microreduced=microreduced, auto = auto, unc= unc)
        
        final_parameters=findtable(parameterstable,path=savepath,format=format_i,delimiter=delimiter_i)

        if microreduced != "fixed":

            if "micro"+el_input in parameterstab.colnames:
                final_micro = final_parameters['micro'+el_input][0]
            else:
                final_micro = final_parameters['micro'][0]

            final_micro = microlist[min(range(len(microlist)),key=lambda i : abs(int(microlist[i])-float(final_micro)))]

            redo = 'n'

        else: 
            redo = "n"

    print('Done')
    
    return None
