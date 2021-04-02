'''=============================================================================
Script to make different plots for paper PerOB1: Kinematics
============================================================================='''
import sys
sys.path.append('../')

import matplotlib.path as mpath
from matplotlib.ticker import AutoMinorLocator
from astropy.table import join, setdiff
from scipy.stats import gaussian_kde

from spec import *
from RV import *
from mist import *

from tools_plt import *


'''======================== Find and load the tables ========================'''
table = findtable('IACOB_O9BAs_SNR20.fits') # file where quality flags are
table_REF = findtable('O9BAs_RVEWFWs.fits') # file where RVs, EWs and FWs are
table_REF.remove_columns(['mySpC','SpC','SpT_code','LC_code'])
table_IB = findtable('IB_results_ver5.txt') # file where vsini and vmac are
table_IB.remove_columns(['filename','line','snr'])
results = findtable('zsummary_results.txt') # file with output from MAUI
results.remove_columns(['vsini','vmac','SpC','FW3414'])
gonzalo_raw = findtable('table_gon_vsini_paper.txt')
gaia = findtable('zsummary_results_Gaia.txt') # file Gaia for stars from MAUI

table['Name'] = [i.strip() for i in table['Name']]
table_REF['Name'] = [i.strip() for i in table_REF['Name']]

table_f = join(table,table_REF,keys='Name')
table_f = join(table_f,table_IB,keys='Name')
table_f = join(table_f,results,keys='Name')
table_gaia = join(table_f,gaia,keys='Name')
gonzalo = setdiff(gonzalo_raw,table_f,keys='Name')

#table_f['FIES'].sum()+table_f['HERMES'].sum()+table_f['FEROS'].sum()
#table['FIES'].sum()+table['HERMES'].sum()+table['FEROS'].sum()

'''=============================== Grids MAUI ==============================='''
names = ['Grids coverage','BSgs_CNOSiMg_old','BDws_CNOSIMg_old','OBSgs_hot_NOSi_new','BSgs_cool_NOSi_new']
box_all = [[4.543,4.290,4.290,4.146,4.146,4.543,4.543],
[2.391,2.391,3.092,3.092,4.391,4.391,2.391]]
# nlte_10.1.6_SOLAR_expoclump_2019-10-24.idl
box1 = [[4.190,4.477,4.477,4.190,4.190],[3.785,3.785,4.391,4.391,3.785]]
# nlte_10.1.6_bdwarfs_SOLAR_2020-01-29.idl
box2 = [[4.290,4.543,4.543,4.290,4.290],[2.391,2.391,3.889,3.889,2.391]]
# nlte_10.4.7_OB.Sg_SOLAR_2021-01-23.idl
box3 = [[4.399,4.544,4.544,4.399,4.399],[3.488,3.488,4.386,4.386,3.488]]
# nlte_10.4.7_late.bsgs_SOLAR_expoclump_NOSi.djl_2021-02-06.idl
box4 = [[4.146,4.322,4.322,4.146,4.146],[3.092,3.092,4.391,4.391,3.092]]
grids = [box_all,box1,box2,box3,box4]

'''=============================== Parameters ==============================='''
units = {'pc' : u.pc, 'kpc' : u.kpc, 'au' : u.au, 'lyr' : u.lyr, 'deg' : u.deg}

#mass_list = [.8,.9,1.0,1.1,1.2,1.3,1.5,1.7,2,2.5,3,4,5,7,9,12,15,20,25,32,40,60,85,120]
#mass_list = [.8,.9,1.,1.1,1.2,1.3,1.5,1.7,2,2.5,3,4,5,7,9,12,15,20,25,32,40,60,85]
mass_list = [7,9,12,15,20,25,32,40,60,85]



'''=========================================================================='''
'''================================= Fig. 1 ================================='''
fig_1, ax1 = plt.subplots(figsize=(8,7),tight_layout=True) # 8,7 ; 6,8
ax1_top = ax1.twiny(); ax1_right = ax1.twinx()
xlim = [4.79,4.01]; ylim = [2.2,4.49]

# this, gon, subgroup, grids, masses, density
show = 'this'

# Plot the tracks first:
for i in mass_list:
    mist = trackmist(mass=i,vr=0.0); mist = mist[mist['phase']<=4]
    log_Teff = mist['log_Teff']; log_LLsol = 4*mist['log_Teff']-mist['log_g']-10.61
    s = .5; c = ['gray','k']
    if show == 'density':
        c = ['gray','white']
        if i == 15: s = 1; c = ['white','white']
    ax1.scatter(log_Teff,log_LLsol,s=s,c=c[0],alpha=.5)
    #ax1.scatter(log_Teff,log_LLsol,s=.3,c=mist['surface_he4'],cmap='gnuplot')
    ax1.text(log_Teff[0]+.03,log_LLsol[0]-.04,str(i),c=c[1],fontsize=9)

# Plot the TAMS/ZAMS:
zams = np.arange(2.62,4.19,0.1); tams = np.arange(2.84,3.88,0.1)
polyzams = np.poly1d([-.080,.743,2.973]); polytams = np.poly1d([-.144,.982,2.618])
ax1.plot(polyzams(zams),zams,'--',c='white',alpha=.5)
ax1.plot(polytams(tams),tams,'--',c='white',alpha=.5)

# Plot the results from MAUI creating the sub-tables:
log_Teff = np.asarray(4+np.log10(table_f['Teff']))
log_LLsol = np.asarray(5.39-table_f['lgf'])
points = np.column_stack([log_Teff,log_LLsol])

if show == 'this':
    ax1.scatter(log_Teff,log_LLsol,s=6,c='b',label='This work')
    ax1.scatter(3+np.log10(gonzalo['Teff']),gonzalo['logL'],s=6,c='purple',label='Holgado,G. thesis 2019')

if show == 'gon':
    ax1.scatter(3+np.log10(gonzalo_raw['Teff']),gonzalo_raw['logL'],s=6,c='purple',label='Holgado,G. thesis 2019')

if show == 'subgroup':
    table_i = table_f[(table_f['LC_code']>0) & (table_f['LC_code']<3)]
    #table_i = table_f[(table_f['LC_code']>0) & (table_f['LC_code']<3) & \
    #                  (table_f['SpT_code']>=1.9) & (table_f['SpT_code']<=2.7)]
    #table_i = table_f[(table_f['SpT_code']>=1.9) & (table_f['SpT_code']<=2.15)]

    color = table_i['SpT_code']
    log_Teff = np.asarray(4+np.log10(table_i['Teff']))
    log_LLsol = np.asarray(5.39-table_i['lgf'])

    im = ax1.scatter(log_Teff,log_LLsol,c=color,s=20,lw=.2,ec='k',cmap='gnuplot')
    cax = fig_1.add_axes([0.8,0.33,0.02,0.55])
    cbar = fig_1.colorbar(im,cax=cax,label='SpT')
    cbar.set_ticks([1.9,2.1,2.3,2.5,2.7]); cbar.set_ticklabels(['O9','B1','B3','B5','B7'])


if show == 'grids':

    ax1.scatter(3+np.log10(gonzalo['Teff']),gonzalo['logL'],s=6,c='purple',label='Holgado,G. thesis 2019')

    # Change names[1:],grids[1:] for individual grids
    for name,grid in zip(names[:1],grids[:1]):
        ax1.plot(grid[0],grid[1],lw=.8,ls='--')
        verts = np.array([grid[0],grid[1]]).T
        path = mpath.Path(verts)
        inout = path.contains_points(points)
        log_Teff_in,log_LLsol_in = points[inout].T

        ax1.scatter(log_Teff_in,log_LLsol_in,s=6,label=name)

elif show == 'masses':

    ax1.scatter(3+np.log10(gonzalo['Teff']),gonzalo['logL'],s=6,c='purple',label='Holgado,G. thesis 2019')
    ax1.scatter(log_Teff,log_LLsol,s=6,c='b',label='This work (M<15Msol)')

    # Create a table and plot only the stars above a 15Msol
    mist = trackmist(mass=15,vr=0.0)
    log_Teff = [5]+mist['log_Teff'].tolist()+[4,5,5]
    log_LLsol = [3]+(4*mist['log_Teff']-mist['log_g']-10.61).tolist()+[4.5,4.5,3]
    verts = np.array([log_Teff,log_LLsol]).T
    path = mpath.Path(verts); inout = path.contains_points(points)
    log_Teff_in,log_LLsol_in = points[inout].T
    ax1.scatter(log_Teff_in,log_LLsol_in,s=6,c='limegreen',label='This work (M>15Msol)')
    #ax1.plot(verts.T[0],verts.T[1],c='k')
    results_in = table_f[path.contains_points(points)]['Name']
    table_15 = table_f[[i['Name'].strip() in results_in for i in table_f]]

    # Create a table with only the stars above a centain mass
    mist = trackmist(mass=32,vr=0.0)
    log_Teff = [5]+mist['log_Teff'].tolist()+[4,5,5]
    log_LLsol = [3.6]+(4*mist['log_Teff']-mist['log_g']-10.61).tolist()+[4.5,4.5,3.6]
    verts = np.array([log_Teff,log_LLsol]).T
    path = mpath.Path(verts); inout = path.contains_points(points)
    log_Teff_in,log_LLsol_in = points[inout].T
    ax1.scatter(log_Teff_in,log_LLsol_in,s=6,c='orange',label='This work (M>30Msol)')
    #ax1.plot(verts.T[0],verts.T[1],c='k')
    results_in = table_f[path.contains_points(points)]['Name']
    table_30 = table_f[[i['Name'].strip() in results_in for i in table_f]]

elif show == 'density':

    # Right limit of MAUI grid
    ax1.plot([4.29,4.29,4.146,4.146],[2.391,3.092,3.092,4.391],c='white',lw=.5)

    # Plot Gonzalo's sample
    ax1.scatter(3+np.log10(gonzalo['Teff']),gonzalo['logL'],s=6,
                c='white',label='Holgado,G. thesis 2019')

    # Plot all stars from MAUI
    ax1.scatter(log_Teff,log_LLsol,s=6,c='white',label='This work')

    nbins = 50
    x = np.concatenate((log_Teff,np.asarray(3+np.log10(gonzalo['Teff']))))
    y = np.concatenate((log_LLsol,gonzalo['logL']))
    xy = np.vstack([x,y]); k = gaussian_kde(xy)
    xi, yi = np.mgrid[xlim[0]:xlim[1]:nbins*1j,ylim[0]:ylim[1]:nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    ax1.pcolormesh(xi,yi,zi.reshape(xi.shape),
        shading='gouraud',cmap='gnuplot',zorder=0,vmin=0.05,vmax=4.5)

# Fast rot:
#table_fast = table_f[table_f['vsini'] > 100]
#log_Teff = np.asarray(4+np.log10(table_fast['Teff']))
#log_LLsol = np.asarray(5.39-table_fast['lgf'])
#ax1.scatter(log_Teff,log_LLsol,c='k',marker='+',s=30,alpha=.5,label='Fast rotators')

if show == 'density': c = 'white'
else: c = 'k'
ax1.tick_params(direction='in',top='on',right='on',which='both',length=4,color=c)
#ax1.set_yticks(np.arange(2.3, 4.4, .5))
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator())
ax1.invert_xaxis()
ax1.set_xlabel(r"log(T$_{eff})\,$[K]",size=13)
ax1.set_ylabel(r"log($\mathcal{L}$/$\mathcal{L}_{\odot}$)",size=13)
ax1.set_xlim(xlim); ax1.set_ylim(ylim)
ax1.legend(ncol=1,loc=3,fontsize=9,borderaxespad=1.5)

ax1_top.tick_params(direction='in',top='off')
ax1_top.set_xlim(ax1.get_xlim())
ax1_top.set_xticks(ax1.get_xticks()[1:-1])
ax1_top.set_xticklabels([round(10**i/1000)*1000 for i in ax1.get_xticks()[1:-1]])
ax1_top.set_xlabel(r"T$_{eff}\,$[K]",size=13,labelpad=13)

ax1_right.tick_params(direction='in',right='off')
ax1_right.set_ylim(ax1.get_ylim())
ax1_right.set_yticks(ax1.get_yticks()[1:-1])
ax1_right.set_yticklabels([round(5.39-i,1) for i in ax1.get_yticks()[1:-1]])
ax1_right.set_ylabel(r"log(g$_{F}$)",size=13,rotation=-90,labelpad=15)

ax1.figure.savefig(maindir+'tmp_plots/Bs_Fig1_%s.png' % show,format='png',dpi=300)
plt.show(block=False)


'''=========================================================================='''
'''================================= Fig. 2 ================================='''
fig_2, ax2 = plt.subplots(figsize=(7,6),tight_layout=True) # 7,6 9,3

ax2.scatter(3+np.log10(gonzalo['Teff']),gonzalo['vsini(GOF)'],s=25,ec='k',fc='none',marker='+')

table_i = table_f
all_Teff = (4+np.log10(table_i['Teff'])).tolist()+(3+np.log10(gonzalo['Teff'])).tolist()
all_vsin = table_i['vsini'].tolist()+gonzalo['vsini(GOF)'].tolist()

color = (5.39-table_i['lgf']).tolist()+gonzalo['logL'].tolist()
im = ax2.scatter(all_Teff,all_vsin,c=color,s=25,lw=.2,ec='k',label='M > 15Msol',cmap='gnuplot')

cax = fig_2.add_axes([0.86,0.33,0.02,0.55])

fig_2.colorbar(im,cax=cax,label=r"log($\mathcal{L}$/$\mathcal{L}_{\odot}$)",aspect=50)

# Enable next 4 for only the colors without bar
#table_15d = setdiff(table_15,table_30,keys='Name')
#ax2.scatter(3+np.log10(gonzalo['Teff']),gonzalo['vsini(GOF)'],c='purple',lw=.2,ec='k',s=25,label='Holgado,G. thesis 2019')
#ax2.scatter(4+np.log10(table_15d['Teff']),table_15d['vsini'],c='limegreen',lw=.2,ec='k',s=25,label='M > 15Msol')
#ax2.scatter(4+np.log10(table_30['Teff']),table_30['vsini'],c='orange',lw=.2,ec='k',s=25,label='M > 30Msol')
#ax2.legend(ncol=1,loc=1,fontsize=9,borderaxespad=1.5)

ax2.tick_params(direction='in',top='on',right='on',which='both')
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())
ax2.invert_xaxis()
ax2.set_xlim(4.72,4.14)
ax2.set_xlabel(r"log(T$_{eff})\,$[K]",size=13)
ax2.set_ylabel(r"vsin($i$) [km/s]",size=13)
ax2.set_ylim([-20,480])
ax2.figure.savefig(maindir+'tmp_plots/Bs_Fig2.png',format='png',dpi=300)


'''=========================== 34FW-14FW Hb vs lgf =========================='''
'''================================= Fig. 3a ================================'''
fig_3a, ax3a = plt.subplots(figsize=(7,6),tight_layout=True)

table_i = table_f

color = table_i['LC_code']
im = ax3a.scatter(table_i['FW34Hb']-table_i['FW14Hb'],5.39-table_i['lgf'],
     c=color,cmap='gnuplot_r',alpha=0.7)

cax = fig_3a.add_axes([0.2,0.15,0.02,0.4])

cbar = fig_3a.colorbar(im,cax=cax,label='LC code')
cbar.set_ticks([1,2,3,4,5]); cbar.set_ticklabels(['I','II','III','IV','V'])

ax3a.tick_params(direction='in',top='on')
ax3a.set_xlabel(r"3/4$\,$HM Hb - 1/4$\,$HM Hb",size=13)
ax3a.set_ylabel(r"log($\mathcal{L}$/$\mathcal{L}_{\odot}$)",size=13)

ax3a.figure.savefig(maindir+'tmp_plots/Bs_Fig3a.png',format='png',dpi=300)


'''=============================== EW Hb vs lgf ============================='''
'''================================= Fig. 3b ================================'''
fig_3b, ax3b = plt.subplots(figsize=(7,6),tight_layout=True)

table_i = table_f

color = table_i['LC_code']
im = ax3b.scatter(table_i['EWHb'],5.39-table_i['lgf'],
     c=color,cmap='gnuplot_r',alpha=0.7)

cax = fig_3b.add_axes([0.2,0.15,0.02,0.4])

cbar = fig_3b.colorbar(im,cax=cax,label='LC code')
cbar.set_ticks([1,2,3,4,5]); cbar.set_ticklabels(['I','II','III','IV','V'])

ax3b.tick_params(direction='in',top='on')
ax3b.set_xlabel(r"EW Hb",size=13)
ax3b.set_ylabel(r"log($\mathcal{L}$/$\mathcal{L}_{\odot}$)",size=13)

ax3b.figure.savefig(maindir+'tmp_plots/Bs_Fig3b.png',format='png',dpi=300)


'''========================== sHDR Microturbulence =========================='''
'''================================= Fig. 4a ================================'''
fig_4, ax4 = plt.subplots(figsize=(6,6),tight_layout=True)

table_i = table_f

# Plot the tracks first:
for i in mass_list:
    mist = trackmist(mass=i,vr=0.0)
    mist = mist[mist['phase']<=4]
    log_Teff = mist['log_Teff'][(mist['log_Teff']>4.15) & (mist['log_Teff']<4.55)]
    log_LLsol = (4*mist['log_Teff']-mist['log_g']-10.61)[(mist['log_Teff']>4.15) & (mist['log_Teff']<4.55)]
    ax4.scatter(log_Teff,log_LLsol,s=.3,c='gray')
    ax4.text(log_Teff[0]+.01,log_LLsol[0]-.04,str(i),fontsize=7,clip_on=True)

log_Teff = np.asarray(4+np.log10(table_i['Teff']))
log_LLsol = np.asarray(5.39-table_i['lgf'])

color = table_i['Micro']
im = ax4.scatter(log_Teff,log_LLsol,c=color,s=100,cmap='gnuplot',alpha=.7)
im.set_clim(0,24)
cax = fig_4.add_axes([0.8,0.12,0.03,0.4])

fig_4.colorbar(im,cax=cax,label='Microturbulence [km/s]')

ax4.tick_params(direction='in',top='on')
ax4.invert_xaxis()
ax4.set_xlabel(r"log(T$_{eff})\,$[K]",size=13)
ax4.set_ylabel(r"log($\mathcal{L}$/$\mathcal{L}_{\odot}$)",size=13)
ax4.set_ylim(top=4.49)
ax4.figure.savefig(maindir+'tmp_plots/Bs_Fig4a.png',format='png',dpi=300)


'''========================== sHDR Macroturbulence =========================='''
'''================================= Fig. 4b ================================'''
fig_4, ax4 = plt.subplots(figsize=(6,6),tight_layout=True)

table_i = table_f

# Plot the tracks first:
for i in mass_list:
    mist = trackmist(mass=i,vr=0.0)
    mist = mist[mist['phase']<=4]
    log_Teff = mist['log_Teff'][(mist['log_Teff']>4.15) & (mist['log_Teff']<4.55)]
    log_LLsol = (4*mist['log_Teff']-mist['log_g']-10.61)[(mist['log_Teff']>4.15) & (mist['log_Teff']<4.55)]
    ax4.scatter(log_Teff,log_LLsol,s=.3,c='gray')
    ax4.text(log_Teff[0]+.01,log_LLsol[0]-.04,str(i),fontsize=7,clip_on=True)

log_Teff = np.asarray(4+np.log10(table_i['Teff']))
log_LLsol = np.asarray(5.39-table_i['lgf'])

color = table_i['rt_mac']
im = ax4.scatter(log_Teff,log_LLsol,c=color,s=100,cmap='gnuplot',alpha=0.7)
cax = fig_4.add_axes([0.8,0.12,0.03,0.4])

fig_4.colorbar(im,cax=cax,label='Macroturbulence [km/s]')
im.set_clim(0,99)

ax4.tick_params(direction='in',top='on')
ax4.invert_xaxis()
ax4.set_xlabel(r"log(T$_{eff})\,$[K]",size=13)
ax4.set_ylabel(r"log($\mathcal{L}$/$\mathcal{L}_{\odot}$)",size=13)
ax4.set_ylim(top=4.49)
ax4.figure.savefig(maindir+'tmp_plots/Bs_Fig4b.png',format='png',dpi=300)


'''=============================== sHDR wind Q =============================='''
'''================================= Fig. 4c ================================'''
fig_4, ax4 = plt.subplots(figsize=(6,6),tight_layout=True)

table_i = table_f

# Plot the tracks first:
for i in mass_list:
    mist = trackmist(mass=i,vr=0.0)
    mist = mist[mist['phase']<=4]
    log_Teff = mist['log_Teff'][(mist['log_Teff']>4.15) & (mist['log_Teff']<4.55)]
    log_LLsol = (4*mist['log_Teff']-mist['log_g']-10.61)[(mist['log_Teff']>4.15) & (mist['log_Teff']<4.55)]
    ax4.scatter(log_Teff,log_LLsol,s=.3,c='gray')
    ax4.text(log_Teff[0]+.01,log_LLsol[0]-.04,str(i),fontsize=7,clip_on=True)

log_Teff = np.asarray(4+np.log10(table_i['Teff']))
log_LLsol = np.asarray(5.39-table_i['lgf'])

color = table_i['logQs']-10
im = ax4.scatter(log_Teff,log_LLsol,c=color,s=100,cmap='gnuplot',alpha=0.7)
cax = fig_4.add_axes([0.8,0.12,0.03,0.4])

fig_4.colorbar(im,cax=cax,label=r'$log(Q)$')

ax4.tick_params(direction='in',top='on')
ax4.invert_xaxis()
ax4.set_xlabel(r"log(T$_{eff})\,$[K]",size=13)
ax4.set_ylabel(r"log($\mathcal{L}$/$\mathcal{L}_{\odot}$)",size=13)
ax4.set_ylim(top=4.49)
ax4.figure.savefig(maindir+'tmp_plots/Bs_Fig4c.png',format='png',dpi=300)


'''================================= Helium ================================='''
'''================================= Fig. 4d ================================'''
fig_4, ax4 = plt.subplots(figsize=(6,6),tight_layout=True)

table_i = table_f

# Plot the tracks first:
for i in mass_list:
    mist = trackmist(mass=i,vr=0.0)
    mist = mist[mist['phase']<=4]
    log_Teff = mist['log_Teff'][(mist['log_Teff']>4.15) & (mist['log_Teff']<4.55)]
    log_LLsol = (4*mist['log_Teff']-mist['log_g']-10.61)[(mist['log_Teff']>4.15) & (mist['log_Teff']<4.55)]
    ax4.scatter(log_Teff,log_LLsol,s=.3,c='gray')
    ax4.text(log_Teff[0]+.01,log_LLsol[0]-.04,str(i),fontsize=7,clip_on=True)

log_Teff = np.asarray(4+np.log10(table_i['Teff']))
log_LLsol = np.asarray(5.39-table_i['lgf'])

color = table_i['He']
im = ax4.scatter(log_Teff,log_LLsol,c=color,s=100,cmap='gnuplot',alpha=0.7)
cax = fig_4.add_axes([0.8,0.12,0.03,0.4])

fig_4.colorbar(im,cax=cax,label='He')

ax4.tick_params(direction='in',top='on')
ax4.invert_xaxis()
ax4.set_xlabel(r"log(T$_{eff})\,$[K]",size=13)
ax4.set_ylabel(r"log($\mathcal{L}$/$\mathcal{L}_{\odot}$)",size=13)
ax4.set_ylim(top=4.49)
ax4.figure.savefig(maindir+'tmp_plots/Bs_Fig4d.png',format='png',dpi=300)


'''============================== Micro & Macro ============================='''
'''================================= Fig. 5 ================================='''
fig_5, ax5 = plt.subplots(figsize=(6,6),tight_layout=True)

table_i = table_f

color = table_i['LC_code']
im = ax5.scatter(table_i['rt_mac'],table_i['Micro'],c=color,s=30,lw=.3,ec='k',cmap='gnuplot_r',alpha=0.7)

cax = fig_5.add_axes([0.8,0.12,0.03,0.4])

cbar = fig_5.colorbar(im,cax=cax,label='LC code')
cbar.set_ticks([1,2,3,4,5]); cbar.set_ticklabels(['I','II','III','IV','V'])

ax5.tick_params(direction='in',top='on')
ax5.set_xlabel('Macroturbulence [km/s]',size=13)
ax5.set_ylabel('Microturbulence [km/s]',size=13)
ax5.set_xlim(right=150)
ax5.figure.savefig(maindir+'tmp_plots/Bs_Fig5.png',format='png',dpi=300)


'''============================= Histogram Teff ============================='''
'''================================= Fig. 6 ================================='''
fig_6, ax6 = plt.subplots(figsize=(9,3),tight_layout=True)

table_i = setdiff(table_f,table_15,keys='Name')

N,bins,groups = ax6.hist([3+np.log10(gonzalo['Teff']),4+np.log10(table_15['Teff']),
4+np.log10(table_i['Teff'])],
ec = 'k',lw=.2,stacked=True,bins=np.arange(4.16,4.70,0.02),
label=['Holgado,G. thesis 2019','M > 15Msol','M < 15Msol'])

for group,color in zip(groups,['mediumvioletred','r','b']):
    for patch in group: patch.set_facecolor(color)

ax6.tick_params(direction='in',top='on',length=10)
ax6.invert_xaxis()
ax6.set_xlabel(r"log(T$_{eff})\,$[K]",size=13)
ax6.set_ylabel('#',size=13)
ax6.set_xlim(ax2.get_xlim())
ax6.legend(loc=2)
ax6.figure.savefig(maindir+'tmp_plots/Bs_Fig6.png',format='png',dpi=300)


'''========================= Teff mine vs Teff Gonzalo ======================'''
'''================================= Fig. 7 ================================='''
fig_7, ax7 = plt.subplots(figsize=(7,6),tight_layout=True)

table_i = join(table_f,gonzalo_raw,keys='Name')

ax7.scatter(table_i['Teff_1']*1e4,table_i['Teff_2']*1e3,c='b')

plt.plot([26000,36000],[26000,36000],'k')
plt.plot([26000,36000],[25000,35000],'--k')
plt.plot([26000,36000],[27000,37000],'--k')

ax7.tick_params(direction='in',top='on')
ax7.set_xlim(ax7.get_xlim());ax7.set_ylim(ax7.get_xlim())
#ax7.axis('equal')
ax7.set_xlabel(r"log(T$_{eff})\,$[K] (this work)",size=13)
ax7.set_ylabel(r"log(T$_{eff})\,$[K] (Holgado,G. thesis 2019)",size=13)
ax7.figure.savefig(maindir+'tmp_plots/Bs_Fig7.png',format='png',dpi=300)


'''========================= logg mine vs logg Gonzalo ======================'''
'''================================= Fig. 8 ================================='''
fig_8, ax8 = plt.subplots(figsize=(7,6),tight_layout=True)

table_i = join(table_f,gonzalo_raw,keys='Name')

#ax8.scatter(table_i['lgf'],table_i['logg'],c='b')
ax8.scatter(table_i['lgf'],16+table_i['logg']-4*np.log10(table_i['Teff_2']*1e3),c='b')

#plt.plot([29000,35000],[29000,35000],'k')
#plt.plot([29000,35000],[28000,34000],'--k')
#plt.plot([29000,35000],[30000,36000],'--k')

ax8.tick_params(direction='in',top='on')
ax8.set_xlim(1.01,2.3);ax8.set_ylim(1.01,2.3)
ax8.set_xlabel(r"$log(g)\,[m/s^{2}]$ (this work)",size=13)
ax8.set_ylabel(r"$log(g)\,[m/s^{2}]$ (Holgado,G. thesis 2019)",size=13)
ax8.figure.savefig(maindir+'tmp_plots/Bs_Fig8.png',format='png',dpi=300)


'''========================= Gaia EDR3 good parallaxes ======================'''
'''================================= Fig. 9 ================================='''
fig_9, ax9 = plt.subplots(figsize=(6,8),tight_layout=True)

# Plot the tracks first:
for i in mass_list:
    mist = trackmist(mass=i,vr=0.0)
    mist = mist[mist['phase']<=4]
    log_Teff = mist['log_Teff']
    log_LLsol = 4*mist['log_Teff']-mist['log_g']-10.61
    ax9.scatter(log_Teff,log_LLsol,s=.3,c='gray',alpha=.5)
    ax9.text(log_Teff[0]+.02,log_LLsol[0]-.07,str(i),fontsize=7)

# Plot the results from MAUI creating the sub-tables:
log_Teff = np.asarray(4+np.log10(table_gaia['Teff']))
log_LLsol = np.asarray(5.39-table_gaia['lgf'])

color = table_gaia['parallax_over_error']
im = ax9.scatter(log_Teff,log_LLsol,c=color,s=10,label='M > 15Msol',cmap='jet')
im.set_clim(0,80)

cax = fig_9.add_axes([0.82,0.09,0.04,0.40])

fig_9.colorbar(im,cax=cax,label="Parallax over error",aspect=50)

ax9.tick_params(direction='in',top='on')
ax9.invert_xaxis()
ax9.set_xlabel(r"log(T$_{eff})\,$[K]",size=13)
ax9.set_ylabel(r"log($\mathcal{L}$/$\mathcal{L}_{\odot}$)",size=13)
ax9.set_xlim(4.8,4)
ax9.set_ylim(2.2,4.7)
ax9.figure.savefig(maindir+'tmp_plots/Bs_Fig9.png',format='png',dpi=300)


'''=========================== Histogram magnitudes ========================='''
'''================================= Fig. 10 ================================'''
fig_10, ax10 = plt.subplots(figsize=(6,6),tight_layout=True)

table_i = table_f

ax10.hist(table_i['mag_V'],ec = 'k',bins=np.linspace(1,12,23),color='b',label='This work')

ax10.tick_params(direction='in',top='on')
ax10.invert_xaxis()
ax10.set_xlabel(r"V$_{mag}$",size=13)
ax10.set_ylabel('#',size=13)
ax10.set_xticks(np.linspace(1,12,12))
ax10.set_xlim(1,12)
ax10.legend()
ax10.figure.savefig(maindir+'tmp_plots/Bs_Fig10.png',format='png',dpi=300)


'''========================== Interactive find stars ========================'''
plt.close('all')
fig, ax = plt.subplots(figsize=(6,8),tight_layout=True) # 8,7 ; 6,8
xlim = [4.79,4.01]; ylim = [2.2,4.49]

# Plot the tracks first:
track_0 = []
for i in mass_list:
    mist = trackmist(mass=i,vr=0.0); mist = mist[mist['phase']<=4]
    log_Teff = mist['log_Teff']; log_LLsol = 4*mist['log_Teff']-mist['log_g']-10.61
    s = .5; c = ['gray','k']
    c = ['gray','white']
    if i == 15: s = 1; c = ['white','white']
    track_0.append([log_Teff[0],log_LLsol[0]])
    ax.scatter(log_Teff,log_LLsol,s=s,c=c[0],alpha=.5)
    #ax.scatter(log_Teff,log_LLsol,s=.3,c=mist['surface_he4'],cmap='gnuplot')
    #ax.text(log_Teff[0]+.03,log_LLsol[0]-.04,str(i),c=c[1],fontsize=9)

# Plot the TAMS:
tams = np.arange(2.84,3.88,0.1); polytams = np.poly1d([-.144,.982,2.618])
ax.plot(polytams(tams),tams,'--',c='white',alpha=.5)

# Plot the ZAMS:
ax.plot([i[0] for i in track_0],[i[1] for i in track_0],'-',c='white',alpha=1)

# Plot diagonals with same gravity:
for g in np.arange(4.3,1,-0.3):
    g_x = np.arange(4.05,4.75,0.05)
    g_y = 4*g_x-g-10.61
    ax.plot(g_x,g_y,'--',c='white',alpha=.5,lw=.3)


# Right limit of MAUI grid
ax.plot([4.29,4.29,4.146,4.146],[2.391,3.092,3.092,4.391],c='white',lw=.5)

log_Teff = np.asarray(4+np.log10(table_f['Teff']))
log_LLsol = np.asarray(5.39-table_f['lgf'])

nbins = 50
x = np.concatenate((log_Teff,np.asarray(3+np.log10(gonzalo['Teff']))))
y = np.concatenate((log_LLsol,gonzalo['logL']))
xy = np.vstack([x,y]); k = gaussian_kde(xy)
xi, yi = np.mgrid[xlim[0]:xlim[1]:nbins*1j,ylim[0]:ylim[1]:nbins*1j]
zi = k(np.vstack([xi.flatten(), yi.flatten()]))
ax.pcolormesh(xi,yi,zi.reshape(xi.shape),
    shading='gouraud',cmap='gnuplot',zorder=0,vmin=0.05,vmax=4.5)

pts = ax.scatter(x,y,s=6,c='white',label='This work + Holgado,G.')

ax.tick_params(direction='in',top='on',right='on',which='both',length=4,color='white')
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.invert_xaxis()
ax.set_xlabel(r"log(T$_{eff})\,$[K]",size=13)
ax.set_ylabel(r"log($\mathcal{L}$/$\mathcal{L}_{\odot}$)",size=13)
ax.set_xlim(xlim); ax.set_ylim(ylim)
ax.legend(ncol=1,loc=3,fontsize=9,borderaxespad=1.5)

#af =  AnnoteFinder(x,y,table_f['Name'].tolist()+gonzalo['Name'].tolist(),ax=ax,xtol=.05,ytol=.05)
#fig.canvas.mpl_connect('button_press_event', af)
#plt.show(block=False)

table_i = join(table_f,gonzalo,keys='Name',join_type='outer')

selector = SelectFromCollection(ax, pts)
def accept(event):
    if event.key == "enter":
        selector.disconnect()
        fig.canvas.draw()
        accept.table_selction = table_i[selector.mask]

fig.canvas.mpl_connect("key_press_event", accept)

plt.show(block=False)

hdu = fits.BinTableHDU(data=accept.table_selction.filled(0))
hdu.writeto(maindir+'tables/table_selection.fits',overwrite=True)
