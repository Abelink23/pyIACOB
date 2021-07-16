'''=============================================================================
Script to make different plots related to MAUI emulated spectra
============================================================================='''
import sys; sys.path.append('../')

from spec import *

'''======================== Find and load the tables ========================'''
table = findtable('IACOB_O9BAs_SNR20.fits') # file where quality flags are
table_REF = findtable('RVEWFWs_O9BAs.fits') # file where RVs, EWs and FWs are
table_IB = findtable('IB_results_ver5.txt') # file where vsini and vmac are
table_IB.remove_columns(['filename','line','snr'])
results = findtable('MAUI_results.fits')

table_f = join(table,table_REF,keys='ID')
table_f = join(table_f,table_IB,keys='ID')
table_f = join(table_f,results,keys='ID')
table_f = table_f[[str(i['Teff'])!='nan' and str(i['lgf'])!='nan' for i in table_f]]

emul = findtable('RVEWFW_Hb_func=vrgH_vsini=10-410_emulated.fits')
try: emul.rename_column('vsini','vsini_py')
except: pass
emul['Teff'] = [int(re.findall('[0-9]+',i)[0]) for i in emul['Name']]
emul['lgf'] = [int(re.findall('[0-9]+',i)[1])/100 for i in emul['Name']]
emul['vsini'] = [int(re.findall('[0-9]+',i)[2]) for i in emul['Name']]
emul['R'] = [int(re.findall('[0-9]+',i)[3]) for i in emul['Name']]

vsinis = list(set(emul['vsini']))


'''========================= Emul 34FW-14FW Hb vs lgf ======================='''
fig, ax = plt.subplots(figsize=(7,6),tight_layout=True)

#table_i = emul[emul['R'] == 5000]
table_i = emul[(emul['R'] == 5000) & (emul['Teff'] == 22000)]

im = ax.scatter(table_i['FW34Hb']-table_i['FW14Hb'],5.39-table_i['lgf'],
                c=table_i['vsini'],cmap='cool',alpha=0.7)

cax = fig.add_axes([0.2,0.15,0.02,0.4])

for v in vsinis:
    emul_v = table_i[table_i['vsini'] == v]
    ax.plot(emul_v['FW34Hb']-emul_v['FW14Hb'],5.39-emul_v['lgf'],
            c='gray',ls='--',alpha=0.7,zorder=0)

cbar = fig.colorbar(im,cax=cax,label='vsini')

ax.set_xlim(0.76,14.49); ax.set_ylim(2.26,4.49)
ax.tick_params(direction='in',top='on')
ax.set_xlabel(r"3/4-HM H$_{\beta}$ - 1/4-HM H$_{\beta}$",size=13)
ax.set_ylabel(r"log($\mathcal{L}$/$\mathcal{L}_{\odot}$)",size=13)

ax.figure.savefig(maindir+'tmp_plots/Emul_Fig1.jpg',format='jpg',dpi=300)



'''========================= Abel 34FW-14FW Hb vs lgf ======================='''
fig, ax = plt.subplots(figsize=(7,6),tight_layout=True)

table_i = table_f[(table_f['QIB'] > 3)]

im = ax.scatter(table_i['FW34Hb']-table_i['FW14Hb'],5.39-table_i['lgf'],
                c=table_i['vsini'],cmap='cool',alpha=0.7)

cax = fig.add_axes([0.2,0.15,0.02,0.4])

cbar = fig.colorbar(im,cax=cax,label='vsini')

ax.set_xlim(0.76,14.49); ax.set_ylim(2.26,4.49)
ax.tick_params(direction='in',top='on')
ax.set_xlabel(r"3/4$\,$HM Hb - 1/4$\,$HM Hb",size=13)
ax.set_ylabel(r"log($\mathcal{L}$/$\mathcal{L}_{\odot}$)",size=13)

ax.figure.savefig(maindir+'tmp_plots/Emul_Fig2.jpg',format='jpg',dpi=300)


'''========================== vsini IB vs vsini_py =========================='''
fig, ax = plt.subplots(figsize=(7,6),tight_layout=True)

#table_i = emul
#table_i = emul[(emul['vsini_py'] < 200) & (emul['vsini'] > 200)]
table_i = emul[emul['R'] >= 5000]
#table_i = emul[(emul['vsini_py'] < 140) & (emul['vsini'] > 150) & (emul['R'] > 6000)]


im = ax.scatter(table_i['vsini'],table_i['vsini_py'],
                c=table_i['sig1'],cmap='gnuplot_r',alpha=0.7)

ax.plot([0,410],[0,410],'k')
cax = fig.add_axes([0.8,0.15,0.02,0.4])

cbar = fig.colorbar(im,cax=cax,label='sigma 1')

ax.tick_params(direction='in',top='on')
ax.set_xlabel("vsini [km/s] - IACOB-broad",size=13)
ax.set_ylabel("vsini [km/s] - pyIACOB",size=13)

ax.figure.savefig(maindir+'tmp_plots/Emul_Fig3.jpg',format='jpg',dpi=300)



'''=============================== EW Hb vs lgf ============================='''
fig, ax = plt.subplots(figsize=(7,6),tight_layout=True)

table_i = table_f

color = table_i['LC_code']
im = ax.scatter(table_i['EWHb'],5.39-table_i['lgf'],
     c=color,cmap='gnuplot_r',alpha=0.7)

cax = fig.add_axes([0.2,0.15,0.02,0.4])

cbar = fig.colorbar(im,cax=cax,label='LC code')
cbar.set_ticks([1,2,3,4,5]); cbar.set_ticklabels(['I','II','III','IV','V'])

ax.tick_params(direction='in',top='on')
ax.set_xlabel(r"EW Hb",size=13)
ax.set_ylabel(r"log($\mathcal{L}$/$\mathcal{L}_{\odot}$)",size=13)

ax.figure.savefig(maindir+'tmp_plots/Emul_Fig4.jpg',format='jpg',dpi=300)
