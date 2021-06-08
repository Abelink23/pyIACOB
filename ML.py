from RV import *
import random

from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages(maindir+'tmp_plots/ML_results.pdf')

def gen_ascii_ML(input_table='OBAs_ML_raw.fits',not_do=None,manual=False):

    table = findtable(input_table)
    output = open(maindir+'tmp/results_ML.txt','a')

    for row in table:

        # Skip all sources in the 'not_do' input list :
        if not_do is not None and row['ID'] in not_do: continue
        if row['SNR_best'] < 60: continue
        else: print('Analysing %s' % row['ID'])

        # Determines best list for RV calculation and line for sanity check based on SpT
        if row['SpT_code'] < 2:
            rv_list = 'rv_Os.lst'
            line = 5411.52 # 5592.252

        elif 2 <= row['SpT_code'] < 2.5:
            rv_list = 'rv_Bs.lst'
            line = 4552.622

        elif 2.5 <= row['SpT_code'] < 3:
            rv_list = 'rv_Bs.lst'
            line = 6347.11

        elif row['SpT_code'] >= 3:
            rv_list = 'rv_As.lst'
            line = 4233.129

        # Determines the RV with a default fitting function and width
        fun = 'r'; wid = 15
        best_star = spec(row['ID'],SNR='bestMF')
        best_star.offset = RV0(rv_list,best_star.spectrum,ewcut=50,func=fun,width=wid,tol=150)
        best_star.waveflux(3950,6850) # Applies the offset

        # If the line is >.1A from where should be, you determine new best function, width and line
        RV_A = abs(best_star.fitline(line,func=fun,width=wid,tol=20)['RV_A'])
        if np.isnan(RV_A) or RV_A > 0.1: # 5km/s at 5400A
            print('Fitted line outside offset is %.3f (tolerance is 0.1A)' % RV_A)
            fun,wid,line = func_width(best_star,row['SpT_code'],rv_list,line)

        # For all available spectra withing a limit, create the final output ascii
        goodspec = findstar(row['ID'],SNR=60)
        if len(goodspec) > 10: goodspec = random.sample(goodspec,10)

        for j,n in zip(goodspec,range(len(goodspec))):
            star = spec(j,SNR='best')
            star.offset = RV0(rv_list,star.spectrum,ewcut=50,width=wid,tol=150,func=fun)
            star.waveflux(3950,6850) # Applies the offset

            fig,axs = plt.subplots(3,1,figsize=(14,8),tight_layout=True)
            fig.suptitle(star.filename)
            axs[0].tick_params(direction='in',top='on')
            axs[0].set_xlim(3950,6850)
            axs[0].set_ylim(0.5,1.1)
            axs[1].tick_params(direction='in',top='on')
            axs[1].set_xlim(3950,6850)
            axs[1].set_ylim(0.5,1.1)

            axs[0].plot(star.wave,star.flux,c='orange',lw=.2,alpha=.7,label='original')

            # Remove artifacts in FEROS spectra
            if '_F_' in star.filename:
                i = 0
                std = np.std(star.flux[(star.wave > 4968)&(star.wave < 4985)])
                for win in [
                [4089.9,4090.3],
                [4505.,4508.],
                [4693.0,4695.7],
                [4794.3,4796.5],
                [4900.7,4902.2],
                [6636.,6638.],
                ]:
                    mask = (star.wave > win[0]-star.offset)&(star.wave < win[1]-star.offset)
                    gap = star.flux[mask]
                    if np.std(gap) < 5*std: continue # To skip if the artifact is not there
                    length = len(gap)
                    cont = np.mean(np.concatenate((gap[:5],gap[-5:])))
                    star.flux[mask] = [random.uniform(cont-2*std,cont+2*std) for i in range(length)]
                    axs[0].plot(star.wave[mask],star.flux[mask],c='r',lw=.2,label='FEROS fixed issues'if i==0 else "")
                    i += 1

            star.flux = cleanML(star.wave,star.flux,star.filename,axs,manual=manual)

            star.degrade(resol=5000)

            mask = (star.wave > line-10)&(star.wave < line+10)
            axs[2].plot(star.wave[mask],star.flux[mask],c='b',lw=.2,label='final spectra')
            axs[2].plot([line,line],[min(star.flux[mask]),max(star.flux[mask])],'k',label='ref. line position')
            medval = (max(star.flux[mask]) + min(star.flux[mask]))/2
            medpos = [np.where(star.flux[mask] <= medval)[0][value] for value in (0,-1)]
            center = round((star.wave[mask][medpos[1]]+star.wave[mask][medpos[0]])/2,3)
            axs[2].plot([center,center],[min(star.flux[mask]),max(star.flux[mask])],'g',label='actual position')
            axs[2].set_title('Difference in angstroms is: %.5f' % abs(center - line))
            axs[2].set_ylim()

            star.resamp(10*0.02564975,3950,6850)

            star.export(tail='_RV_ML_%i' % n,extension='.ascii')
            if max(star.flux) > 1.5: output.write(star.filename+' | Max flux: %.3f\n' % max(star.flux))
            if abs(center - line) > 0.2: output.write(star.filename+' | RV difference (A): %.3f\n' % (center - line))

            axs[0].legend(); axs[1].legend(); axs[2].legend()
            pp.savefig(fig)
            plt.close('all')

    output.close()
    pp.close()
    np.savetxt(maindir+'tmp/wavelenghtML.ascii',np.c_[best_star.wave],fmt=('%.4f'))


def func_width(spec,spt_code,rv_list,line):
    next = 'n'
    while next != '':

        line_0 = line
        line = input('Choose line (current is %.3f): ' % line)
        if line.strip() == '': line = line_0
        else: line = float(line)

        spec.plotspec(line-10,line+10)

        fun = '-'
        while fun not in ['g','l','v','r','vr']:
            try:
                fun = input('Choose function to fit between g/l/v/r/vr (default is g): ')
                if fun == '': fun = 'g'
            except: pass
        wid = '-'
        while type(wid) is not float:
            try:
                wid = input('Choose the initial width in angstroms (default is 15): ')
                if wid == '': wid = 15.
                else: wid = float(wid)
            except: pass

        plt.close()

        spec.offset = RV0(rv_list,spec.spectrum,ewcut=50,width=wid,tol=150,func=fun)
        spec.waveflux(3950,6850) # Applies the offset
        spec.cosmic()
        spec.degrade(resol=5000)

        spec.plotspec(line-wid/2,line+wid/2,poslines='OB')

        next = input('Hit return to continue with the rest of spectra or any other key to repeat. ')

        plt.close('all')

    return fun,wid,line


def cosmicML(wave, flux, method='zscore', sigclip=1.5, iter=3, sig_g=None):
    '''
    Parameters
    ----------

    method : str, optional
        Method for the cosmic ray removal strategy. Only zscore (def) or kernel.

    sigclip : float, optional
        Sigma clipping value used to remove rays. Default is 1.5.

    iter : int, optional
        Number of iterations of the sigma clipping to remove cosmic rays.

    sig_g : float, optional
        Sigma of the gaussian function used to construct the kernel.
        Default is the theoretical sigma based on wavelenght and resolution.

    Returns: None (but the flux is replaced and cleaned from rays).
    '''

    if method == 'zscore':
        #www.towardsdatascience.com/removing-spikes-from-raman-spectra-8a9fdda0ac22

        # First we calculated ∇x(i):
        delta_flux = [flux[i+1] - flux[i] for i in np.arange(len(flux)-1)]

        median_int = np.median(delta_flux)
        mad_int = np.median([np.abs(delta_flux - median_int)])
        modified_z_scores = 0.6745 * (delta_flux - median_int) / mad_int
        # The multiplier 0.6745 is the 0.75th quartile of the standard normal
        # distribution, to which the median absolute deviation converges to.

        flux_norm =  np.concatenate(([1],abs(modified_z_scores)))

    elif method == 'kernel':

        resolution = 5000
        dlam = 0.2564975

        if sig_g == None:
            lambda0 = np.mean(wave)
            # Two times the theoretical sigma offers better results
            sig_g = 2*lambda0/(2.35482*float(resolution))
        else: sig_g = float(sig_g)

        x = np.arange(-5*sig_g,5*sig_g+dlam,dlam)
        gauss = f_gaussian(x,sig_g)
        kernel = gauss/np.trapz(gauss)

        convoluted = 1 + convolve(flux-1,kernel,mode='same')

        flux_norm = flux/convoluted

    # Iterative sigma clipping replacing outliers (cosmics) with nan
    for i in range(iter):
        flux_norm = np.where(abs(flux_norm-1) > sigclip*np.nanstd(flux_norm),np.nan,flux_norm)

    flux_clean = np.where(np.isnan(flux_norm),np.nan,flux)

    nans = np.isnan(flux_clean); x = lambda z: z.nonzero()[0]
    flux_clean[nans]= np.interp(x(nans),x(~nans),flux_clean[~nans])

    flux_clean = np.where((flux_clean > flux)|(abs(flux_clean-flux)<0.05),flux,flux_clean)

    return flux_clean


def cleanML(wave, flux, filename, axs, manual=False):

    # If manual, it allows to visually play with the sig clipping to remove cosmics
    if manual == True:
        next = 'n'; sigclip1 = 2.5; sigclip2 = 1.3 # 3.0 1.3
        print('Default sigma clipping values are 2.5, 1.3.')
        print('Current max/min:',round(flux.max(),1),round(flux.min(),1))

        while next == 'n':
            plt.figure(figsize=(13,4))
            plt.plot(wave,flux,c='orange',lw=.5,label='original')
            flux_iter1 = cosmicML(wave,flux,      method='zscore',sigclip=sigclip1,iter=3)
            flux_iter2 = cosmicML(wave,flux_iter1,method='zscore',sigclip=sigclip2,iter=3)
            new_flux = np.where((flux_iter1 > 1.01)&((wave < 6550)|(wave > 6576)),flux_iter2,flux_iter1)

            plt.plot(wave,new_flux,c='b',lw=.5,label='final')
            plt.legend(); plt.tight_layout(); plt.ylim(-0.1,2.1)#; plt.show(block=False)

            print('New max/min:',round(new_flux.max(),1),round(new_flux.min(),1))

            inp = input('Enter a new sigma clips e.g. "3.0 2" or hit return to continue: ')
            if inp == '': next = 'y'
            else:
                try: sigclip1,sigclip2 = [float(i) for i in inp.split()]
                except: print('Input is not valid')

            plt.close()
    else:
        sigclip1 = 2.5; sigclip2 = 1.3
        flux_iter1 = cosmicML(wave,flux,      method='zscore',sigclip=sigclip1,iter=3)
        flux_iter2 = cosmicML(wave,flux_iter1,method='zscore',sigclip=sigclip2,iter=3)
        new_flux = np.where((flux_iter1 > 1.01)&((wave < 6550)|(wave > 6576)),flux_iter2,flux_iter1)

    axs[1].plot(wave,flux,c='orange',lw=.2,alpha=.7,label='original*')
    i = 0
    for wl_em in [3967.79,4958.911,5006.843,6300.304,6548.04,6583.46,6716.44,6730.82]:
        mask = (wave > wl_em-0.8)&(wave < wl_em+0.8)
        new_flux[mask] = flux[mask]
        axs[1].plot(wave[mask],new_flux[mask],c='g',lw=.2,label='Em. lines'if i==0 else "",zorder=10)
        i += 1

    if new_flux.min() < 0: new_flux = np.where(new_flux < 0,0.0,new_flux)

    axs[1].plot(wave,new_flux,c='b',lw=.2,alpha=.7,label='final')

    return new_flux


def remove_wave(path=maindir+'tmp/',only_list='to_correct.txt'):
    if only_list != None:
        table = findtable(only_list, delimiter=' ')

    for file in os.listdir(path):
        if file.endswith('.ascii'):

            if only_list != None and not file in table['File']:
                continue

            data = Table.read(path+file,format='ascii',delimiter=' ')
            try:
                if max(data['col2']) > 1.5:
                    plt.plot(data['col1'],data['col2'],lw=.5)
                    plt.plot([3950,6850],[2,2],'k',lw=.5)
                    plt.show(block=False)
                    inp = input('Wanna print data for %s it? [y/ ]: ' % file)
                    if inp == 'y':
                        print(file,'>1.5',max(data['col2']))
                    plt.close()
                if min(data['col2']) < 0: print(file,'<0',min(data['col2']))
                data.remove_columns(['col1'])
                data.write(maindir+'tmp/new/'+file,format='ascii.no_header')
            except:
                print(file)

        else: continue


#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# Update tables with new data | / IMPLEMENT IN DB AT SOME POINT!!
#table1,table2 = findtable('OBAs_ML_raw.fits'),findtable('MAUI_results.fits')
#columns_to_update = table2.colnames[1:-9]+[table2.colnames[-1]]
#for row,i in zip(table1,range(len(table1))):
#    updated = False
#    if row['ID'] not in table2['ID']: continue
#    if row['Teff'] == table2[table2['ID'] == row['ID']]['Teff'].data[0]: continue
#    for col_name in columns_to_update:
#        table1[i][col_name] = table2[table2['ID'] == row['ID']][col_name].data[0]
#        if updated == False: print(row['ID']); updated = True
#table1.write(maindir+'tables/table1_updated.fits',format='fits',overwrite=True)

# This one is to empty bad data
#table1 = findtable('OBAs_ML_raw.fits')
#columns_to_update = [i for i in table1.columns[36:-1]]
#for row,i in zip(table1,range(len(table1))):
#    for col_name in columns_to_update:
#        if row[col_name] in ['d','<','>']:
#            for j in columns_to_update[columns_to_update.index(col_name)+1:columns_to_update.index(col_name)+4]:
#                table1[i][j] = np.nan
#table1.write(maindir+'tables/OBAs_ML_ver1.fits',format='fits',overwrite=True)

#table1 = findtable('OBAs_ML_ver1b.fits')
#for row,i in zip(table1,range(len(table1))):
#    if row['QSiIII'] < 3:
#        table1[i]['EWSiIII1'] = table1[i]['FWSiIII1'] = table1[i]['depSiIII1'] = np.nan
#        table1[i]['EWSiIII2'] = table1[i]['FWSiIII2'] = table1[i]['depSiIII2'] = np.nan
#        table1[i]['EWSiIII3'] = table1[i]['FWSiIII3'] = table1[i]['depSiIII3'] = np.nan
#    if row['QSiII'] < 3:
#        table1[i]['EWSiII'] = table1[i]['FWSiII'] = table1[i]['depSiII'] = np.nan
#    if row['QHb'] < 3:
#        table1[i]['EWHb'] = table1[i]['FWHb'] = table1[i]['FW14Hb'] = table1[i]['FW34Hb'] = table1[i]['depHb'] = table1[i]['gamma'] = np.nan
#table1.write(maindir+'tables/OBAs_ML_ver1.fits',format='fits',overwrite=True)

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# Replace values in table
#table = findtable('OBAs_ML_raw.fits')
#for row,i in zip(table,range(len(table))):
#    for column,j in zip(row,range(len(row))):
#        #if column == 'N': table[i][j] = '='
#        if str(column) == str(1e+20): table[i][j] = np.nan
#table.write(maindir+'tables/OBAs_ML_raw_.fits',format='fits',overwrite=True)
