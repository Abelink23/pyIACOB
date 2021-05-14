from db import *

from math import fsum

import scipy.constants as cte
from scipy.special import wofz,erf
from scipy.optimize import curve_fit
from scipy.signal import convolve
from scipy.interpolate import interp1d

from astropy.time import Time
from astropy.stats import sigma_clip

import matplotlib.pyplot as plt


class spec():
    def __init__(self,spectrum,SNR=None,offset=0,txt=False):
        '''
        Parameters
        ----------

        spectrum : str
            Enter the input spectrum, either name(s) of the star(s), the fits files
            separated by coma, a .txt/.lst file containing the filenames, or '*'
            if you want to select all the fits files inside the working folder.

        SNR : str, optional
            If 'best' as input, it finds the best SNR spectrum for the given name.
            If 'bestMF' same as 'best' but prioritizing spectra from HERMES/FEROS.

        offset : float, optional
            Enter the offset in wavelength [AA] of the spectrum to plot. Default is 0.

        txt : boolean, optional
            If True, it assumes spectrum from a two-columns file with wavelenght and flux.

        '''

        if type(spectrum) == list:
            if len(spectrum) > 1:
                print('Error: More than one spectrum selected.\nExitting...')
                return None
            else: spectrum = spectrum[0]

        if SNR in ['best','bestMF'] and not '.fits' in spectrum and not txt == True:
            try: self.spectrum = findstar(spectrum,SNR=SNR)[0]
            except: self.spectrum = None
        else: self.spectrum = spectrum

        self.filename = self.spectrum.split('/')[-1]
        self.id_star = self.spectrum.split('/')[-1].split('_')[0]
        self.resolution = int(re.split('(\d*\d+)',self.spectrum)[-2])
        self.offset = offset # Note, run self.waveflux to apply offset.
        if txt == False: self.waveflux()
        elif txt == True: self.txtwaveflux()
        self.txt = txt


    def spc(self):
        try:
            query = Simbad.query_object(self.id_star)

            if query == None and 'HD' in self.id_star:
                new_id_star = self.id_star.replace('HD','HD ')
                query = Simbad.query_object(new_id_star)
            self.SpC = query['SP_TYPE'][0]
        except:
            self.SpC = ''


    def waveflux(self,lwl=None,rwl=None,width=0,helcorr='hel'):
        '''
        Parameters
        ----------
        lwl : float, optional
            Sets the start wavelenght of the spectrum.

        rwl : float, optional
            Sets the end wavelenght of the spectrum.

        width : int, optional
            Sets the width in [AA] where the line fits well in. Default is 10.

        helcorr : str, optional
            If 'hel' as input (default), it applies the heliocentric correction.

        '''

        # Retrieve the key values fron the fits header
        hdu = fits.open(self.spectrum)  # Open the fits image file
        header0 = hdu[0].header         # Read header of primary extension

        instrum = header0['INSTRUME']   # Instrument

        x0 = header0['CRVAL1']          # Get the wavelenght of the first pixel
        dx = header0['CDELT1']          # Step of increase in wavelength
        pix0 = header0['CRPIX1']        # Reference pixel (generally 1, FEROS -49)
        spec_length = header0['NAXIS1'] # Length of the spectrum
        # Alternatively use len(hdu[0].data[0]) (NOT/MERCATOR) or len(hdu[0].data)

        # Correct Mercator CRVAL1 20101018-19:
        if any(bad in self.spectrum for bad in ['_20101018_','_20101019_']) \
               and x0 == 3763.9375: x0 = 3763.61

        try: vbar = header0['I-VBAR'] # [km/s] Barycent. rv correction at midpoint
        #    vbar = header0['BVCOR']  # [km/s] Barycent. rv correction at midpoint | MERCATOR
        #    vbar = header0['VHELIO'] # [km/s] Barycent. rv correction at midpoint | NOT
        except: print('No helio/bary-centric correction applied to' + self.spectrum); vbar = 0

        self.vbar = vbar

        try: hjd = header0['I-HJD']   # Heliocentric Julian date at midpoint
        except: hjd = Time(header0['DATE'],scale='utc').jd

        self.hjd = hjd

        try: I_SNR = header0['I-SNR']   # SNR from header
        except: I_SNR = np.nan

        self.snr = I_SNR

        # Make lists with wavelenght and flux for each spectrum
        if width >= 200: width = 200; \
            print('\nWARNING: Width value %f is too large, setting it to 200. ' %width)

        wave = x0 + dx*(np.arange(spec_length) - pix0 + 1)
        if '_log' in self.spectrum: wave = np.exp(wave)
        elif helcorr == 'hel' and not instrum == 'FEROS':
            wave = wave*(1 + 1000*vbar/cte.c)
        # Those with log and those from FEROS are already corrected from helcorr

        wave = wave - self.offset

        try: flux = hdu[0].data[0]
        except: flux = hdu[0].data

        if lwl != None and rwl != None:
            if wave[0] > lwl+dx or wave[-1] < rwl-dx:
                print('WARNING: Wavelenght limits outside spectrum wavelenght range.')
            flux = flux[(wave >= lwl-width/2.)&(wave <= rwl+width/2.)]
            wave = wave[(wave >= lwl-width/2.)&(wave <= rwl+width/2.)]

        if '_log' in self.spectrum:
            self.dx = (wave[-1]-wave[0])/(len(wave)-1)
        else: self.dx = dx

        hdu.close()

        self.wave_0 = wave; self.flux_0 = flux
        self.wave = wave; self.flux = flux

        return wave,flux,hjd


    def txtwaveflux(self,lwl=None,rwl=None,width=0):
        '''
        Parameters
        ----------
        See help for spec.waveflux

        '''
        data = findtable(self.spectrum,path=datadir+'ASCII/')

        try: wave = np.asarray(data['wavelenght']); flux = np.asarray(data['flux'])
        except: wave = np.asarray(data['col1']); flux = np.asarray(data['col2'])

        wave = wave - self.offset

        if lwl != None and rwl != None:
            dx = (wave[-1]-wave[0])/(len(wave)-1)
            if wave[0] > lwl+dx or wave[-1] < rwl-dx:
                print('WARNING: Wavelenght limits outside spectrum wavelenght range.')
            flux = flux[(wave >= lwl-width/2.)&(wave <= rwl+width/2.)]
            wave = wave[(wave >= lwl-width/2.)&(wave <= rwl+width/2.)]

        self.dx = (wave[-1]-wave[0])/len(wave); self.vbar = 0; self.hjd = 0

        self.wave_0 = wave; self.flux_0 = flux
        self.wave = wave; self.flux = flux

        return wave,flux,0


    def fitline(self,line,width=15,tol=150.,func='g',iter=3,info=False,outfit=False,plot=False):
        '''
        Parameters
        ----------
        line : float
            Sets the central wavelenght of the line to search and fit.

        width : int, optional
            Sets the width in [AA] where the line fits well in. Default is 15.

        tol : int, optional
            Sets the tolerance [km/s] to shifting the spectrum in order to fit the line.

        func : str, optional
            Choose the function to fit the line:
            'g' Gaussian (default); 'l' Lorentzian; 'v' Voigt; 'r' Rotational.

        iter : int, optional
            Number of iterations to optimize window width. Default is 3.

        info : boolean, optional
            If 'True', it will print information for each line fitting.

        plot : boolean, optional
            If 'True', it will create and show the plots.

        Returns
        -------
        Dictionary containing the results, fitting parameters and optionally the
        original line, normalized flux and fitted flux.

        Notes
        -----
        Emission lines are excluded, see "Filtering emission lines" section.
        Returns: Parameters from the fitted line and a last value containing data
        with input wavelenght and flux, plus the flux normalized, the flux of the
        fitted line, and the fitting parameters.
        '''

        fitsol = {'sol': 0, 'line': np.nan, 'RV_A': np.nan, 'RV_kms': np.nan,
                  'EW': np.nan, 'FWHM': np.nan, 'depth': np.nan, 'q_fit': np.nan}

        '''============================ Parameters =========================='''
        # Catch line input containing more than one line
        if type(line) == str and (',' in line or ' ' in line.strip()):
            print('Error: More than one line selected.\nExitting...')
            return fitsol

        line = float(line)

        dlamb = line/self.resolution

        # Maximum shift between the minimum of the fitted line and the tabulated value
        tol_aa = float(tol)*(line)*1000/cte.c  # Changes km/s to angstroms

        # Maximum FWHM allowed (should be up to 18 for H lines)
        FWHM_max = 17

        #=======================================================================
        #=========== Dictionary and boundary limits for the functions ==========
        fit_dic = {'g': ('A','x0','sigma'),
                   'l': ('A','x0','gamma','y'),
                   'v': ('A','x0','sigma','gamma','y'),
                   'r': ('A','x0','sigma','vsini'),
                  'vr': ('A','x0','sigma','gamma','vsini','y'),
                 'vrg': ('A','x0','sigma','gamma','vsini','A2','sigma2','y')}

        # Fitting function: Gaussian | A,x0,sig
        if func == 'g':
            fitfunc = f_gaussian1
            bounds  = ([-1,line-tol_aa,0],
                       [ 0,line+tol_aa,6])

        # Fitting function: Lorentzian | A,x0,gamma,y
        elif func == 'l':
            fitfunc = f_lorentzian
            bounds  = ([-1,line-tol_aa, 0,1. ],
                       [ 0,line+tol_aa,10,1.01])

        # Fitting function: Voigt profile | A,x0,sigma,gamma,y
        elif func == 'v':
            fitfunc = f_voigt
            bounds  = ([-10,line-tol_aa,0. ,0. ,1.  ], #'A' ~15 for As
                       [  0,line+tol_aa,7.5,8.5,1.01])

        # Fitting function: Rotational profile | A,x0,sigma,vsini
        elif func == 'r':
            fitfunc = f_rot
            bounds  = ([.0,line-tol_aa,0. ,  1],
                       [.3,line+tol_aa,2.5,410])

        # Fitting function: Voigt x Rotational profile | A,x0,sigma,gamma,vsini,y
        elif func == 'vr_H':
            fitfunc = f_voigtrot
            bounds  = ([-.3,line-tol_aa,0. ,0,  1,.0 ],
                       [ .0,line+tol_aa,4.5,8,410,.01])
        elif func == 'vr_Z':
            fitfunc = f_voigtrot
            bounds  = ([-.3,line-tol_aa,0. ,0,  1,.0 ],
                       [ .0,line+tol_aa,1.5,1,410,.01])

        # Fitting function: Voigt x Rotational + Gaussian profile
        # A,x0,sigma,gamma,vsini,A2,sigma2,y
        elif func == 'vrg_H':
            fitfunc = f_vrg
            bounds  = ([-.5,line-tol_aa, 0, 0,  1,-.07,0,.0 ],
                       [ .0,line+tol_aa,10,10,410, .0 ,4,.01])
        elif func == 'vrg_Z':
            fitfunc = f_vrg
            bounds  = ([-.1,line-tol_aa,0. ,0,  1,-1.3,0,-.01], # y=-.01 = larger EWs
                       [ .0,line+tol_aa,1.5,1,410, 0. ,2, .01])

        #=======================================================================
        #============================= Line fitting ============================
        i = 0; width_i = width
        while i < iter:

            # Extracting the window of the spectrum
            window = (self.wave >= line-width_i/2) & (self.wave <= line+width_i/2)
            if not any(window):
                print('Line %sA not in spectra.\n' % line)
                return fitsol
            flux = self.flux[window]; wave = self.wave[window]

            dx_mean = (wave[-1]-wave[0])/(len(wave)-1)
            # Auto-resampling
            #if dx_mean >= 0.025 and not 'log' in star.filename:
            #    factor = dx_mean/0.025; star.resamp(factor)

            # Find regions of the continuum to use during normalization (4 iter)
            for j in range(4):
                if j == 0:
                    mask_i = ~np.isnan(flux)
                else:
                    mask_i = ~sigma_clip(flux/continuum_i,maxiters=None,
                              sigma_lower=1.5,sigma_upper=2.5).mask #1.4 before

                c_fit = np.poly1d(np.polyfit(wave[mask_i],flux[mask_i],1))
                continuum_i = c_fit(wave)

            # Final normalization of the iteration
            flux_norm_i = flux / continuum_i


            #===================================================================
            #========================= Fitting the line ========================
            try:
                popt_i = curve_fit(fitfunc,wave,flux_norm_i,bounds=bounds)[0]
                flux_fit_i = fitfunc(wave,*popt_i)

                # Calculate the empirical approximate FWHM
                medval = (max(flux_fit_i) + min(flux_fit_i))/2
                medpos = [np.where(flux_fit_i <= medval)[0][value] for value in (0,-1)]
                FWHM = round(wave[medpos[1]]-wave[medpos[0]],2)

                # Checking step results
                # The min FWHM will be defined by either 3 times the dx, or 3/4 of the
                # minimum theoretical FWHM given the input line and resolution
                FWHM_min = np.max([3*dx_mean,3/4*dlamb])
                if FWHM_min < FWHM < FWHM_max:
                    flux_norm = flux_norm_i; continuum = continuum_i; mask = mask_i
                    flux_fit = flux_fit_i; popt = popt_i; width = width_i
                    i = i + 1; width_i = FWHM*7
                elif FWHM < FWHM_min:
                    print('WARNING: FWHM(%d) < minimum FWHM.' % FWHM); break
                elif FWHM > FWHM_max:
                    print('WARNING: FWHM(%d) < maximum FWHM.' % FWHM); break

            except: break

        #=======================================================================
        #======================== Checking final results =======================
        window = (self.wave >= line-width/2.) & (self.wave <= line+width/2.)
        flux = self.flux[window]; wave = self.wave[window]

        if i == 0:
            if info is True:
                print('Problem in spectrum %s' % self.filename)
                print('Line %sA could not be fitted or does not exist.\n' % line)
            return fitsol

        line_f = wave[flux_fit == min(flux_fit)][0]
        if abs(line - line_f) > tol_aa:
            if info is True:
                print('Line %sA found outside tolerance.\n' % line)
            return fitsol

        line_f = line_f + self.offset
        RV_A   = round((line_f - line),3)
        RV_kms = round(((line_f - line)/line)*cte.c/1000,1)
        line_f = round(line_f,3)

        if info is True:
            print('Line %sA found at %.3fA -> RV: %.1fkm/s\n' % (line,line_f,RV_kms))

        #=======================================================================
        #=========================== Calculate the EW ==========================
        # stackoverflow.com/questions/34075111/calculate-equivalent-width-using-python-code
        EW = .5*abs(fsum((wave[wl-1]-wave[wl])*((1-flux_fit[wl-1]) +
             (1-flux_fit[wl])) for wl in range(1,len(flux_fit))))
        EW = round(1000*EW)

        #=======================================================================
        #====================== Calculate the final FWHM =======================
        medval = (max(flux_fit) + min(flux_fit))/2
        medpos = [np.where(flux_fit <= medval)[0][value] for value in (0,-1)]
        try: l_val = np.interp(medval,[flux_fit[medpos[0]],flux_fit[medpos[0]-1]],
                                      [wave[medpos[0]],wave[medpos[0]-1]])
        except: l_val = wave[medpos[0]]
        try: r_val = np.interp(medval,[flux_fit[medpos[1]],flux_fit[medpos[1]+1]],
                                      [wave[medpos[1]],wave[medpos[1]+1]])
        except: r_val = wave[medpos[1]]
        FWHM = round(r_val-l_val,2)

        #=======================================================================
        #======================= Calculate the line depth ======================
        depth = round(1-min(flux_fit),2)


        #=======================================================================
        #===================== Calculate the SNR continuum =====================
        sigma_cont = np.std(flux_norm[mask])
        snr = int(1/sigma_cont)

        #=======================================================================
        #============================= Quality value ===========================
        q_fit = 1/np.std(flux_norm[flux_fit<.995]/flux_fit[flux_fit<.995]) #simple
        q_fit = round(q_fit,3)

        #=======================================================================
        #================================ Plot =================================
        if plot is True:

            fig, ax = plt.subplots()

            ax.plot(wave,flux,'orange',lw=.5)
            ax.plot(wave,continuum,'r',lw=.5)
            ax.plot(wave,flux_norm,'b',lw=.5)

            ax.plot(wave,flux_fit,'g',lw=.5)

            ax.plot(wave,np.where(mask==False,1,np.nan)+0.01,'k',lw=.5)

            ax.set_title('%s | %.2f | RV: %d | EW: %d | FWHM: %.2f' %
                     (self.id_star,line_f,RV_kms,EW,FWHM))

            ax.set_yticks([])
            ax.set_xlabel('$\lambda$ $[\AA]$',size=13)
            ax.set_ylabel('Normalized flux',size=13)
            ax.tick_params(direction='in',top='on')
            ax.figure.subplots_adjust(top=.9,bottom=.12,right=.88,left=.08)

        plt.show(block=False)

        #=======================================================================
        #================= Packing the results in a dictionary =================
        fitsol = {'sol': 1, 'line': line_f, 'RV_A': RV_A, 'RV_kms': RV_kms,
                   'EW': EW, 'FWHM': FWHM, 'depth': depth, 'q_fit': q_fit}
        for f_par,par in zip(fit_dic[func.split('_')[0]],popt):
            fitsol[f_par] = round(par,3)

        if outfit == True:
            fitsol['wave'] = wave; fitsol['flux'] = flux
            fitsol['flux_norm'] = flux_norm; fitsol['flux_fit'] = flux_fit

        return fitsol

        # Theorerical FWHM:
        #if   func == 'g': jFWHM = 2*np.sqrt(2*np.log(2))*popt[2]
        #elif func == 'l': jFWHM = 2*abs(popt[2])
        #elif func == 'v': jFWHM = 2*(.5346*popt[3]+np.sqrt(.2166*(popt[3]**2)+popt[2]**2))
        #elif func == 'r': jFWHM = 1.7*popt[3]*line*1000/cte.c
        #jFWHM = round(jFWHM, 2)


    def snrcalc(self,zone='v'):
        '''
        Parameters
        ----------
        zone : str, optional
            Select the zone to calculate the spectra.
                'b'/'B'     -> 4000-5000 A
                'v'/'V'     -> 5000-6000 A
                'r'/'R'     -> 6000-7000 A
                'all'/'ALL' -> 4000-7000 A

        Returns: Measured signal-to-noise ratio value.
        '''

        if zone in ['b','B']:
            mask = (self.wave > 4000) & (self.wave < 5000)
        elif zone in ['v','V']:
            mask = (self.wave > 5000) & (self.wave < 6000)
        elif zone in ['r','R']:
            mask = (self.wave > 6000) & (self.wave < 7000)
        elif zone in ['all','ALL']:
            mask = (self.wave > 4000) & (self.wave < 7000)

        lambda0 = np.mean(self.wave[mask]); resol = 10000

        sigma = lambda0/(2.35482*float(resol))

        gauss = f_gaussian(np.arange(-5*sigma,5*sigma,self.dx),sigma)
        kernel = gauss/np.trapz(gauss)

        convoluted = 1 + convolve(self.flux[mask]-1,kernel,mode='same')

        flux_norm = self.flux[mask]/convoluted

        snr_all = []
        for gap in findlist('snr_gaps.txt'):
            lwl,rwl = [float(i) for i in gap.split('-')]

            flux_norm_i = flux_norm[(self.wave[mask] >= lwl) & (self.wave[mask] <= rwl)]

            std = np.std(flux_norm_i); sig_clip = 3
            flux_cleaned = np.where(abs(flux_norm_i-1) > sig_clip*std,np.nan,flux_norm_i)

            snr_all.append(1/np.nanstd(flux_cleaned))

        self.snr = np.nanmean(snr_all)

        return int(round(np.nanmean(snr_all)))


    def cosmic(self,sigclip=1.5):
        '''
        Parameters
        ----------

        sigclip : float, optional
            Sigma clipping value used to remove rays. Default is 1.5.

        Returns: None (but the flux is replaced and cleaned from rays).
        '''

        lambda0 = np.mean(self.wave)

        # Two times the theoretical sigma offers better results
        sigma = 2*lambda0/(2.35482*float(self.resolution))

        x = np.arange(-5*sigma,5*sigma+self.dx,self.dx)
        gauss = f_gaussian(x,sigma)
        kernel = gauss/np.trapz(gauss)

        convoluted = 1 + convolve(self.flux-1,kernel,mode='same')

        flux_norm = self.flux/convoluted

        std = np.nanstd(flux_norm)
        flux_cleaned = np.where(flux_norm > 1 + sigclip*std,np.nan,self.flux)

        nans = np.isnan(flux_cleaned); x = lambda z: z.nonzero()[0]
        flux_cleaned[nans]= np.interp(x(nans),x(~nans),flux_cleaned[~nans])

        self.flux = flux_cleaned

        return None


    def degrade(self,resol=None,profile='g',vsini=None,vmac=None):
        '''
        Parameters
        ----------
        resol : int/float, optional
            Resolution of the gaussian profile used to degrade the spectrum.

        profile : str
            Use 'g' for gaussian profile convolution (Default).
            Use 'rotmac' for rotational+macroturbulence profile convolution.

        vsini : int/float, optiomal
            Value of vsini. Only valid for rotational+macroturbulence profile.

        vmac : int/float, optiomal
            Value of vmac. Only valid for rotational+macroturbulence profile.

        Returns: None (but the flux is replaced by the degraded one).
        '''

        lambda0 = np.mean(self.wave)

        if profile == 'g':
            sigma = lambda0/(2.35482*float(resol))

            x = np.arange(-10*sigma,10*sigma+self.dx,self.dx)
            gauss = f_gaussian(x,sigma)
            kernel = gauss/np.trapz(gauss)
            self.resolution = resol

        elif profile == 'rotmac':
            x = np.arange(-9,9+self.dx,self.dx)
            rotmac = f_rotmac(x,lambda0,vsini,vmac)
            kernel = rotmac/np.trapz(rotmac)

        convoluted = 1 + convolve(self.flux-1,kernel,mode='same')

        self.flux = convoluted


    def resamp(self,dx,lwl=None,rwl=None,method='linear'):
        '''
        Parameters
        ----------

        dx : float/int
            New delta lambda to be used for the output spectra.

        lwl : float/int, optional
            Enter the forced initial wavelenght to be used during interpolation.
            If None, the original initial wavelenght will be used.

        rwl : float/int, optional
            Enter the forced final wavelenght to be used during interpolation.
            If None, the original final wavelenght will be used.

        method : str, optional
            Enter the interpolation method to be used. See doc for np.interp1d.
            Default is 'linear'.

        Returns: None (but the spectrum (wavelenght,flux) is resampled).
        '''

        try: float(dx)
        except: print('Input should be float or integrer.'); return None

        self.dx = dx

        if dx > np.mean(self.wave)/self.resolution/3:
            # It is divided by 3 to at least have 3 pixels in a gaussian
            print('WARNING: The new delta lambda implies lossing information...')

        if lwl == None and rwl == None:
            lwl = self.wave[0]; rwl = self.wave[-1]

        f = interp1d(self.wave,self.flux,kind=method,fill_value="extrapolate")
        self.wave = np.arange(lwl,rwl+self.dx,self.dx)
        self.flux = f(self.wave)

        return None


    def export(self,tail='',extension='.dat'):
        filename = self.filename.replace('.fits','')
        np.savetxt(maindir+'tmp/%s' % (filename + tail + extension),
                   np.c_[self.wave,self.flux],fmt=('%.4f','%.6f'))


    def plotline(self,lines,width=10,ylim=None):
        '''
        Parameters
        ----------

        lines : float, str
            Enter the wavelenght(s) of the line(s) to plot, either in a coma-separated
            string, or in a .txt/.lst file containing the lines.

        width : int, optional
            Sets the width in [AA] where the line fits well in. Default is 10.

        ylim : tuple/list, optional
            Sets the y-limits for the plot. Input must be like "[ymin,ymax]".

        Returns: None (but the plots are generated).
        '''

        self.spc()

        lines,elements,_ = findlines(lines)
        if len(lines) > 1: nrows = ncols = int(round(np.sqrt(len(lines)),0))
        else: nrows = ncols = 1

        for line,element,nplot in zip(lines,elements,range(len(lines))):

            mask = (self.wave > line - width/2) & (self.wave < line + width/2)

            if len(lines) > 1:
                plt.subplot(nrows,ncols,nplot+1)
                plt.xticks([round(line-width/3,1),round(line,1),round(line+width/3,1)])
                plt.title(element,fontsize=6,pad=1)

            plt.plot(self.wave[mask],self.flux[mask],lw=.3,label=self.id_star+' '+self.SpC)
            plt.tick_params(direction='in',top='on')

            if ylim is not None and (type(ylim) is list or type(ylim) is tuple): plt.ylim(ylim)

            if len(lines) == 1:
                plt.xlabel('$\lambda$ $[\AA]$',size=13)
                plt.ylabel('Normalized flux',size=13)

            plt.tight_layout()

        plt.legend(); plt.show(block=False)

        return None


    def plotspec(self,lwl=3800,rwl=8000,poslines=None,ylim=None):
        '''
        Parameters
        ----------
        lwl : float, optional
            Sets the start wavelenght of the spectrum.

        rwl : float, optional
            Sets the end wavelenght of the spectrum.

        poslines : str, optional
            If 'all' or 'OB', it will overplot position of spectral lines.

        ylim : tuple/list, optional
            Sets the y-limits for the plot. Input must be like "[ymin,ymax]".

        Returns: None (but the plots are generated).
        '''

        self.spc()

        if lwl < min(self.wave): lwl = min(self.wave)
        if rwl > max(self.wave): rwl = max(self.wave)

        mask = (self.wave > lwl) & (self.wave < rwl)

        if poslines in ['all','OB']:
            if poslines == 'all': table = findtable('ALL_all.txt',delimiter=',')
            elif poslines == 'OB': table = findtable('ALL_OBs_n4+.txt',delimiter=',')

            synlines = table['wl_air']; elements = table['spc']; gfs = table['-lg(gf)']

            # Aqui falta definir mejor los constrains para plotear lineas loggf por ejemplo
            for synline,element,gf in zip(synlines,elements,gfs):

                if synline + self.offset < lwl or \
                   synline + self.offset > rwl or gf <= -0.5: continue
                else:
                    try: depth = max(self.flux[mask]) - min(self.flux[mask]) # or 1-min
                    except: print('Problem finding max/min in masked flux.'); return None
                    # depth line mask = depth deepest line
                    plt.text(synline-.1,1.025,element,size=6,rotation=75)
                    plt.plot([synline,synline],[np.mean(self.flux[mask])-depth,
                             np.mean(self.flux[mask])],'k',lw=10**gf/5)
                    # 10**gf/5 empiric way to draw thicker lines for instense lines

        plt.plot(self.wave[mask],self.flux[mask],lw=.3,label=self.id_star+' '+self.SpC)
        plt.tick_params(direction='in',top='on')

        if ylim is not None and (type(ylim) is list or type(ylim) is tuple): plt.ylim(ylim)

        plt.xlabel('$\lambda$ $[\AA]$',size=13)
        plt.ylabel('Normalized flux',size=13)

        plt.legend(); plt.tight_layout(); plt.show(block=False)

        return None


def f_gaussian(x,sigma):
    return np.exp(-(x/sigma)**2/2)

def f_gaussian1(x,A,x0,sigma):
    # A -> Amplitude;  x0 -> center
    return A*np.exp(-(x-x0)**2/(2*sigma**2)) + 1

def f_lorentzian(x,A,x0,gamma,y):
    return A*gamma**2/((x-x0)**2 + gamma**2) + y

def f_voigt(x,A,x0,sigma,gamma,y):
    # sigma -> gaussian width; gamma -> lorentzian width
    # sigma = alpha / sqrt(2 * np.log(2))
    return A*np.real(wofz((x-x0+1j*gamma)/sigma/np.sqrt(2)))/sigma/np.sqrt(2*np.pi) + y

def f_rot(x,A,x0,sigma,vsini):
    G = A*np.exp(-(x - x0)**2/(2*sigma**2))

    # Default value: beta=1.5 (epsilon=0.6) beta= epsilon/(1-epsilon)
    eps = 0.6; delta = 1000*x0*vsini/cte.c; doppl = 1 - ((x-x0)/delta)**2
    R = A*(2*(1-eps)*np.sqrt(doppl) + np.pi*eps/2.*doppl)/(np.pi*delta*(1-eps/3))
    R = np.nan_to_num(R)

    return 1-convolve(G,R,mode='same')

def f_voigtrot(x,A,x0,sigma,gamma,vsini,y):
    V = A*np.real(wofz((x-x0+1j*gamma)/sigma/np.sqrt(2)))/sigma/np.sqrt(2*np.pi) + y

    eps = 0.6; delta = 1000*x0*vsini/cte.c; doppl = 1-((x-x0)/delta)**2
    R = A*(2*(1-eps)*np.sqrt(doppl) + np.pi*eps/2.*doppl)/(np.pi*delta*(1-eps/3))
    R = np.nan_to_num(R)

    return 1-convolve(V,R,mode='same')

def f_vrg(x,A,x0,sigma,gamma,vsini,A2,sigma2,y):
    VG = A*np.real(wofz((x-x0+1j*gamma)/sigma/np.sqrt(2)))/sigma/np.sqrt(2*np.pi) + y\
         + A2*np.exp(-(x-x0)**2/(2*sigma2**2))

    eps = 0.6; delta = 1000*x0*vsini/cte.c; doppl = 1-((x-x0)/delta)**2
    R = A*(2*(1-eps)*np.sqrt(doppl)+np.pi*eps/2.*doppl)/(np.pi*delta*(1-eps/3))
    R = np.nan_to_num(R)

    return 1-convolve(VG,R,mode='same')

def f_rotmac(x,x0,vsini=None,vmac=None):

    if vsini != None:
        # Rotational function:
        delta_R = 1000*x0*vsini/cte.c
        doppl = 1 - (x/delta_R)**2

        eps = 0.6
        R = (2*(1-eps)*np.sqrt(doppl) + np.pi*eps/2.*doppl)/(np.pi*delta_R*(1-eps/3))
        R = np.nan_to_num(R)

        if vmac == None: return R

    if vmac != None:
        # Macroturbulence function:
        delta_M = 1000*x0*vmac/cte.c
        A = 2/np.sqrt(np.pi)/delta_M

        x_2 = x[len(x)//2:]
        x_d = x_2/delta_M

        M_T = A*x_d*(-np.sqrt(np.pi)+np.exp(-x_d**2)/x_d+np.sqrt(np.pi)*erf(x_d))

        M = M_T # + M_R

        M = np.concatenate((M[::-1],M[1:]))

        if vsini == None: return M

    if vsini != None and vmac != None: return convolve(R,M,mode='same')
