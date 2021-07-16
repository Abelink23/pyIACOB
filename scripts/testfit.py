import sys; sys.path.append('../')

from spec import *

width = 60
rv0 = 0.416
tol = 150
func = 'vrg_H'
plot = 'y'
iter = 1

# SiIII 4552.622 | SiII 6347.11 | SiII 6371.37 | CII 6578.05 | MgII 4481.13 | Hb 4861.325
line = 4861.325

errors = False

#star = spec('maui_T20000lgf130v110_V25000.txt',SNR='best',txt=True)
star = spec('HD36959',SNR='best',rv0=rv0)

# Test contaminations:
# The size of the line-core must be manually set as # times the FWHM of the first
# iteration which is highly dependent on the function used and sometimes shifts
# the fitted line. Number of iterations must be 2 or more.
test = 0 # 1/2 Tested for HD2905 - Hg 4340.2
nFWHM = 1 # def is 1

fitsol = {'sol': 0, 'line': None, 'RV_A': None, 'RV_kms': None,
          'EW': None, 'FWHM': None, 'depth': None, 'q_fit': None}

t0 = time.time()
'''============================ Parameters =========================='''
# Catch line input containing more than one line
if type(line) == str and (',' in line or ' ' in line.strip()):
    print('Error: More than one line selected.\nExitting...')

line = float(line)

dlamb = line/star.resolution

# Maximum shift between the minimum of the fitted line and the tabulated value
tol_aa = float(tol)*(line)*1000/cte.c  # Changes km/s to angstroms

# Maximum FWHM allowed (should be up to 18 for H lines)
FWHM_max = 17


'''======== Dictionary and boundary limits for the functions ========'''
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


'''========================== Line fitting =========================='''
i = 0; width_i = width
while i < iter:

    # Extracting the window of the spectrum
    window = (star.wave >= line-width_i/2.) & (star.wave <= line+width_i/2.)
    if not any(window):
        print('Line %sA not in spectra.\n' % line)
        break
    flux = star.flux[window]; wave = star.wave[window]

    # The min FWHM will be defined by either 3 times the dlam, or 3/4 of the
    # minimum theoretical FWHM given the input line and resolution
    dlam_mean = (wave[-1]-wave[0])/(len(wave)-1)
    FWHM_min = np.max([3*dlam_mean,3/4*dlamb])

    # Auto-resampling
    #if dlam_mean >= 0.025 and not 'log' in star.file_name:
    #    factor = dlam_mean/0.025; star.resamp(factor)

    # Find regions of the continuum to use during normalization (4 iter)
    for j in range(4):
        if j == 0:
            mask_i = ~np.isnan(flux)
        else:
            mask_i = ~sigma_clip(flux/continuum_i,maxiters=None,
                      sigma_lower=1.5,sigma_upper=2.5).mask #1.4

        c_fit = np.poly1d(np.polyfit(wave[mask_i],flux[mask_i],1))
        continuum_i = c_fit(wave)

    # Final normalization of the iteration
    flux_norm_i = flux / continuum_i


    '''====================== Fitting the line ===================='''
    popt_i, pcov_i = curve_fit(fitfunc,wave,flux_norm_i,bounds=bounds)
    try:
        popt_i, pcov_i = curve_fit(fitfunc,wave,flux_norm_i,bounds=bounds)
        flux_fit_i = fitfunc(wave,*popt_i)

        if test in [1,2]:
            try:
                plt.cla(); plt.figure(tight_layout=True)
                plt.plot(wave,flux_fit_i,'g',lw=.5,label='current fitting')
                plt.plot(wave,flux_norm_i,'b',lw=.1,label='normalized flux')
                # Center of fitted line:
                line_f = wave[flux_fit_i == min(flux_fit_i)][0]
                # Define how much of the core you want to mask (at the beginning)
                core = nFWHM*FWHM
                mask_l = (wave >= line_f-core/2) & (wave <= line_f+core/2)

                if test == 1:
                    resid = flux_fit_i/flux_norm_i
                    plt.plot(wave,resid+.05,'grey',lw=.3)
                    mask_c = ~sigma_clip(resid,maxiters=None,
                              sigma_lower=2.2,sigma_upper=1.3).mask #2.5/1.3
                    # Add the continuum regions and such line core together
                    mask_l = np.where(mask_l==False,mask_c,True)

                elif test == 2:
                    # Add the continuum regions and such line core together
                    mask_l = np.where(mask_l==False,mask_i,True)

                plt.plot(wave,np.where(mask_l==False,1,np.nan)+0.01,'k',lw=.3)

                # a) Use the mask
                popt_a, pcov_a = curve_fit(fitfunc,wave[mask_l],flux_norm_i[mask_l],bounds=bounds)
                flux_fit_a = fitfunc(wave,*popt_a)
                plt.plot(wave,flux_fit_a,'-.r',lw=.3,label='masked fitting')

                # b) Use the cropped flux where masked values are interpolated
                flux_norm_b = np.where(mask_l==True,flux_norm_i,np.nan)
                plt.plot(wave,flux_norm_b,'-.k',lw=.4,label='flux to fit')
                nans = np.isnan(flux_norm_b); x = lambda z: z.nonzero()[0]
                flux_norm_b[nans]= np.interp(x(nans),x(~nans),flux_norm_b[~nans])
                popt_b, pcov_b = curve_fit(fitfunc,wave,flux_norm_b,bounds=bounds)
                flux_fit_b = fitfunc(wave,*popt_b)
                plt.plot(wave,flux_fit_b,'-.g',lw=.3,label='interp. fitting')

                plt.legend(fontsize=9)
                plt.ylim([0.6,1.1])
                completeName = os.path.join(maindir+"/tmp_plots/fitted_0.jpg")
                plt.savefig(completeName,format='jpg',dpi=300)
            except: pass

        # Error fittings
        if errors is True:
            sig_errEW = np.std(flux_norm_i[mask_i])
            # Upper error:
            try:
                popt_up = curve_fit(fitfunc,wave,flux_norm_i+sig_errEW,bounds=bounds)[0]
                flux_fit_up_i = fitfunc(wave,*popt_up)
            except: print('Bad up flux')
            # Lower error:
            try:
                popt_dw = curve_fit(fitfunc,wave,flux_norm_i-sig_errEW,bounds=bounds)[0]
                flux_fit_dw_i = fitfunc(wave,*popt_dw)
            except: print('Bad dw flux')

        # Calculate the empirical approximate FWHM
        medval = (max(flux_fit_i) + min(flux_fit_i))/2
        medpos = [np.where(flux_fit_i <= medval)[0][value] for value in (0,-1)]
        FWHM = round(wave[medpos[1]]-wave[medpos[0]],2)

        # Checking step results
        if FWHM_min < FWHM < FWHM_max:
            flux_norm = flux_norm_i; continuum = continuum_i; mask = mask_i
            flux_fit = flux_fit_i; popt = popt_i; pcov = pcov_i; width = width_i
            print(popt)
            i = i + 1; width_i = FWHM*7
            if errors is True:
                flux_fit_up = flux_fit_up_i; flux_fit_dw = flux_fit_dw_i
        elif FWHM < FWHM_min:
            print('WARNING: FWHM(%d) < minimum FWHM.' % FWHM); break
        elif FWHM > FWHM_max:
            print('WARNING: FWHM(%d) < maximum FWHM.' % FWHM); break

    except: break


'''======================= Checking final results ======================='''
window = (star.wave >= line-width/2.) & (star.wave <= line+width/2.)
flux = star.flux[window]; wave = star.wave[window]

if i == 0:
    sys.exit('Problem in spectrum %s\nLine %sA could not be fitted or does not exist.\n'
    % (star.file_name,line))

line_f = wave[flux_fit == min(flux_fit)][0]
if abs(line - line_f) > tol_aa:
    sys.exit('Line %sA found outside tolerance.\n' % line)

# Add the rv0 on top on the position of the line
#line_f = line_f + float(star.rv0)*float(line)*1000/cte.c

RV_A = round((line_f - line),3)
RV_kms = round(((line_f - line)/line)*cte.c/1000,1)
line_f = round(line_f,3)

print('Line %sA found at %.3fA -> RV: %.1fkm/s\n' % (line,line_f,RV_kms))


'''========================== Calculate the EWs ========================='''
# stackoverflow.com/questions/34075111/calculate-equivalent-width-using-python-code
EW = .5*abs(fsum((wave[wl-1]-wave[wl])*((1-flux_fit[wl-1]) +
     (1-flux_fit[wl])) for wl in range(1,len(flux_fit))))
EW = round(1000*EW)


'''======================== Calculate the EW ========================'''
if i > 0 and errors is True:
    EW_up = .5*abs(fsum((wave[wl-1]-wave[wl])*((1-flux_fit_up[wl-1]) \
               +(1-flux_fit_up[wl])) for wl in range(1,len(flux_fit_up))))
    EW_up = round(1000*EW_up - EW)

    EW_dw = .5*abs(fsum((wave[wl-1]-wave[wl])*((1-flux_fit_dw[wl-1]) \
               +(1-flux_fit_dw[wl])) for wl in range(1,len(flux_fit_dw))))
    EW_dw = round(1000*EW_dw - EW)
    print('Error in EW is:',EW_up,'/',EW_dw)


'''==================== Calculate the final FWHM ===================='''
medval = (max(flux_fit) + min(flux_fit))/2
medpos = [np.where(flux_fit <= medval)[0][value] for value in (0,-1)]

try: l_val = np.interp(medval,[flux_fit[medpos[0]],flux_fit[medpos[0]-1]],
                              [wave[medpos[0]],wave[medpos[0]-1]])
except: l_val = wave[medpos[0]]
try: r_val = np.interp(medval,[flux_fit[medpos[1]],flux_fit[medpos[1]+1]],
                              [wave[medpos[1]],wave[medpos[1]+1]])
except: r_val = wave[medpos[1]]
FWHM = round(r_val-l_val,2)


'''===================== Calculate the line depth ==================='''
depth = round(1-min(flux_fit),2)


'''=================== Calculate the SNR continuum =================='''
sigma_cont = np.std(flux_norm[mask])
snr = int(1/sigma_cont)


'''=========================== Quality value ========================'''
q_fit = 1/np.std(flux_norm[flux_fit<.995]/flux_fit[flux_fit<.995]) #simple
q_fit = round(q_fit,3)

q_fit = np.sqrt(np.diag(pcov))[0] # Estimate of error for the Area


'''========================== Find more lines ======================='''
flux_norm_2 = flux_norm/flux_fit; mirror_line = line - RV_A
if func == 'g':
    bounds = ([-1,mirror_line-2*tol_aa,0],
              [ 0,mirror_line+2*tol_aa,6])
try:
    popt,pcov = curve_fit(fitfunc,wave,flux_norm_2,bounds=bounds)
    flux_fit_2 = fitfunc(wave,*popt); amp = popt[0]
    line_fit_2 = (wave,flux_fit_2)

    EW_2 = .5*abs(fsum((line_fit_2[0][wl-1]-line_fit_2[0][wl])*((1-line_fit_2[1][wl-1]) \
                    +(1-line_fit_2[1][wl])) for wl in range(1,len(line_fit_2[1]))))
    EW_2 = round(1000*EW_2)

    line_f_2 = line_fit_2[0][line_fit_2[1].tolist().index(min(line_fit_2[1]))]
    line_f_2 = line_f_2 + float(star.rv0)*float(line)*1000/cte.c
    RV_A_2 = round((line_f_2-line),3)
    RV_kms_2 = round(((line_f_2-line)/line)*const.c/1000,3)
    line_f_2 = round(line_f_2,3)
    print(line_f_2,RV_A_2,'A ',RV_kms_2,'km/s ',EW_2,'mA')
except: None

print(time.time() - t0)
'''================================ Plot ================================'''
if plot in ['y','yes']:

    fig, ax = plt.subplots()

    ax.scatter([l_val,r_val],[medval,medval],c='b',marker='+')

    ax.plot(wave,flux,'orange',lw=.3)
    ax.plot(wave,continuum,'r',lw=.3)
    ax.plot(wave,flux_norm,'b',lw=.3)

    # The three fittings
    if i > 0 and errors is True:
        ax.plot(wave,flux_fit_up,'k',lw=.3)
        ax.plot(wave,flux_fit_dw,'k',lw=.3)
    ax.plot(wave,flux_fit,'g',lw=1)

    try:
        ax.plot(wave,np.where(mask   == False,1,np.nan)+0.01,'k',lw=.3)
        ax.plot(wave,np.where(mask_l == False,1,np.nan)+0.02,'k',lw=.3)
    except: pass

    # For more lines
    #ax.plot(wave,flux_norm/flux_fit,'grey',lw=.3)
    #try: ax.plot(line_fit_2[0],flux_fit_2,'purple',lw=1)
    #except: pass

    ax.set_title('%s | %.2f | RV: %d | EW: %d | FWHM: %.2f' %
                (star.id_star,line_f,RV_kms,EW,FWHM))

    #ax.set_yticks([])
    ax.set_xlabel('$\lambda$ $[\AA]$',size=13)
    ax.set_ylabel('Normalized flux',size=13)
    ax.tick_params(direction='in',top='on')
    ax.figure.subplots_adjust(top=.9,bottom=.12,right=.9,left=.1)

    completeName = os.path.join(maindir+"/tmp_plots/fitted.jpg")
    ax.figure.savefig(completeName,format='jpg',dpi=300)
    plt.show(block=False)


fitsol = {'sol': 1, 'line': line_f, 'RV_A': RV_A, 'RV_kms': RV_kms,
           'EW': EW, 'FWHM': FWHM, 'depth': depth, 'q_fit': q_fit}
for f_par,par in zip(fit_dic[func.split('_')[0]],popt): fitsol[f_par] = par

#fitsol['wave'] = wave; fitsol['flux'] = flux
#fitsol['flux_norm'] = flux_norm; fitsol['flux_fit'] = flux_fit

print(fitsol)
100*(0.1-0.013)

# Theorical FHWM:
#if   func == 'g': FWHM = 2*np.sqrt(2*np.log(2))*popt[2]
#elif func == 'l': FWHM = 2*abs(popt[2])
#elif func == 'v': FWHM = 2*(.5346*popt[3]+np.sqrt(.2166*(popt[3]**2)+popt[2]**2))
#elif func == 'r': FWHM = 1.7*popt[3]*line*1000/const.c
#FWHM = round(FWHM, 2)

popt
