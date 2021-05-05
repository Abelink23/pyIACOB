import sys; sys.path.append('../')

from spec import *

width = 15
offset = 0
tol = 150
func = 'vr_Z'
vsini = 120
plot = 'y'
iter = 3


# SiIII 4552.622 | SiII 6347.11 | SiII 6371.37 | CII 6578.05 | MgII 4481.13 | Hb 4861.325
line = 4552.622

errors = False

star = spec('maui_T24000lgf130v210_V25000.txt',SNR='best',txt=True)
#star = spec('HD55857',SNR='best')


t0 = time.time()
'''============================= Parameters ============================='''
try:
    if len(line.split(',')) > 1:
        print('More than one line was used as input.'); line = 0
except: None
if type(line) == str: line = float(line)

# Maximum shift between the minimum of the fitted line and the tabulated value
tol = float(tol)
tol_aa = tol*(line)*1000/cte.c  # Changes km/s to angstroms

dlamb = line/star.resolution
sigma,sig_min,sig_max = [0.8,dlamb/2/np.sqrt(2*np.log(2)),4]
# sig_min should be larger the minimum theoretical value dlam/2*sqrt(2log2)

# Maximum FWHM allowed (should be ~20 for H lines, ~13 for metalic lines)
FWHM_max = 20

'''=========== Set initial parameters for the chosen function ==========='''
inf = np.inf
# Fitting function: Gaussian | A,x0,sig
if func == 'g':
    fitfunc = f_gaussian1
    guess   =  [-0.1,line       ,sigma  ]
    bounds  = ([-inf,line-tol_aa,sig_min],
               [ 0. ,line+tol_aa,sig_max])

# Fitting function: Lorentzian | A,x0,gamma,y
elif func == 'l':
    fitfunc = f_lorentzian
    guess   =  [-0.1,line       , 0.5,1. ]
    bounds  = ([-inf,line-tol_aa,-inf,1. ],
               [ 0. ,line+tol_aa, inf,1.01])

# Fitting function: Voigt profile | A,x0,sigma,gamma,y
elif func == 'v':
    fitfunc = f_voigt
    bounds  = ([-.5,line-tol_aa,0,.0,1.  ],
               [ .0,line+tol_aa,2,.5,1.01])

# Fitting function: Rotational profile | A,x0,sigma,vsini
elif func == 'r':
    fitfunc = f_rot
    bounds  = ([.0,line-tol_aa,0. ,  1],
               [.3,line+tol_aa,2.5,410])

# Fitting function: Voigt x Rotational profile | A,x0,sigma,gamma,vsini,y
elif func == 'vr_H':
    fitfunc = f_voigtrot
    bounds  = ([-.5,line-tol_aa,  0,  0,  1,.0 ],
               [ .0,line+tol_aa,1.5,1.5,410,.01])
elif func == 'vr_Z':
    fitfunc = f_voigtrot
    bounds  = ([-.1,line-tol_aa,0. ,0,  1,.0 ],
               [ .0,line+tol_aa,1.5,1,410,.01])

# Fitting function: Voigt x Rotational + Gaussian profile
# A1,x0,sigma1,gamma,vsini,A2,sigma2,y
elif func == 'vrg_H':
    fitfunc = f_vrg
    bounds  = ([-.4,line-tol_aa, 0, 0,  1,-.07,0,.0 ],
               [ .0,line+tol_aa,10,10,410, .0 ,4,.01])
elif func == 'vrg_Z':
    fitfunc = f_vrg
    bounds  = ([-.1,line-tol_aa,0. ,0,  1,-1.3,0,-.01], # y=-.01 = larger EWs
               [ .0,line+tol_aa,1.5,1,410, 0. ,2, .01])
    #bounds  = ([-.3,line-tol_aa,0., 0.,  1,-2. ,0.,.0 ],
    #           [ .0,line+tol_aa,8.,10.,410, 0. ,4.,.01])

#popt

'''========================== Line fitting =========================='''
iterations = iter; i = 0; width_i = width
while i < iterations:

    '''============ Extracting the window of the spectrum ==========='''
    window = (star.wave >= line-width_i/2.) & (star.wave <= line+width_i/2.)
    if not any(window): print('Line %sA not in spectra.\n' % line); break
    flux = star.flux[window]; wave = star.wave[window]

    '''====================== Auto-resampling ======================='''
    #if star.dx >= 0.025 and not 'log' in star.file_name:
    #    factor = star.dx/0.025; star.resamp(factor)

    '''======== Find regions to exclude during normalization ========'''
    iter_norm = 4; mask_i = ~np.isnan(flux)
    for j in range(iter_norm):
        c0_fit = np.poly1d(np.polyfit(wave[mask_i],flux[mask_i],1))
        continuum_i = c0_fit(wave)

        if j < iter_norm: mask_i = ~sigma_clip(flux/continuum_i,\
            sigma_lower=1.4,sigma_upper=2.5,maxiters=None).mask

    '''===================== Final normalization ===================='''
    flux_norm_i = flux / continuum_i

    '''====================== Fitting the line/s ===================='''
    popt_i,pcov = curve_fit(fitfunc,wave,flux_norm_i,bounds=bounds)
    try:
        popt_i,pcov = curve_fit(fitfunc,wave,flux_norm_i,bounds=bounds)
        flux_fit_i = fitfunc(wave,*popt_i)

        #resid = flux_norm_i/flux_fit_i
        #std = np.std(resid); sigma = 0.6
        #mask_l = abs(resid-1)<sigma*std

        #popt_i,pcov = curve_fit(fitfunc,wave[mask_l],flux_norm_i[mask_l],guess,bounds=bounds)
        #flux_fit_i = fitfunc(wave,*popt_i)

        '''======================== Error fittings ======================'''
        if errors is True:
            sigma_errEW = abs(np.std(flux_norm_i[mask_i]))
            # Upper error:
            try: popt_up,pcov_up = curve_fit(fitfunc,wave,flux_norm_i+sigma_errEW,
                 guess,bounds=bounds); flux_fit_up_i = fitfunc(wave,*popt_up)
            except: print('Bad up flux')
            # Lower error:
            try: popt_dw,pcov_dw = curve_fit(fitfunc,wave,flux_norm_i-sigma_errEW,
                 guess,bounds=bounds); flux_fit_dw_i = fitfunc(wave,*popt_dw)
            except: print('Bad dw flux')

        '''====================== Calculate the FWHM ===================='''
        # Empirical approximate FWHM:
        medval = (max(flux_fit_i) + min(flux_fit_i))/2
        medpos = [np.where(flux_fit_i <= medval)[0][value] for value in (0,-1)]
        FWHM = round(wave[medpos[1]]-wave[medpos[0]],2)

        '''====================== Checking step results ====================='''
        if dlamb < FWHM < FWHM_max:

            flux_norm = flux_norm_i; continuum = continuum_i; mask = mask_i
            flux_fit = flux_fit_i; popt = popt_i; width = width_i; print(popt)

            if errors is True:
                flux_fit_up = flux_fit_up_i; flux_fit_dw = flux_fit_dw_i

            width_i = FWHM*7; i = i + 1

        else:
            if FWHM < dlamb: print('WARNING: FWHM<dlam')
            if FWHM > FWHM_max: print('WARNING: FWHM>%i' % FWHM_max)
            break

    except: break


'''======================= Checking final results ======================='''
window = (star.wave >= line-width/2.) & (star.wave <= line+width/2.)
flux = star.flux[window]; wave = star.wave[window]

if i == 0:
    sys.exit('Problem in spectrum %s\nLine %sA could not be fitted or does not exist.\n'
    % (star.file_name,line))

if FWHM >= 2 and not func in ['r','vr','vrg_H','vrg_Z']:
    print('FWHM > 2, consider switching to a model with rotation for',line)

fitted_line = wave[np.where(flux_fit == min(flux_fit))][0]
if abs(line - fitted_line) > tol_aa:
    sys.exit('Line %sA found outside tolerance.\n' % line)

fitted_line = fitted_line + star.offset
RV_angs = round((fitted_line - line),3)
RV_lamb = round(((fitted_line - line)/line)*cte.c/1000,3)
fitted_line = round(fitted_line,3)

print('Line %sA found at ' % (line) + str(round(fitted_line,2)) + \
      'A | RV: ' + str(RV_lamb) + ' [km/s] \n')


'''========================== Calculate the EWs ========================='''
# stackoverflow.com/questions/34075111/calculate-equivalent-width-using-python-code
EW = .5*abs(fsum((wave[wl-1]-wave[wl])*((1-flux_fit[wl-1]) \
            +(1-flux_fit[wl])) for wl in range(1,len(flux_fit))))
EW = round(1000*EW,2)


'''======================= Calculate the EW errors ======================'''
if i > 0 and errors is True:
    EW_up = .5*abs(fsum((wave[wl-1]-wave[wl])*((1-flux_fit_up[wl-1]) \
               +(1-flux_fit_up[wl])) for wl in range(1,len(flux_fit_up))))
    EW_up = round(1000*EW_up - EW,2)

    EW_dw = .5*abs(fsum((wave[wl-1]-wave[wl])*((1-flux_fit_dw[wl-1]) \
               +(1-flux_fit_dw[wl])) for wl in range(1,len(flux_fit_dw))))
    EW_dw = round(1000*EW_dw - EW,2)
    print('Error in EW is: ' + str(EW_up) + ' / ' + str(EW_dw))


'''================= Calculate the final FWHM interpolating ================='''
fwplt = [1,1]
medval = (max(flux_fit) + min(flux_fit))/2
medpos = [np.where(flux_fit <= medval)[0][value] for value in (0,-1)]

try: l_val = np.interp(medval,[flux_fit[medpos[0]],flux_fit[medpos[0]-1]],
                              [wave[medpos[0]],wave[medpos[0]-1]])
except: l_val = wave[medpos[0]]; fwplt[0] = 0
try: r_val = np.interp(medval,[flux_fit[medpos[1]],flux_fit[medpos[1]+1]],
                              [wave[medpos[1]],wave[medpos[1]+1]])
except: r_val = wave[medpos[1]]; fwplt[1] = 0
FWHM = round(r_val-l_val,2)


'''===================== Calculate the line depth ==================='''
depth = round(1-min(flux_fit),2)


'''===================== Calculate the SNR continuum ===================='''
sigma_cont = np.std(flux_norm[mask])
snr = int(1/sigma_cont)


'''============================ Quality value ==========================='''
q_fit = np.std(flux_norm[flux_fit<.995]/flux_fit[flux_fit<.995]) #simple
#q_fit = np.std(flux_norm[flux_fit<.995]/flux_fit[flux_fit<.995])/sigma_cont
q_fit = round(q_fit,3)


'''============================ Find more lines ========================='''
flux_norm_2 = flux_norm/flux_fit; mirror_line = line - RV_angs
if func == 'g':
    bounds = ([-inf,mirror_line-2*tol_aa,sig_min],[0,mirror_line+2*tol_aa,sig_max])
    guess = [-.1,mirror_line,sigma]
try:
    popt,pcov = curve_fit(fitfunc,wave,flux_norm_2,guess,bounds=bounds)
    flux_fit_2 = fitfunc(wave,*popt); amp = popt[0]
    line_fit_2 = (wave,flux_fit_2)

    EW_2 = .5*abs(fsum((line_fit_2[0][wl-1]-line_fit_2[0][wl])*((1-line_fit_2[1][wl-1]) \
                    +(1-line_fit_2[1][wl])) for wl in range(1,len(line_fit_2[1]))))
    EW_2 = round(1000*EW_2, 2)

    fitted_line_2 = line_fit_2[0][line_fit_2[1].tolist().index(min(line_fit_2[1]))]
    fitted_line_2 = fitted_line_2 + offset
    RV_angs_2 = round((fitted_line_2-line),3)
    RV_lamb_2 = round(((fitted_line_2-line)/line)*const.c/1000,3)
    fitted_line_2 = round(fitted_line_2,3)
    print(fitted_line_2,RV_angs_2,'A ',RV_lamb_2,'km/s ',EW_2,'mA')
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

    ax.plot(wave,np.where(mask == False,1,np.nan)+0.01,'k',lw=.3)
    #ax.plot(wave,np.where(mask_l == False,1,np.nan)+0.03,'k',lw=.4)

    #For more lines
    #ax.plot(wave,flux_norm/flux_fit,'grey',lw=.3)
    #try: ax.plot(line_fit_2[0],flux_fit_2,'purple',lw=1)
    #except: None

    ax.set_title(star.name_star + ' | ' + str(line) + ' | ' + 'RV: ' + str(RV_lamb)
     + ' | ' + 'EW: ' + str(EW) + ' | ' + 'FWHM: ' + str(FWHM))

    #ax.set_yticks([])
    ax.set_xlabel('$\lambda$ $[\AA]$',size=13)
    ax.set_ylabel('Normalized flux',size=13)
    ax.tick_params(direction='in',top='on')
    ax.figure.subplots_adjust(top=.9,bottom=.12,right=.9,left=.1)

    completeName = os.path.join(maindir+"/tmp_plots/fitted.jpg")
    ax.figure.savefig(completeName,format='jpg',dpi=300)
    plt.show(block=False)


print(fitted_line,RV_angs,RV_lamb,EW,FWHM,depth,snr,q_fit)


# Theorical FHWM:
#if   func == 'g': FWHM = 2*np.sqrt(2*np.log(2))*popt[2]
#elif func == 'l': FWHM = 2*abs(popt[2])
#elif func == 'v': FWHM = 2*(.5346*popt[3]+np.sqrt(.2166*(popt[3]**2)+popt[2]**2))
#elif func == 'r': FWHM = 1.7*popt[3]*line*1000/const.c
#FWHM = round(FWHM, 2)
