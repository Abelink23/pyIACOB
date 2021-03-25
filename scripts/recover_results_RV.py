'''Script to calculate the RV mean values and their errors'''

import numpy as np

with open('/home/abelink/PhD/radial_velocity/HD36861_RV.csv','r') as file:
    list = file.read().splitlines()

mean_rv = []; mean_rv_error = []; mbjd = [] #, num_lines, snr, snr_std, spectrum
for spectra in list[1:]:
    spectra = spectra.split(',')
    mean_rv.append(float(spectra[0]))
    mean_rv_error.append(float(spectra[1]))
    mbjd.append(float(spectra[2]))


RVs_mean_all_means = np.mean(mean_rv)    # Mean of all spectra
RVs_mean_all_means_std = np.std(mean_rv) # Std of the mean of all spectra

peak_to_peak = abs(max(mean_rv) - min(mean_rv))
peak_to_peak_err = np.sqrt(mean_rv_error[mean_rv.index(max(mean_rv))]**2 + \
                        mean_rv_error[mean_rv.index(min(mean_rv))]**2)

RVs_mean_all_means
RVs_mean_all_means_std
peak_to_peak
peak_to_peak_err

max(mbjd)-min(mbjd)
