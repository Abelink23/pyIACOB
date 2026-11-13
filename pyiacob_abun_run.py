"""
@author: A. Casasbuenas Corral

GOAL:
    Quantitative spectroscopic analysis of B-type stars, mainly within the IACOB database. 
"""

###############################################################################################
#%% 
###############################################################################################

print("----------------------------------------------------------------------------------------------------------------")
print("FIRST STEP: FIT THE SILICON LINES OF A B-TYPE STAR INTO GAUSSIANS TO GET THE EQUIVALENT WIDTH MEASUREMENTS.")
print("----------------------------------------------------------------------------------------------------------------")

stop = input("Press return to continue.")

from measure import *

# MEDIR LAS EW DE LAS ESTRELLAS 
el_input = "Si"
lines_list = "%s_lines_FASTWIND.lst"%el_input

star_db = "WS_parameters.fits"


rvlines = 'rv_Bs.lst'
table_measure = "WS_measure_%s.fits"%el_input
orig = "txt"  # if .csv

measure(lines_list, star_db, output_table=table_measure, 
            orig = orig,col_line_names='label', rv_lines=rvlines, 
            show_plot=True) #, redo = "HD36822" if you want to repeat it

###############################################################################################
#%% 
###############################################################################################


print("----------------------------------------------------------------------------------------------------------------")
print("SECOND STEP: QUANTITATIVE SPECTROSCOPIC ANALYSIS OF A B-TYPE STAR: TEFF, LOGG, LOGQ, AB(SI), MICROTURBULENCE")
print("----------------------------------------------------------------------------------------------------------------")

stop = input("Press return to continue.")

from pyiacob_abun_functions import * 

savepath=maindir+'pyiacobabun_results/' # DEFINE HERE YOUR PATH FOR STORING THE RESULTS

el_input='Si'

fastwinddir = maindir + "FASTWIND/"
FW_EW_measure = "measure_FASTWIND_iac_bmsgrid_HHeCNOSimu_Z100_He10_240131_xdr_%s.fits"%el_input

# star_db = "BMS_STD_parameters.fits"
star_db = "WS_parameters.fits"

star_measure = "WS_measure_%s.fits"%el_input

try:
    os.mkdir(savepath)
except:
    print("Directory already exists.")

analysis_allstars = 'analysis_WS.csv'

linelist='%s_lines_FASTWIND.lst'%(el_input.capitalize())
finallinelist='lines_finalplot.lst'
lines='initial_lines_%s.lst'%(el_input.capitalize())

linestab = '%s_lines_star.csv' %(el_input.capitalize())
if findtable(linestab, path = savepath, format = "csv") is None: 
    lines_full = findlines(linelist)
    t_lines = Table([ np.array(lines_full[1]), np.array(lines_full[0]), ["green"]*len(lines_full[1]) ], 
                    names = ["label", "WL", "flag"]
                    #, dtype= ["S20", "float64", "S10"]
                    )
    t_lines.write(savepath + linestab, format = "ascii.csv", overwrite=True, delimiter=(','))     

db_STD = findtable(star_db)

# ANÁLISIS DE UNA ESTRELLA: 

star = "HD36822" # for example

ref_file = db_STD["Ref_file"][db_STD["ID"] == "HD36822"][0]

if '.fits' in ref_file:
    orig = 'IACOB'
else:
    orig = 'txt'     

t=parameters(ref_file,el_input,star_db,star_measure,finallinelist,linestab,lines,
            FW_EW_measure,analysis_allstars,savepath, auto="n", show = "y") #parámetros + ab(Si)

make_summary_pdf(ref_file,savepath)

print("DONE!")

###############################################################################################

#%%
###############################################################################################
print("----------------------------------------------------------------------------------------------------------------")
print("THIRD STEP: FIT THE OXYGEN LINES OF A B-TYPE STAR INTO GAUSSIANS TO GET THE EQUIVALENT WIDTH MEASUREMENTS.")
print("----------------------------------------------------------------------------------------------------------------")

stop = input("Press return to continue.")

from measure import *

# MEDIR LAS EW DE LAS ESTRELLAS 
el_input = "O"
lines_list = "%s_lines_FASTWIND.lst"%el_input

star_db = "WS_parameters.fits"

rvlines = 'rv_Bs.lst'
table_measure = "WS_measure_%s.fits"%el_input
orig = "txt"  # if .csv

measure(lines_list, star_db, output_table=table_measure, 
            orig = orig,col_line_names='label', rv_lines=rvlines, 
            show_plot=True)

# ##############################################################################################
# %%
# ##############################################################################################

from pyiacob_abun_functions import * 

print("----------------------------------------------------------------------------------------------------------------")
print("FOURTH STEP: OXYGEN SURFACE ABUNDANCE DETERMINATION FOR A B-TYPE STAR GIVEN THE PREVIOUSLY DETERMINED PARAMETERS")
print("----------------------------------------------------------------------------------------------------------------")

stop = input("Press return to continue.")


# ABUNDANCIAS DE UN ELEMENTO SABIENDO LOS PARÁMETROS: 

el_input = "O"

savepath=maindir+'pyiacobabun_results/' # DEFINE HERE YOUR PATH FOR STORING THE RESULTS

star_db = "WS_parameters.fits"
star_measure = "WS_measure_%s.fits"%el_input
linelist='%s_lines_FASTWIND.lst'%el_input.capitalize()  
FW_EW_measure = "measure_FASTWIND_iac_bmsgrid_HHeCNOSimu_Z100_He10_240131_xdr_%s.fits"%el_input
analysis_allstars = 'analysis_WS.csv'

linestab = '%s_lines_star.csv' %(el_input.capitalize()) 

if findtable(linestab, path = savepath, format = "csv") is None: 
    lines_full = findlines(linelist)
    t_lines = Table([ np.array(lines_full[1]), np.array(lines_full[0]), ["green"]*len(lines_full[1]) ], 
                    names = ["label", "WL", "flag"]
                    #, dtype= ["S20", "float64", "S10"]
                    )
    t_lines["flag"] = ["red" if "b" in t_lines["label"][i]
                       else "green" for i in range(len(t_lines))] # marks a blend in the line list
    t_lines.write(savepath + linestab, format = "ascii.csv", overwrite=True, delimiter=(','))     

# ANÁLISIS DE UNA ESTRELLA: 
db_STD = findtable(star_db)
star = "HD36822" # for example
ref_file = db_STD["Ref_file"][db_STD["ID"] == "HD36822"][0]

abundances_element(ref_file,el_input,star_db,star_measure,linestab,
                FW_EW_measure,analysis_allstars,
                savepath,mult_use="all",mult_highlight=None,interpolate=False,
                microreduced="fixed", unc = "n", auto="n")

wave, flux = plot_best_model_element(ref_file,el_input,star_db,star_measure,
                linestab,FW_EW_measure,analysis_allstars,savepath)

###############################################################################################
# %%