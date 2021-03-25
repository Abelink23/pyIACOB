'''=============================================================================
Script to make different plots using fits tables
============================================================================='''
from db import *
from tools import *

from matplotlib.patches import Rectangle,Ellipse

ftable =
id =
x = ; y =
e_x = 0; e_y = 0,
coords = None
dist = None
distunit = None
ratio = None
outlicol = None
sig = 2
pm = None
names = False
title = None
grid = False
xlim = None; ylim = None
invert = None

'''===================== Find and load the table ========================'''
with fits.open(findtable(ftable),mode='readonly') as hdu_list:
    data = hdu_list[1].data


'''=============================== ID name =============================='''
try: data[id]
except: id = input('Error: Enter a valid name for the ID names: ')


'''============================= Get x and y ============================'''
try: data[x]
except: x = input('Error: Enter a valid name for the x axis to plot: ')

if type(e_x) == str:
    try: data[e_x]
    except: e_x = input('Error: Enter a valid name for the x axis error: ')

try: data[y]
except: y = input('Error: Enter a valid name for the y axis to plot: ')

if type(e_y) == str:
    try: data[e_y]
    except: e_y = input('Error: Enter a valid name for the y axis error: ')


'''============================= Guess RADEC ============================'''
RADEC = None
if 'RA' in x and 'DEC' in y:
    RADEC = x + ',' + y; RA = RADEC.split(',')[0]; DEC = RADEC.split(',')[1]
else:
    RADEC = input("If RA/DEC are used, provide their column names, e.g. 'RA,DEC': ")
    if RADEC != '': RA = RADEC.split(',')[0]; DEC = RADEC.split(',')[1]
    else: RADEC = None


'''========================== RA/DEC to distance ========================'''
units = {'pc' : u.pc, 'kpc' : u.kpc, 'au' : u.au, 'lyr' : u.lyr, 'deg' : u.deg}
if (dist != None or distunit != None) and RADEC != None:

    if dist == None:
        print('Provide the distance as input parameter.')
    if not distunit in units:
        print('The distance unit is not valid.')
    if coords == None:
        print('Provide the central coordinates in deg.')

    RA_0 = float(coords.split()[0]); DEC_0 = float(coords.split()[1])

    d = dist*units[distunit]
    center = SkyCoord(coords,unit=u.deg,distance=d,frame='icrs')
    if ':' in str(data[0][RA]) and ':' in str(data[0][DEC]):
        RA_1 = SkyCoord(data[RA],DEC_0,unit=(u.hourangle,u.deg),distance=d,frame='icrs')
        RADEC_1 = SkyCoord(data[RA],data[DEC],unit=(u.hourangle,u.deg),distance=d,frame='icrs')
    else:
        RA_1 = SkyCoord(data[RA],DEC_0,unit=(u.deg,u.deg),distance=d,frame='icrs')
        RADEC_1 = SkyCoord(data[RA],data[DEC],unit=(u.deg,u.deg),distance=d,frame='icrs')
    DEC_1 = SkyCoord(RA_0,data[DEC], unit=(u.deg,u.deg),distance=d,frame='icrs')

    newxy = input('Provide a column name to replace by the total distance: ')

    if distunit != 'deg':
        data[RA] = center.separation_3d(RA_1).value
        data[DEC] = center.separation_3d(DEC_1).value
        data[newxy] = center.separation_3d(RADEC_1).value
    else:
        data[RA] = center.separation(RA_1).value
        data[DEC] = center.separation(DEC_1).value
        data[newxy] = center.separation(RADEC_1).value


'''===================== SpC and color classification ==================='''
spc_column = input('Enter SpC column name SpT groups. Otherwise hit return: ')

if spc_column != '':
    try: data[spc_column]
    except: spc_column = input('Error: Enter a valid name for the SpC column name: ')
    spc_list = np.asarray([i[0] for i in data[spc_column]])

    stars_O = data[data[spc_column].startswith('O')]
    stars_B = data[data[spc_column].startswith('B')]
    stars_A = data[data[spc_column].startswith('A')]
    stars_M = data[data[spc_column].startswith('M')|\
                   data[spc_column].startswith('F')|\
                   data[spc_column].startswith('G')]
    stars_spt = [stars_O,stars_B,stars_A,stars_M]

else: spc_list = np.asarray(['-']*len(data))


'''=============================== Outliers ============================='''
bads = []
if outlicol != None:
    sclip = [];
    if spc_column != '': spc_classif = input('Remove outliers by spectral type? [y,n] ')
    if spc_column != '' and spc_classif == 'y':
        for group,spt in zip(range(len(stars_spt)),['O','B','A','M']):
            try:
                # This has to be modified...
                sclip = outliers(group[outlicol],siglo=sig,sigup=sig)
                stars_spt[group] = sclip[0]; bads.append(sclip[1][id])
                print('Clipped group %s has a mean value for %s column of: ' %(spt,outlicol))
                print(str(np.mean(sclip[0]['RV_best'])) + ' +- ' + str(np.std(sclip[0]['RV_best'])))
            except: print('Outliers removal failed, check outliers column name. ')

    else:
        try: sclip.append(outliers(data[outlicol]))
        except: print('Outliers removal failed, check outliers column name. ')
    bads = [j for i in bads for j in i]


'''================================= Plot ==============================='''
fig, ax = plt.subplots(figsize=ratio,tight_layout=True)
colorcode = {'O':'purple', 'B':'b', 'A':'c', 'F':'r', 'K':'r', 'M':'r', '-':'g'}

for source,cc_value in zip(data,spc_list):

    #if xlim != None and type(xlim) == tuple:
    #    if not (xlim[0] < source[x] < xlim[1]): continue
    #if ylim != None and type(ylim) == tuple:
    #    if not (ylim[0] < source[y] < ylim[1]): continue


    '''====================== Manual selection rules ===================='''
    #if source['RV_all'] == 0: continue

    '''======================== Add proper motions ======================'''
    # https://matplotlib.org/api/_as_gen/matplotlib.pyplot.quiver.html
    if pm != None:
        try: pmRA,pmDEC = pm.split(',')
        except: print('Wrong input for proper motion.')

        pm_x = source[pmRA]+0.483; pm_y = source[pmDEC]+0.9058
        if invert != None and 'x' in invert: pm_x = -pm_x
        if invert != None and 'y' in invert: pm_y = -pm_y
        ax.quiver(source[x],source[y],pm_x,pm_y,
                    color=colorcode[cc_value], width = .02, scale=2, units='xy')

    try:
        if source[e_x] < 1: e_x_p = 0
        else: e_x_p = source[e_x]
    except: e_x_p = e_x
    try:
        if source[e_y] < 1: e_y_p = 0
        else: e_y_p = source[e_y]
    except: e_y_p = e_y

    if source[id] in bads:
        ax.errorbar(source[x],source[y],xerr=e_x_p,yerr=e_y_p,marker='o',
        elinewidth=.5,markersize=6,mfc='none',color=colorcode[cc_value])

    else: ax.errorbar(source[x],source[y],xerr=e_x_p,yerr=e_y_p,fmt='*',
          elinewidth=.5,markersize=6,color=colorcode[cc_value])

    if names == True: ax.text(source[x]+.5,source[y]+.5,source[id],fontsize=6)


'''============= Plot rectangles and mean with dashed lines ============='''
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Verdana'
plt.rcParams['mathtext.default'] = 'regular'


'''======================== Manual plotting rules ======================='''


'''=========================== Plot parameters =========================='''
ax.set_title(title)
if grid == True: ax.grid()
if invert != None and 'x' in invert: ax.invert_xaxis()
if invert != None and 'y' in invert: ax.invert_yaxis()
if xlim != None and type(xlim) == tuple: ax.set_xlim(xlim)
if ylim != None and type(ylim) == tuple: ax.set_ylim(ylim)
ax.set_xlabel(x); ax.set_ylabel(y)
ax.tick_params(labeltop=True,labelright=True,direction='in',top='on',right='on')
ax.figure.subplots_adjust(top=.9,bottom=.1,right=.9,left=.1)

#ax.yaxis.set_label_position("right")
#ax.xaxis.set_label_position("top")

completeName = os.path.join(maindir+'tmp_plots/plotable.jpg')
ax.figure.savefig(completeName,format='jpg',dpi=300)

plt.show(block=False)

print(bads)
