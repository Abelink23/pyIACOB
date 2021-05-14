import sys; sys.path.append('../')

from spec import *


lwl = 4510; rwl = 4600

stars = [
'HD30614',
'HD24431',
'HD216898',
'HD195592',
'HD189957',
'HD36512',
'HD37128',
'HD48434',
'HD149438',
'HD41117',
'HD35468',
'HD144218',
'HD14134',
'HD36212',
'HD161573',
'HD183143',
'HD23408',
'HD90994']

nrows = int(len(stars)/3)
if len(stars) % 3 != 0.0: nrows += 1
if len(stars) < 3: ncols = len(stars)
else: ncols = 3

for star,nplot in zip(stars,range(len(stars))):

    star = spec(star,SNR='best'); star.spc()

    star.cosmic()
    star.degrade(resol=50)

    mask = (star.wave > lwl) & (star.wave < rwl)

    plt.subplot(nrows,ncols,nplot+1)

    plt.plot(star.wave[mask],star.flux[mask],linewidth=.3)
    plt.tick_params(direction='in',top='on')
    plt.ylim(ymin=0.9)
    plt.title(star.name_star+' '+star.SpC,size=9)
    plt.subplots_adjust(hspace=0,wspace=0)

plt.subplots_adjust(hspace=0,wspace=0)

plt.tight_layout(); plt.show(block=True)
