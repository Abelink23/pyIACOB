import os

import numpy as np
import matplotlib.pyplot as plt

import lightkurve as lk
import scipy.signal as sig

from astropy import units as u
from astropy.table import Table


SB1 = Table.read('/Users/abelink/MEGA/PhD/tables/peridiograms/HD199579', format = 'ascii')
SB1.rename_column('col6','t')
SB1.rename_column('col7','vr')

# https://radvel.readthedocs.io/en/latest/quickstartcli.html#example-fit
# https://docs.scipy.org/doc/scipy/reference/signal.html

def liku(table):
    #lc = lk.LightCurve(time = SB1['time'][SB1['time'] > 2.4558e6], flux = SB1['vrad'][SB1['time'] > 2.4558e6])
    lc = lk.LightCurve(time = table['t'], flux = table['vr'])

    lc.scatter()
    #pg = lc.to_periodogram(oversample_factor=1)
    pg = lc.to_periodogram(minimum_period=1*u.day, maximum_period=100*u.day, oversample_factor=10)
    pg.plot()
    print(pg.period_at_max_power)
    #[pg.period[i] for i in range(pg) if pg.power[i] > 100]

    lc.fold(pg.period_at_max_power).scatter()


def lombsc(table):
    periods = 10**np.linspace(-3,3,10000)
    omega = 2*np.pi/periods

    power = sig.lombscargle(SB1['t'],SB1['vr'],omega)
    max_pers = np.argsort(power)

    print('Best f =', round(omega[max_pers[-1]],3))
    print('Best w =', round(2*np.pi/omega[max_pers[-1]],3))

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.semilogx(periods,power)
    ax.set_xlabel('Period (days)',fontsize=15)
    ax.set_ylabel('Power',fontsize=20)

    plt.show(); return None

liku(SB1)
lombsc(SB1)
