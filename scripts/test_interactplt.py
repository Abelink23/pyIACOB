import sys
sys.path.append('../')

from db import *

import matplotlib; matplotlib.use('QT5Agg')
import matplotlib.pyplot as plt

from tools_plt import *

def on_lclick(event):
    if event.inaxes and event.button is MouseButton.LEFT:
        clickX = event.xdata
        clickY = event.ydata
        print('Data coords %f %f' % (event.xdata, event.ydata))

def on_rclick(event):
    if event.inaxes and event.button is MouseButton.RIGHT:
        print('Disconnecting callback')
        plt.disconnect(binding_id)

#def on_key(event):
#    print('you pressed', event.key, event.xdata, event.ydata)
#cid = fig.canvas.mpl_connect('key_press_event', on_key)

plt.close('all')
fig, ax = plt.subplots(tight_layout=True)

table = findtable('O9BAs_RVEWFWs.fits')
x = table['EWSiIII1']; y = table['EWSiIII2']; names = table['Name']

test = 1
if test == 1:

    pts = ax.scatter(x,y,s=10,c='b')
    selector = SelectFromCollection(ax, pts)

    def accept(event):
        if event.key == "enter":
            print("Selected points:")
            print(selector.xys[selector.ind])
            selector.disconnect()
            ax.set_title("")
            fig.canvas.draw()
            accept.table_f = table[selector.mask]

    fig.canvas.mpl_connect("key_press_event", accept)
    ax.set_title("Press enter to accept selected points.")

    plt.show(block=False)

elif test == 2:

    ax.scatter(x,y,s=10,c='b')

    binding_id = plt.connect('button_press_event', on_lclick)

    plt.connect('button_press_event', on_rclick)
    #binding_id = fig.canvas.mpl_connect('button_press_event', on_lclick)
    #fig.canvas.mpl_connect('button_press_event', on_rclick)

    plt.show(block=False)

accept.table_f
