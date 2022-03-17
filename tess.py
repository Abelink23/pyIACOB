from db import *

from astroquery.mast import Observations

from lightkurve import search_targetpixelfile, search_lightcurve

# https://docs.lightkurve.org/
# https://docs.lightkurve.org/tutorials/
# https://docs.lightkurve.org/reference/api/lightkurve.KeplerTargetPixelFile.html

class LC():

    def __init__(self, ID, kic=None):

        self.id_star = ID
        if kic is not None:
            self.kic = kic


    def query_lc(self, mission=[], quarter=None, sector=None):

        if hasattr(self, 'kic'):
            lc = search_lightcurve(self.kic, mission=mission, quarter=quarter)
        else:
            lc = search_lightcurve(self.id_star, mission=mission, quarter=quarter, sector=sector)

        print(lc)
        self.lc = lc.download()


    def query_pixfile(self, mission=[], quarter=None, sector=None):

        if hasattr(self, 'kic'):
            tpf = search_targetpixelfile(self.kic, mission=mission, quarter=quarter)
        else:
            tpf = search_targetpixelfile(self.id_star, mission=mission, quarter=quarter, sector=[])

        print(tpf)
        self.tpf = tpf.download()


    def change_aperture(self):

        '''
        Function to visually change the tpf mask.

        Parameters
        ----------
        tpf : lightkurve.targetpixelfile.KeplerTargetPixelFile
            Enter the tpf.

        Returns
        -------
        New mask for the fpf.
        '''

        if hasattr(self, 'tpf'):

            new_mask = self.tpf.pipeline_mask; change = 'y'
            while change == 'y':

                print('Showing original mask...')
                self.tpf.plot(aperture_mask=new_mask)
                print([[1 if i == True else 0 for i in j] for j in new_mask.tolist()])

                change = input('Do you want to change it? [n/y]: ')
                if change == 'y':

                    raw,col,val = input('Use raw,column,0/1: ').split(',')
                    if val == '1':
                        val = True
                    elif val == '0':
                        val = False

                    new_mask[int(raw)][int(col)] = val

            print('Showing final mask...')
            self.tpf.plot(aperture_mask=new_mask)

            self.tpf.new_mask = new_mask


    def tpf_to_lc(self, mask='default'):

        if hasattr(self, 'tpf'):
            if mask == 'default':
                mask = self.tpf.pipeline_mask
            self.lc = self.tpf.to_lightcurve(aperture_mask=mask)
        else:
            print('No tpf file link to the class...')


    def get_mag(self):

        if hasattr(self.lc, 'remove_nans'):
            flux = self.lc.remove_nans().flux.value
        else:
            flux = self.flux.value

        mag = -2.5 * np.log10(flux)
        mag -= np.average(mag)
        self.lc.magnitude = mag


    def export(self, output_path=maindir+'lightcurves/', format='ascii'):

        np.savetxt(output_path+"%s.txt" % self.id_star, np.array([self.lc.time.value,self.lc.flux.value]).T, delimiter=",")


search_lightcurve('KIC 10963065')
search_lightcurve('HD 42608')

search_lightcurve('HD 42608', quarter=None, sector=None, mission=[])

search_lightcurve(self.kic)
help(search_lightcurve)

test = LC(ID='HD 42608')
test.query_lc()
test.lc.target_name
test.query_lc(mission='TESS', sector=6)

pg = test.lc.normalize(unit='ppm').to_periodogram()
pg.plot()
test.tpf_to_lc()
test.get_mag()
test.lc.magnitude

test.lc.PDCSAP_FLUX.remove_nans()
test.export()


lcf = search_lightcurvefile("kic 6922244", quarter=4).download()
tpf = search_targetpixelfile('KIC 6922244', mission='TESS', quarter=4).download()
