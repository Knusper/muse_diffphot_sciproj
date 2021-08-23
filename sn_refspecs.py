import pickle
import numpy as np
from matplotlib import pyplot as plt
from wavel import wavel
from astropy.io import fits
from astropy.table import Table

wdir = '/media/ehere/work/muse_diffphot/red1/'
filter_curve_fname = './Ibessel.ascii'

ref_specexdict_filename = wdir + 'reference_wfm-ao-n_specexdict.bin' # 2019-03-10 - AM: 1.0145

header =  fits.getheader(wdir + 'LTT_3218wfm-ao-n_2019-03-10.fits', 1)
xax = wavel(header)

filter_curve_tab = Table.read('./Ibessel.ascii', format='ascii')
filter_curve_xax = np.interp(xax,
                             filter_curve_tab['lam'], filter_curve_tab['throu'])

with open(ref_specexdict_filename,'rb') as f1:
    ref_dic = pickle.load(f1)

star_id_max = 25
df=np.dtype(float)
di=np.dtype(int)
sn_refspec_table = Table(names=('id', 1.5, 2.5, 3.5, 4.5, 5.5),
                         dtype=(di, df, df, df, df, df))
for star_id in range(1,star_id_max+1):
    sns_iband = []
    for i,target_radius in enumerate([1.5, 2.5, 3.5, 4.5, 5.5]):
        ref_spec = ref_dic[star_id][i]
        ref_var = ref_dic[star_id, 'var'][i]
        s_iband = np.sum(filter_curve_xax * ref_spec) 
        n_iband = np.sqrt(np.sum(filter_curve_xax**2 * ref_var))
        sn_iband = s_iband / n_iband 
        sns_iband.append(sn_iband)
    
    sn_refspec_table.add_row(tuple([star_id]) + tuple(sns_iband))

sn_refspec_table.write('sn_refspecs_table.ascii',format='ascii.fixed_width')
