import itertools
import dateutil.parser as dparser
from glob import glob
import re
import numpy as np
import pickle
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.table import Table
from tqdm import tqdm
from wavel import wavel


wdir = '/media/ehere/work/muse_diffphot/red1/'  # working directro
target_radius = 4.5 # aperture radius (spectra extracted for 
                    # [1.5,2.5, 3.5, 4.5, 5.5])

# IBessel filter curve on MUSE grid                    
header =  fits.getheader(wdir + 'LTT_3218wfm-ao-n_2019-03-10.fits', 1)
xax = wavel(header)
filter_curve_tab = Table.read('./Ibessel.ascii', format='ascii')
filter_curve_xax = np.interp(xax,
                             filter_curve_tab['lam'],
                             filter_curve_tab['throu'])

# sn refspec table
snr_tab = Table.read('sn_refspecs_table.ascii',
                            format='ascii.fixed_width')
# files
glob_list = glob(wdir + 'LTT_3218*_????-??-??_specexdict.bin')
# all 300 combinations of stellar
pair_list = list(itertools.chain(*[[(i,j)
                                    for j in range(i+1,26)]
                                   for i in range(1,26)])) 
    
snrsh = []
snrsg = []
stds = []
for pair in tqdm(pair_list):
    # calc avg sn
    id1, id2 = pair[0], pair[1]
    if id1 == 25 or id2 == 25:
        continue #  star to close to the edge
    snr1 = snr_tab[str(target_radius)][snr_tab['id'] == id1]
    snr2 = snr_tab[str(target_radius)][snr_tab['id'] == id2]
    snr_hmean = 2/((1/snr1) + 1/snr2)
    snr_gmean = (snr1 * snr2)**(1/2)
    # calc std
    diffs = []
    for ana_specexdict_filename in glob_list:
        with open(ana_specexdict_filename,'rb') as f1:
            ref_dic = pickle.load(f1)

        pix_rad = ref_dic['pix_rad']
        target_index = pix_rad.index(target_radius)
            
        ref_specs = np.asarray([ref_dic[star_id][target_index] * filter_curve_xax
                                for star_id in [id1, id2]])
        ref_sum = [np.sum(rs) for rs in ref_specs]
        diff = ref_sum[0] / ref_sum[1]
        if not diff != diff:
            diffs.append(diff)

    diffs = np.asarray(diffs)
    mean_diff = np.mean(diffs)
    diffs_final = (diffs - mean_diff) / mean_diff

    # crop the 2 worst outliers in std
    diffs_final = np.sort(diffs_final)
    std = np.std(diffs_final[1:-1])

    if not std != std:
        snrsh.append(snr_hmean)
        snrsg.append(snr_gmean)
        stds.append(std)

plt.scatter(snrsg, stds)
plt.xscale('log')
plt.xlim((100,10000))
plt.xlabel('SNR_gmean', fontsize=15)
plt.ylabel('std. dev.', fontsize=15)

plt.tight_layout()

plt.show()
