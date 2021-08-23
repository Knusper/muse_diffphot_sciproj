import sys
import warnings
import datetime
import dateutil.parser as dparser
from glob import glob
import re
import numpy as np
import pickle
from matplotlib import pyplot as plt
from astropy.convolution import Gaussian1DKernel, convolve
from astropy.io import fits
from astropy.table import Table
from wavel import wavel
from brokenaxes import brokenaxes  # https://github.com/bendichter/brokenaxes
warnings.filterwarnings('ignore')


wdir = '/media/ehere/work/muse_diffphot/red1/'  # working directroy


#star_ids = [22,4]  # which stars to use
star_ids = [int(sys.argv[1]), int(sys.argv[2])]  # supply star IDs as
                                                 # command line
                                                 # arguments
target_radius = 4.5 # aperture radius (spectra extracted for 
                    # [1.5,2.5, 3.5, 4.5, 5.5])
lim = 0.2 # range of plot (value outside of plot will be shown as arrows)


header =  fits.getheader(wdir + 'LTT_3218wfm-ao-n_2019-03-10.fits', 1)
xax = wavel(header)
filter_curve_tab = Table.read('./Ibessel.ascii', format='ascii')

filter_curve_xax = np.interp(xax,
                             filter_curve_tab['lam'], filter_curve_tab['throu'])

modes = ['wfm-noao-n', 'wfm-ao-n', 'wfm-noao-e', 'wfm-ao-e']  # now different modes
colors = ['blue', 'green', 'cyan', 'y']

fig = plt.figure(figsize=(0.8*15,0.8*5))

# adjust here where the axes should be broken - depends on data
bax = brokenaxes(
    xlims=((datetime.datetime(2019, 1, 1),
            datetime.datetime(2019, 6, 1),),
           (datetime.datetime(2019, 11, 1),
            datetime.datetime(2019, 12, 31),)
    ))

diff_dict = {}
dates_dict = {}
for mode in modes:
    # filenames of files storing extracted spectra are in the form
    # LTT_3218wfm-ao-n_2019-01-09_specexdict.bin
    glob_list = glob(wdir + 'LTT_3218' + mode + '_' + '????-??-??' + '_specexdict.bin')
    date_strs = [re.search(r'\d{4}-\d{2}-\d{2}', s)[0]
                 for s in glob_list]
    dates_dict[mode] = np.asarray([dparser.parse(s) for s in date_strs])

    diffs = []
    for ana_specexdict_filename in glob_list:
        with open(ana_specexdict_filename,'rb') as f1:
            ref_dic = pickle.load(f1)

        pix_rad = ref_dic['pix_rad']
        target_index = pix_rad.index(target_radius)
            
        ref_specs = np.asarray([ref_dic[star_id][target_index] * filter_curve_xax
                                for star_id in star_ids])

        ref_sum = [np.sum(rs) for rs in ref_specs]
        diff = ref_sum[0] / ref_sum[1]
        diffs.append(diff)

    diff_dict[mode] = np.asarray(diffs)

# mean over all measurements (all modes)    
diffs_mean = np.mean(np.concatenate([diff_dict[mode] for mode in modes]))

diffs_final_d = []
for mode, color in zip(modes, colors):
    dates = dates_dict[mode]
    diffs_final = (diff_dict[mode] - diffs_mean)/diffs_mean

    diffs_final_d.extend(diffs_final)
    
    outlier_sel = np.abs(diffs_final) >= lim

    bax.plot_date(dates[~outlier_sel], diffs_final[~outlier_sel], xdate=True,
                  color=color, label=mode)

    # plot outliers
    bax.plot_date(dates[diffs_final >= lim],
                  lim * 0.95 * np.ones_like(dates[diffs_final >= lim]),
                  fmt='^', color=color, markersize=9)
    bax.plot_date(dates[diffs_final <= -lim],
                  -lim * 0.95 * np.ones_like(dates[diffs_final <= -lim]),
                  fmt='v', color=color, markersize=9)

diffs_final_d = np.asarray(diffs_final_d)
diffs_final_d.sort()

bax.tick_params(axis='both', which='major', labelsize=10)
bax.legend()

bax.set_ylim((-lim, lim))

fig.autofmt_xdate()


plt.title('Stars: ' + str(star_ids[0]) + ' & ' + str(star_ids[1]) + \
          ' ; std = ' + str(np.round(np.std(diffs_final_d[1:-1]),2)))


bax.set_xlabel(r'Date', fontsize=15)
bax.set_ylabel(r'$(F_1 / F_2) / < F_1 / F_2 >_t - 1$', fontsize=15, labelpad=50)

#plt.tight_layout()
plt.subplots_adjust(top=0.925, bottom=0.240)
[x.remove() for x in bax.diag_handles]
bax.draw_diags()

#plt.show()

outfigname = 'spec_spat_diff_s' + \
            str(star_ids[0]).zfill(2) + 's' + \
            str(star_ids[1]).zfill(2) + '.png'
plt.savefig(outfigname)
print(outfigname)


