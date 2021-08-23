import pickle
from os import path
from glob import glob
from astropy.io import fits
from astropy.table import Table
import numpy as np

def extract_specs(x, y, cube, pix_rad=[1.5, 2.5, 3.5, 4.5, 5.5]):
    Y,X = np.indices(cube.shape[1:])
    rad = np.sqrt((x - X)**2 + (y - Y)**2)
    apps = [(rad <= pix_rad_i).astype('float')
            for pix_rad_i in pix_rad]
    specs = [np.nansum(np.nansum(app * cube,
                                 axis=1),
                    axis=1)
             for app in apps]
    return specs

wdir = '/media/ehere/work/muse_diffphot/red1/'
pix_rad=[1.5, 2.5, 3.5, 4.5, 5.5]

# mode = 'WFM-AO-N'
# cube_names = glob(wdir + 'LTT_3218' + mode.lower() + '_????-??-??.fits')
# cube_names.sort()

# REMEMBA reference_starcat_xy.txt -> LTT_3218wfm-ao-n_2019-03-10_starcat_xy.txt
star_obs_name = 'LTT_3218wfm-ao-n_2019-03-10'
ref_cat_name = wdir + star_obs_name + '_starcat_xy.txt'
cube_name = wdir + 'LTT_3218wfm-ao-n_2019-03-10.fits'
outfile_name = wdir + 'reference_wfm-ao-n_specexdict.bin'

hdu = fits.open(cube_name)
cube = hdu['DATA'].data
var_cube = hdu['STAT'].data

ref_cat = Table.read(ref_cat_name, format='ascii.no_header', names=['x','y'])
x1 = ref_cat['x']
y1 = ref_cat['y']
ref_stars_ids = np.arange(start=1, stop=len(x1)+1)
    
spec_dict = {'pix_rad': pix_rad}

for x1_i, y1_i, ref_star_id in zip(x1, y1, ref_stars_ids):
    print(ref_star_id)
    specs = extract_specs(x1_i, y1_i, cube,
                          pix_rad=pix_rad)
    var_specs = extract_specs(x1_i, y1_i, var_cube,
                             pix_rad=pix_rad)
    
    spec_dict[ref_star_id] = specs
    spec_dict[(ref_star_id,'var')] = var_specs

outfile = open(outfile_name, 'wb')
pickle.dump(spec_dict, outfile)
outfile.close()
