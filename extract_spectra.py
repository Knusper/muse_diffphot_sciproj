import pickle
import os
from os import path, popen
import sys
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

def coord_trafo(cat_name, reference_cat_name):
    # uses kd-match from Heyl 2013
    # http://ubc-astrophysics.github.io/kd-match/
    kd_match_path = '~/compiled_projects/kd-match/'  # adjust 
    command_name = kd_match_path + 'quad_kd -p 5 ' + \
                              reference_cat_name + ' ' + cat_name + ' | tail -1'
    print(command_name)
    coord_trafo_output = os.popen(command_name)
    coord_trafo_output = coord_trafo_output.read()
    print(coord_trafo_output)
    coord_trafo_output = coord_trafo_output.split('\n')[0]
    a,b,c,d,e,f = tuple(coord_trafo_output.split()[1:])
    a,b,c,d,e,f = float(a), float(b), float(c), float(d), float(e), float(f)

    ref_cat = Table.read(reference_cat_name, format='ascii.no_header', names=['x','y'])
    x1 = ref_cat['x']
    y1 = ref_cat['y']

    x2 = a * x1 + b * y1 + c
    y2 = d * x1 + e * y1 + f

    return x2, y2
    
wdir = '/media/ehere/work/muse_diffphot/red1/'
pix_rad=[1.5, 2.5, 3.5, 4.5, 5.5]

ref_cat_name = wdir + 'reference_starcat_xy.txt'  # -> LTT_3218wfm-ao-n_2019-03-10_starcat_xy.txt
cube_name = sys.argv[1]

print(cube_name)
star_obs_name = path.basename(cube_name)[:-5] #'LTT_3218wfm-ao-n_2019-03-10'
star_cat_name = wdir + star_obs_name + '_starcat_xy.txt'

try:
    x2, y2 = coord_trafo(star_cat_name, ref_cat_name)
except:
    sys.exit(0)

#cube_name = wdir + 'LTT_3218wfm-ao-n_2019-03-10.fits'

hdu = fits.open(cube_name, memmap=False)
cube = hdu['DATA'].data
varcube = hdu['STAT'].data


star_cat = Table.read(star_cat_name, format='ascii.no_header', names=['x','y'])
ref_stars_ids = np.arange(start=1, stop=len(x2)+1)

spec_dict = {'pix_rad': pix_rad, 'x2': x2, 'y2': y2}

for x2_i, y2_i, ref_star_id in zip(x2, y2, ref_stars_ids):
    print(ref_star_id)
    specs = extract_specs(x2_i, y2_i, cube,
                          pix_rad=pix_rad)
    var_specs = extract_specs(x2_i, y2_i, varcube,
                             pix_rad=pix_rad)

    spec_dict[ref_star_id] = specs
    spec_dict[(ref_star_id,'var')] = var_specs


outfile_name = wdir + star_obs_name + '_specexdict.bin'
outfile = open(outfile_name, 'wb')
pickle.dump(spec_dict, outfile)
outfile.close()
print(outfile_name)
