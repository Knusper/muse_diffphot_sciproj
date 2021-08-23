import sys
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.visualization import simple_norm
from photutils.detection import DAOStarFinder

incube_name = sys.argv[1]
#incube_name = 'LTT_3218wfm-ao-n_2019-05-03.fits'
print(incube_name)

win_width = 40

data = fits.getdata(incube_name, 1)
z_len = data.shape[0]
z_start = int(z_len // 2 - 2.5 * win_width)

cbe = np.asarray([np.nansum(data[z_start + i * win_width:z_start + (i + 1) * win_width,
                                 :,:], axis=0) for i in range(3)])
med_cbe = np.nanmedian(cbe, axis=0)


mea, med, std = sigma_clipped_stats(med_cbe, sigma=3.0)

det_im = med_cbe - med
# mask region around standard star
mask = np.zeros_like(det_im, dtype=bool)
mask[125:175,125:250] = True

outcat_name = incube_name[:-5] + '_starcat.fits'
daofind = DAOStarFinder(threshold=6*std, fwhm=5) # 5px = 1''
sources = daofind(det_im, mask=mask)
sources.write(outcat_name, format='fits', overwrite=True)
print(outcat_name)

outcat_xy_name = incube_name[:-5] + '_starcat_xy.txt'
xy_sources = sources['xcentroid','ycentroid']
xy_sources.write(outcat_xy_name, format='ascii.no_header')
print(outcat_xy_name)

# image for inspection
outim_name = incube_name[:-5] + '_starcat.png'
norm = simple_norm(det_im, percent=90)
fig = plt.figure(1)
ax = plt.subplot(111)
im = ax.imshow(det_im, norm=norm)
sc = ax.scatter(sources['xcentroid'], sources['ycentroid'], marker='x')

for i, (x,y) in enumerate(zip(sources['xcentroid'], sources['ycentroid'])):
    print(str(i+1).zfill(2))
    ax.annotate(str(i+1), (x+3, y+3), color='red')

fig.savefig(outim_name)
print(outim_name)

