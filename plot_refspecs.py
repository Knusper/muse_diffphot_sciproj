import pickle
from matplotlib import pyplot as plt
from astropy.convolution import Gaussian1DKernel, convolve
from astropy.io import fits

gk_width = 3.
gk = Gaussian1DKernel(gk_width)

wdir = '/media/ehere/work/muse_diffphot/red1/'
star_id_max = 25

ref_specexdict_filename = wdir + 'reference_wfm-ao-n_specexdict.bin' # 2019-03-10 - AM: 1.0145

with open(ref_specexdict_filename,'rb') as f1:
    ref_dic = pickle.load(f1)


for star_id in range(1,star_id_max+1):
    for i,target_radius in enumerate([1.5, 2.5, 3.5, 4.5, 5.5]):
        ref_spec = ref_dic[star_id][i]
        ref_spec_conv = convolve(ref_spec, gk)
        line = plt.plot(ref_spec_conv, label='rad = '+str(target_radius))
        plt.plot(ref_spec, alpha=0.3, color=line[0].get_color())

    plt.title('star_id: '+str(star_id))    
    plt.legend()
    outfig_name = 'plot_refspec_'+str(star_id).zfill(2)+'.png'
    plt.savefig(outfig_name)
    print(outfig_name)
    plt.close()
