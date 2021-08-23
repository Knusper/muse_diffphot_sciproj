import sys
import os
from astroquery.eso import Eso
# python query_standard_and_calib.py 2019-03-09 WFM-NOAO-N outdir
# download standard star dataets given mode / given date - if science
# target was observed during that night

eso_username = 'echerenzops'
date = sys.argv[1]  # '2019-03-09'
mode = sys.argv[2] # 'WFM-NOAO-N'
star = 'LTT 3218'
eso_cache_location = '/homelocal2/eherenz/aq_cache/'

outdir = '/homelocal2/eherenz/muse_diffphot/rawdata_' + date + '_' + mode + '/'
if not os.path.isdir(outdir):
    os.mkdir(outdir)

eso = Eso()

eso.cache_location = eso_cache_location

std_table = eso.query_instrument('muse',
                                 column_filters={'target': star,
                                                 'ins_mode': mode,
                                                 'night': date,
                                                 'dp_type': 'STD'})

if type(std_table) == type(None):
    raise Exception('STD '+star+' not observed on '+date)


obj_table = eso.query_instrument('muse',
                                 column_filters={'ins_mode': mode,
                                                 'night': date,
                                                 'dp_type': 'OBJECT'})

if type(obj_table) == type(None):
    raise Exception('STD '+star+' observed in '+mode+' on '+date+\
                    ', but no OBJECT observed in that mode during that night.')

# select only the first dataset - we are only interested in the
# associated calibration files
dpid = obj_table['DP.ID'][0]

eso.login(eso_username)

eso.retrieve_data(dpid,
                  destination=outdir, with_calib='raw', unzip=False)


