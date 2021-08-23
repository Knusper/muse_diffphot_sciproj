import sys
import warnings
from datetime import datetime, timedelta
from calendar import isleap

from tqdm import tqdm  

from astropy.config import set_temp_cache
from astropy.table import Table
from astroquery.eso import Eso

from astroquery.exceptions import NoResultsWarning
warnings.filterwarnings("ignore", category=NoResultsWarning)

''' 
Query ESO Archive for standard star exposures of a defined standard.
Check whether at a given night where the standard was observed also a
science exposure was taken with the given setting.  Only for those
cases a simple archive query can deliver the standard star observation
and the associated calibration data.
'''

eso_username = 'echerenzops'
year = sys.argv[1]  # e.g. 2019
star = 'LTT 3218'

eso = Eso()

num_days = 366 if isleap(int(year)) else 365
date_list = [(datetime.strptime(year, '%Y') + \
              timedelta(days=i)).strftime('%Y-%m-%d')
             for i in range(num_days)] #
wfm_noao_n = []
wfm_noao_e = []
wfm_ao_n = []
wfm_ao_e = []
date_out_list = []

for date in tqdm(date_list):
    std_table_noao_n = eso.query_instrument('muse',
                                column_filters={'target': star,
                                                'ins_mode': 'WFM-NOAO-N',
                                                'night': date,
                                                'dp_type': 'STD'})
    std_table_noao_e = eso.query_instrument('muse',
                                column_filters={'target': star,
                                                'ins_mode': 'WFM-NOAO-E',
                                                'night': date,
                                                'dp_type': 'STD'})
    std_table_ao_n = eso.query_instrument('muse',
                                column_filters={'target': star,
                                                'ins_mode': 'WFM-AO-N',
                                                'night': date,
                                                'dp_type': 'STD'})
    std_table_ao_e = eso.query_instrument('muse',
                                column_filters={'target': star,
                                                'ins_mode': 'WFM-AO-E',
                                                'night': date,
                                                'dp_type': 'STD'})
    
    if type(std_table_noao_n) != type(None):
        obj_table = eso.query_instrument('muse',
                                         column_filters={'ins_mode': 'WFM-NOAO-N',
                                                         'night': date,
                                                         'dp_type': 'OBJECT'})

        if type(obj_table) == type(None):
            wfm_noao_n.append('NOSCI')
        else:
            wfm_noao_n.append('SCI')
    else:
        wfm_noao_n.append('')
            
    if type(std_table_noao_e) != type(None):
        obj_table = eso.query_instrument('muse',
                                         column_filters={'ins_mode': 'WFM-NOAO-E',
                                                         'night': date,
                                                         'dp_type': 'OBJECT'})
        if type(obj_table) == type(None):
            wfm_noao_e.append('NOSCI')
        else:
            wfm_noao_e.append('SCI')
    else:
        wfm_noao_e.append('')

    if type(std_table_ao_n) != type(None):
        obj_table = eso.query_instrument('muse',
                                         column_filters={'ins_mode': 'WFM-AO-N',
                                                         'night': date,
                                                         'dp_type': 'OBJECT'})
        if type(obj_table) == type(None):
            wfm_ao_n.append('NOSCI')
        else:
            wfm_ao_n.append('SCI')
    else:
        wfm_ao_n.append('')

    if type(std_table_ao_e) != type(None):
        obj_table = eso.query_instrument('muse',
                                         column_filters={'ins_mode': 'WFM-AO-E',
                                                         'night': date,
                                                         'dp_type': 'OBJECT'})
        if type(obj_table) == type(None):
            wfm_ao_e.append('NOSCI')
        else:
            wfm_ao_e.append('SCI')
    else:
        wfm_ao_e.append('')

    if wfm_noao_n[-1] == wfm_noao_e[-1] == wfm_ao_n[-1] == wfm_ao_e[-1] == '':
        del wfm_noao_n[-1]
        del wfm_noao_e[-1]
        del wfm_ao_n[-1]
        del wfm_ao_e[-1]
    else:
        date_out_list.append(date)
        
out_table = Table([date_out_list,
                   wfm_noao_n, wfm_noao_e, wfm_ao_n, wfm_ao_e],
                  names=['DATE',
                         'WFM_NOAO_N', 'WFM_NOAO_E', 'WFM_AO_N', 'WFM_AO_E'])


out_table.write('standard_overview_table_' + str(year) + '.ascii',
                format='ascii.fixed_width')

