# python standard library
import sys
import time
import os
from glob import glob as g
import xml.etree.ElementTree as ET
import logging
# astroquery http://ascl.net/1708.004
# cpl https://ascl.net/1612.001
from astroquery.eso import Eso
import cpl


# path to MUSE DRS recipes
cpl.Recipe.path = '/usr/lib64/esopipes-plugins/muse-2.8.3/'  # chapman2

# params
eso_username = 'echerenzops'  # ESO User portal username
# script requires authentication with the ESO User portal
# see https://astroquery.readthedocs.io/en/latest/eso/eso.html#authentication-with-eso-user-portal
date =  sys.argv[1]  # '2019-05-04'
mode =  sys.argv[2]  # 'WFM-NOAO-N'
star = 'LTT 3218'  # for now only using LTT 3218
outdir_prefix = '/homelocal2/eherenz/muse_diffphot/'
remove_and_zip = True # if true all immediate big files are removed, only
                      # pixtable (zipped) and cube and skyspectra remain on disk

# code 
outdir = outdir_prefix + 'rawdata_' + date + '_' + mode + '/'

eso = Eso()

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
    raise Exception('STD ' + star + ' observed in ' + mode + ' on ' +date + \
                    ', but no OBJECT observed in that mode during that night.')

dpid = obj_table['DP.ID'][0]

archive_xml_filename = outdir + dpid + '.xml'

# below will only work after download with query_standard_and_calib.py
assert os.path.isfile(archive_xml_filename)

# we parse the XML
tree = ET.parse(archive_xml_filename)
root = tree.getroot()
std_tree_xpath = "./associatedFiles/association/[@category='STD']"
file_xpath = "mainFiles/file"
assocfile_xpath = 'associatedFiles/'
assert len(root.findall(std_tree_xpath)) == 1  # there should be only one associated std!
std_tree = root.find(std_tree_xpath)

# the standard star exposure
std_exp_file = std_tree.find(file_xpath).get('name')
std_exp_category = std_tree.find(file_xpath).get('category')
# this dict we will use with python cpl then
dict_for_redux = {std_exp_category: ([std_exp_file + '.fits.fz'],
                                     None)}  # see next for loop for rationale

all_files_list = [std_exp_file + 'fits.fz']  # files not in this list
                                             # will be delted

# associated calibration data
std_assoc_tree = std_tree.findall(assocfile_xpath)

def file_list(elem, file_xpath):
    """make list of proper filenames from entries stored in XML"""
    def create_filename(e):
        d = {'MUSE': '.fits.fz', 'M': '.fits'}
        name = e.get('name')
        suffix = d[name[:name.find('.')]]
        return name + suffix
    return [create_filename(e) for e in elem.findall(file_xpath)]


for elem in std_assoc_tree:
    category = elem.get('category')  # ARC / BIAS etc.
    # filenames of the ARCs, the BIASEs etc. in a list
    main_files = file_list(elem, file_xpath)

    all_files_list += main_files
    # take care associated calibs for that category
    # (e.g. relevant for skyflat taken on different day)
    assoc_assoc = elem.findall(assocfile_xpath)
    if len(assoc_assoc) == 0:
        # output dict - values are tuples, second elem is None if
        # associated calibration data is standalone ...
        dict_for_redux[category] = (main_files, None)
    else:
        assoc_assoc_dict = {}
        for a in assoc_assoc:
            assoc_assoc_files = file_list(a, file_xpath)
            assoc_assoc_dict[a.get('category')] = assoc_assoc_files
            all_files_list += assoc_assoc_files

        # ... otherwise we store the association dictonary in the
        # second elem
        dict_for_redux[category] = (main_files,
                                    assoc_assoc_dict)


# cpl reduction
def append_outdir(file_list):
    return [outdir + filename for filename in file_list]


def logfile(logfilename):
    log = logging.getLogger()
    log.setLevel(logging.INFO)
    ch = logging.FileHandler(logfilename)
    ch.setLevel(logging.INFO)
    fr = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', '%H:%M:%S')
    ch.setFormatter(fr)
    log.addHandler(ch)


def get_timestring(starttime=None):
    assert type(starttime) == type(time.time())
    if starttime == None:
        return ''
    else:
        return '('+str(round(time.time()-starttime,3))+'s)'


# MASTER_BIAS
# -----------
def run_bias(master_bias_name, bias_files, badpix_table):
    muse_bias = cpl.Recipe('muse_bias', threaded=True)
    # non-default (but for us default) parameters
    muse_bias.param.nifu = -1
    muse_bias.param.merge = True
    # calibration frames
    muse_bias.calib['BADPIX_TABLE'] = badpix_table
    
    # call recipe
    print('muse_bias start for '+master_bias_name)
    starttime = time.time()
    logfile(master_bias_name+'.log')
    muse_bias_result = muse_bias(bias_files)  # actual call
    muse_bias_result.MASTER_BIAS.writeto(master_bias_name+'.fits')
    print('muse_bias done... Execution time '+get_timestring(starttime))
    print('muse_bias output written to '+master_bias_name+'.fits')


#  files
bias_files = append_outdir(dict_for_redux['BIAS'][0])
badpix_table = append_outdir(dict_for_redux['BIAS'][1]['BADPIX_TABLE'])
master_bias_name = 'MASTER_BIAS_'+date
print('Creating MASTER BIAS for '+date+' ...')
run_bias(master_bias_name, bias_files, badpix_table)


# MASTER_FLAT & TRACE_TABLE
# -------------------------
def run_flat(master_flat_name, trace_table_name, lamp_flat_files,
             master_bias, badpix_table):
    muse_flat = cpl.Recipe('muse_flat', threaded=True)
    # non-default (but for us default) parameters
    muse_flat.param.nifu = -1
    muse_flat.param.merge = True
    # calibration frames
    muse_flat.calib['BADPIX_TABLE'] = badpix_table
    muse_flat.calib['MASTER_BIAS'] = master_bias
    # call recipe
    print('muse_flat start for '+master_flat_name)
    starttime = time.time()
    logfile(master_bias_name+'.log')
    muse_flat_result = muse_flat(lamp_flat_files)
    muse_flat_result.MASTER_FLAT.writeto(master_flat_name+'.fits')
    muse_flat_result.TRACE_TABLE.writeto(trace_table_name+'.fits')
    print('muse_flat done... Execution time '+get_timestring(starttime))
    print('muse_flat output written to '+master_flat_name+'.fits & '+\
          trace_table_name+'.fits')
    
master_flat_name = 'MASTER_FLAT_'+date
trace_table_name = 'TRACE_TABLE_'+date

# lamp flat is associated to a different set of bias frames
# probably unnecessary? 
lamp_flat_bias_files = append_outdir(dict_for_redux['LAMP_FLAT'][1]['BIAS'])
print('Creating MASTER BIAS to be used for for LAMP FLAT '+date+' ..')
lamp_flat_bias_name = 'LAMP_FLAT_MASTER_BIAS_'+date
run_bias(lamp_flat_bias_name, lamp_flat_bias_files, badpix_table)
# lamp flat files
lamp_flat_files = append_outdir(dict_for_redux['LAMP_FLAT'][0])
print('Creating MASTER_FLAT & TRACE_TABLE for '+date+' ...')
run_flat(master_flat_name, trace_table_name,
         lamp_flat_files,
         lamp_flat_bias_name+'.fits', badpix_table)


# WAVECAL_TABLE
# -------------
def run_wavecal(wavecal_table_name, arc_files,
                 master_bias, trace_table, line_catalog, badpix_table):
    muse_wavecal = cpl.Recipe('muse_wavecal', threaded=True)
    # non-default (but for us default) parameters
    muse_wavecal.param.nifu = -1
    muse_wavecal.param.merge = True
    # calibration frames
    muse_wavecal.calib = {'BADPIX_TABLE': badpix_table,
                          'LINE_CATALOG': line_catalog,
                          'TRACE_TABLE': trace_table,
                          'MASTER_BIAS': master_bias}

    # call recipe
    print('muse_wavecal start for '+wavecal_table_name)
    starttime = time.time()
    logfile(wavecal_table_name+'.log')
    muse_wavecal_result = muse_wavecal(arc_files)
    muse_wavecal_result.WAVECAL_TABLE.writeto(wavecal_table_name+'.fits')
    print('muse_wavecal done... Execution time '+get_timestring(starttime))
    print('muse_wavecal output written to '+wavecal_table_name+'.fits')

wavecal_table_name = 'WAVECAL_TABLE_'+date
    
arc_files = append_outdir(dict_for_redux['ARC'][0])
line_catalog = append_outdir(dict_for_redux['ARC'][1]['LINE_CATALOG'])
print('Creating WAVECAL_TABLE for '+date+' ...')
run_wavecal(wavecal_table_name, arc_files,
            master_bias_name+'.fits', trace_table_name+'.fits', line_catalog, badpix_table)


# TWILIGHT_CUBE
# -------------
def run_twilight(twilight_cube_name, skyflat_files, skyflat_illum_file,
                 master_bias, master_flat, trace_table, wavecal_table,
                 geometry_table, badpix_table):
    muse_twilight = cpl.Recipe('muse_twilight', threaded=True)
    # calibration frames
    muse_twilight.calib = {'BADPIX_TABLE': badpix_table,
                           'GEOMETRY_TABLE': geometry_table,
                           'WAVECAL_TABLE': wavecal_table,
                           'TRACE_TABLE': trace_table,
                           'MASTER_FLAT': master_flat,
                           'MASTER_BIAS': master_bias}
    # call recipe
    print('muse_twilight start for '+twilight_cube_name)
    starttime = time.time()
    logfile(twilight_cube_name+'.log')
    muse_twilight_result = muse_twilight(raw={'SKYFLAT': skyflat_files,
                                             'ILLUM': skyflat_illum_file})
    muse_twilight_result.DATACUBE_SKYFLAT.writeto(twilight_cube_name+'.fits')
    print('muse_twilight done... Execution time '+get_timestring(starttime))
    print('muse_twilight output written to '+twilight_cube_name+'.fits')


twilight_cube_name = 'TWILIGHT_CUBE_'+date
skyflat_files = append_outdir(dict_for_redux['SKY_FLAT'][0])
skyflat_illum_file = append_outdir(dict_for_redux['SKY_FLAT'][1]['ILLUM'])
geometry_table = append_outdir(dict_for_redux['SKY_FLAT'][1]['GEOMETRY_TABLE'])
print('Creating SKY FLAT for '+date+' ....')

# bias, flat, trace table, and wavecal table could be generated from files associated
# - but this is usually stable enough 
run_twilight(twilight_cube_name, skyflat_files, skyflat_illum_file,
             master_bias_name+'.fits', master_flat_name+'.fits',
             trace_table_name+'.fits', wavecal_table_name+'.fits',
             geometry_table, badpix_table)

# SCIBASIC
# --------
def run_scibasic(std_pixtable_name, std_exp_filename, std_illum_file,
                 master_bias, master_flat, trace_table, wavecal_table, twilight_cube,
                 geometry_table, badpix_table):
    muse_scibasic = cpl.Recipe('muse_scibasic', threaded=True)
    # confirmed bug by PW, pixtables do not get combined automatically
    combine_pixtables = cpl.Recipe('muse_scipost_combine_pixtables', threaded=True)
    # non-default (but for us default) parameters
    muse_scibasic.param.nifu = -1
    muse_scibasic.param.merge = True
    # calibration frames
    muse_scibasic.calib = {'BADPIX_TABLE': badpix_table,
                           'GEOMETRY_TABLE': geometry_table,
                           'TWILIGHT_CUBE': twilight_cube,
                           'WAVECAL_TABLE': wavecal_table,
                           'TRACE_TABLE': trace_table,
                           'MASTER_FLAT': master_flat,
                           'MASTER_BIAS': master_bias}
    # call recipe
    print('muse_scibasic start for '+std_pixtable_name)
    starttime = time.time()
    logfile(std_pixtable_name+'.log')
    muse_scibasic_result = muse_scibasic(raw={'STD': std_exp_filename,
                                              'ILLUM': std_illum_file})
    
    pixtable_combined = combine_pixtables(muse_scibasic_result.PIXTABLE_STD)

    pixtable_combined.PIXTABLE_COMBINED.writeto('combined_'+std_pixtable_name+'.fits')
    print('muse_scibasic done... Execution time '+get_timestring(starttime))
    print('muse_scibasic output written to '+std_pixtable_name+'.fits')

    return pixtable_combined.PIXTABLE_COMBINED
    
std_pixtable_name = star.replace(" ", "_") + '_' + mode.lower() + '_' + date + '_PIXTABLE'
std_exp_filename = append_outdir(dict_for_redux['STD'][0])
std_illum_file = append_outdir(dict_for_redux['ILLUM'][0])
geometry_table = append_outdir(dict_for_redux['GEOMETRY_TABLE'][0])

print("Creating PIXTABLE for "+std_pixtable_name)
pixtable_combined = run_scibasic(std_pixtable_name, std_exp_filename, std_illum_file,
                                 master_bias_name+'.fits',
                                 master_flat_name+'.fits',
                                 trace_table_name+'.fits',
                                 wavecal_table_name+'.fits', twilight_cube_name+'.fits',
                                 geometry_table, badpix_table)


# arbitrary FLUX CALIBRATION / EXTINCTION CORRECTION / SKY SUBTRACTION
# --------------------------------------------------------------------
def run_scipost(outcube_name, pixtable_combined_name):
    muse_scipost = cpl.Recipe('muse_scipost', threaded=True)
    static_cal_dir = '/usr/share/esopipes/datastatic/muse-2.8.3/'
    muse_scipost.calib = \
        {'EXTINCT_TABLE': static_cal_dir + 'extinct_table.fits',
         'SKY_LINES': static_cal_dir + 'sky_lines.fits',
         'FILTER_LIST': static_cal_dir + 'filter_list.fits',
         'STD_RESPONSE': static_cal_dir + 'std_response_' + mode.lower() + '.fits',
         'LSF_PROFILE': static_cal_dir + 'lsf_profile_slow_' + mode.lower() + '.fits',
         'ASTROMETRY_WCS': static_cal_dir + 'astrometry_wcs_wfm.fits'}#,
#         'OUTPUT_WCS': './static_cal_flx/LTT_3218_WCS_TEMPLATE_DATACUBE.fits'}
    print('muse_scipost start for ' + outcube_name)
    starttime = time.time()
    logfile(outcube_name+'.log')
    muse_scipost_result = muse_scipost(raw={'PIXTABLE_OBJECT': pixtable_combined_name})
    muse_scipost_result.DATACUBE_FINAL.writeto(outcube_name+'.fits')
    muse_scipost_result.IMAGE_FOV.writeto('IMAGE_FOV_'+outcube_name+'.fits')
    muse_scipost_result.SKY_MASK.writeto('SKY_MASK_'+outcube_name+'.fits')
    muse_scipost_result.SKY_SPECTRUM.writeto('SKY_SPECTRUM_'+outcube_name+'.fits')
    muse_scipost_result.SKY_LINES.writeto('SKY_LINES_'+outcube_name+'.fits')
    muse_scipost_result.SKY_CONTINUUM.writeto('SKY_CONTINUUM_'+outcube_name+'.fits')
    print('muse_scipost done... Execution time '+get_timestring(starttime))
    print('muse_scipost output written to '+outcube_name+'.fits')
    
# MAKE non-fluxcal / non-skysub CUBE
# ----------------------------------
outcube_name = star.replace(' ', '_') + mode.lower() + '_' + date
run_scipost(outcube_name, 'combined_'+std_pixtable_name+'.fits')

print("Resampling non-flux calibrated / non sky-subtracted pixtable to datacube ... ")
make_cube = cpl.Recipe('muse_scipost_make_cube', threaded=True)
#make_cube.calib['OUTPUT_WCS'] = './static_cal_flx/LTT_3218_WCS_TEMPLATE_DATACUBE.fits'
logfile(outcube_name + '_basic.log')
c = make_cube(pixtable_combined)
c.DATACUBE_FINAL.writeto(outcube_name +  '_basic.fits')
print("DONE!!!! Written "+outcube_name + '_basic.fits to disk!')

if remove_and_zip:
    # remove (or zip) intermediate files
    os.remove(master_bias_name+'.fits')
    os.remove(master_flat_name+'.fits')
    os.remove(trace_table_name+'.fits')
    os.remove(lamp_flat_bias_name+'.fits')
    os.remove(wavecal_table_name+'.fits')
    os.remove(twilight_cube_name+'.fits')
    os.system('pigz -9 combined_'+std_pixtable_name+'.fits')

# END
