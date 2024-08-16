from jwstspec import (S0_download, S1_stage1_process, 
			S1_bkg_subtract, S2_stage2_process, S3_stage3_process, 
			S3_stage3tso_process, aux)
import jwst
print('Using JWST pipeline version {}'.format(jwst.__version__))
import os

# Specify CRDS context
os.environ['CRDS_CONTEXT'] = 'jwst_1261.pmap'

# Initiate pipeline parameter object
params = aux.params()

## ==== Program and observation information ===
params.prog_id    = '04496'		# JWST program number (must be five-digit string)
params.obs_numb   = '017'		# Observation number (must be three-digit string)
params.bkg_obs_numb = None		# Observation number of dedicated background files (carry out nod subtraction if ==obs_numb). Use None if no background.
params.instrument = 'miri'		# 'miri'
params.obs_type   = 'slit'		# slit' or 'slitless' of 'ifu'
params.tso_observation = False 	# Time-series observation
## ============ Directory settings ============
params.dwnld_dir  = '/Users/iwong/Downloads/JWSTData/'								# User-set staging area for MAST downloads
params.data_dir   = '/Users/iwong/Documents/Astro/STScI/pipeline_testing/data/'	# Top-level directory for stored data and outputs
params.download   = True		# If True, download data from MAST
params.dwnld_all  = False		# If True, download all available data products. If False, only uncals are downloaded
## ========== Pipeline run settings ===========
params.bkg_subtract = None		# Dedicated background subtraction (None or 'pixel' or 'asn')
params.stage2_suffix = '_test0'	# Ending of Stage 2 file directory, used to distinguish between different data processing versions
params.stage3_suffix = '_test0'	# Ending of Stage 3 file directory, used to distinguish between different data processing versions
## ============================================

# Choose which steps to run
run_stage0 = False 			# Download data from MAST and initiate file directories
run_stage1 = False			# Produce countrate files, clean readnoise, handle background subtraction
run_stage2 = True			# Produce calibrated products for individual dithers/nods
run_stage3 = True			# Produce combined calibrated products

# Define custom rules for the various JWST pipeline stages

params.stage1_rules['jump'] = {'rejection_threshold' : 5.0}
# params.stage1_rules['emicorr'] = {'skip' : True}

params.stage2_rules['photom'] = {'override_photom' : '/Users/iwong/Documents/Astro/STScI/pipeline_testing/diagnostics/jwst_miri_photom_slit_240815.fits'}	
# params.stage2_rules['photom'] = {'override_photom' : '/Users/iwong/Documents/Astro/STScI/pipeline_testing/diagnostics/jwst_miri_photom_slitless_240815.fits'}
# # params.stage2_rules['photom'] = {'skip' : True}
# params.stage2_rules['assign_wcs'] = {'override_specwcs' : '/Users/iwong/Documents/Astro/STScI/pipeline_testing/diagnostics/MIRI_FM_MIRIMAGE_P750L_SLITLESSPRISM_DISTORTION_10.00.01.fits'}
# # params.stage2_rules['extract_1d'] = {'subtract_background' : True, "override_extract1d" : 'custom_lrs_extract1d.json'}	
params.stage2_rules['extract_1d'] = {'override_apcorr' : '/Users/iwong/Documents/Astro/STScI/pipeline_testing/diagnostics/MIRI_FM_MIRIMAGE_P750L_APCORR_240815_OPSF_flight_corr.fits'}
# # params.stage2_rules['pixel_replace'] = {'skip' : False}

# params.stage3_rules['outlier_detection'] = {'skip' : True}	
# params.stage3_rules['outlier_detection'] = {'save_intermediate_results' : True}	
# # params.stage3_rules = {'outlier_detection' : {'snr' : '15.0 15.0'}}	# Raise outlier detection thresholds
# # params.stage3_rules['extract_1d'] = {'subtract_background' : True, "override_extract1d" : 'custom_lrs_extract1d.json'}	
params.stage3_rules['extract_1d'] = {'override_apcorr' : '/Users/iwong/Documents/Astro/STScI/pipeline_testing/diagnostics/MIRI_FM_MIRIMAGE_P750L_APCORR_240815_OPSF_flight_corr.fits'}


## ============================================
## ============================================
## ============================================
print('## ============================================')
print('START PIPELINE PROCESSING')
print('## ============================================')

## Do parameter sanity check
params.sanity_check()

## STAGE 0
if run_stage0:
	# Download data from MAST
	params = S0_download.run(params)

## STAGE 1
if run_stage1:
	# Run calwebb_detector1 and produce countrate files
	params = S1_stage1_process.run(params)

	# Handle background subtraction, if needed
	if params.bkg_subtract is not False:
		params = S1_bkg_subtract.run(params)

## STAGE 2
if run_stage2:
	# Run calwebb_spec2 and produce calibrated products
	params = S2_stage2_process.run(params)

## STAGE 3
if run_stage3:
	# Run calwebb_spec3 or calwebb_tso3 and produce combined calibrated products
	if not params.tso_observation:
		params = S3_stage3_process.run(params)
	else:
		params = S3_stage3tso_process.run(params)

print('## ============================================')
print('END PIPELINE PROCESSING')
print('## ============================================')
