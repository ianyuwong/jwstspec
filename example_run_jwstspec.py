from jwstspec import (S0_download, S1_stage1_process, S1_noise_correct, 
			S1_bkg_subtract, S2_stage2_process, S3_stage3_process, 
			S3_stage3tso_process, S3_extract, S3_combine, 
			S3_special_defringe, S3_miri_spectral_leak, aux)
import jwst
print('Using JWST pipeline version {}'.format(jwst.__version__))
import os

# # Specify CRDS context
# os.environ['CRDS_CONTEXT'] = 'jwst_1261.pmap'

## ==== Program and observation information ===
prog_id    = '04496'		# JWST program number (must be five-digit string)
obs_numb   = '017'			# Observation number (must be three-digit string)
bkg_obs_numb = None			# Observation number of dedicated background files (carry out nod subtraction if ==obs_numb). Use None if no background.
instrument = 'miri'			# 'nirspec' or 'miri'
obs_type   = 'slit'			# 'ifu' or 'slit' or 'slitless'
tso_observation = False 	# Time-series observation
## ============ Directory settings ============
dwnld_dir  = '/Users/iwong/Downloads/JWSTData/'								# User-set staging area for MAST downloads
data_dir   = '/Users/iwong/Documents/Astro/STScI/pipeline_testing/data/'	# Top-level directory for stored data and outputs
download   = True			# If True, download data from MAST
dwnld_all  = False			# If True, download all available data products. If False, only uncals are downloaded
## ========== Pipeline run settings ===========
readnoise_correct = False	# If True, correct for 1/f detector readnoise in NIRSpec data
destriping = 'nsclean'		# 'nsclean' is the pipeline-internal method, which is selected by default; otherwise, this can be None, 'constant' or 'moving_median'
bkg_subtract = 'asn'		# Dedicated background subtraction (None or 'pixel' or 'asn')
cube_align  = None			# If None, use default (sky coordinates). If 'ifu', use detector coordinates. If 'internal', keep cube aligned with IFU slices
stage2_suffix = '_test0'	# Ending of Stage 2 file directory, used to distinguish between different data processing versions
stage3_suffix = '_test0'	# Ending of Stage 3 file directory, used to distinguish between different data processing versions
## ==== Spectral extraction settings ==========
extract_stage = 'Stage2' 	# Stage of data files on which spectral extraction is performed
extr_method   = 'PSF'		# PSF fitting ('PSF'), box ('box', only relevant for IFU data), or aperture ('aperture') extraction
extr_suffix	  = ''	 		# Custom output directory suffix to distinguish between different extraction versions
extr_aper_rad = 3			# Radius/half-length of circular/square extraction aperture centered around centroid 
bkg_aper_in   = 10			# Inner dimension of background region
bkg_aper_out  = 10			# Outer dimension of background region
window_width  = 10  		# Half-width of sliding window for PSF fitting (10 is default)
pix_sig_clip  = 5.0 		# Sigma clip threshold for IFU cube pixel array outlier flagging (should generally be >4 to be safe)
bkg_sig_clip  = 5.0 		# Sigma clip threshold for background outlier flagging (5.0 is the default)
fix_centroid  = None 		# Specified centroid position (None is default; use strings for relative offsets)
save_cleaned  = False		# If True, save cleaned data cubes (only relevant for IFU data)
## ====== Post-processing parameters ==========
spec_bkg_sub     = True		# Choose whether to use background-subtracted fluxes (should always be True for PSF fitting)
spec_sig_clip    = 5.0  	# Sigma clip threshold for spectrum trimming
spec_window_half = 10      	# Half-width of moving median filter window
special_defringe = False	# Carry out custom spectrum-level defringing (only applicable to MIRI MRS data)
## ============================================

# Choose which steps to run
run_stage0 = False 			# Download data from MAST and initiate file directories
run_stage1 = False			# Produce countrate files, clean readnoise, handle background subtraction
run_stage2 = True			# Produce calibrated products for individual dithers/nods
run_stage3 = True			# Produce combined calibrated products

run_extract = False			# Extract spectra
run_combine = False			# Combine spectra

# Initiate pipeline parameter object
params = aux.params()

# Define custom rules for the various JWST pipeline steps
params.stage1_rules['jump'] = {'rejection_threshold' : 5.0}
# params.stage1_rules['emicorr'] = {'skip' : False}

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

# Populate parameter values and do sanity checks
params.add_params(prog_id, obs_numb, bkg_obs_numb, instrument, obs_type, tso_observation, dwnld_dir, data_dir, download, dwnld_all,
			readnoise_correct, destriping, bkg_subtract, cube_align, stage2_suffix, stage3_suffix, extract_stage, extr_method, extr_suffix, extr_aper_rad,
			bkg_aper_in, bkg_aper_out, window_width, pix_sig_clip, bkg_sig_clip, fix_centroid, save_cleaned, spec_bkg_sub, spec_sig_clip, spec_window_half, special_defringe)

####### RUN PIPELINE #######
print('## ============================================')
print('START PIPELINE PROCESSING')
print('## ============================================')

## STAGE 0
if run_stage0:
	# Download data from MAST
	params = S0_download.run(params)

## STAGE 1
if run_stage1:
	# Run calwebb_detector1 and produce countrate files
	params = S1_level1_process.run(params)

	# Correct for readnoise, if needed
	if params.instrument == 'nirspec' and params.readnoise_correct == True:
		params = S1_noise_correct.run(params)

	# Handle background subtraction, if needed
	if params.bkg_subtract is not False:
		params = S1_bkg_subtract.run(params)

## STAGE 2
if run_stage2:
	# Run calwebb_spec2 and produce calibrated products
	params = S2_level2_process.run(params)

## STAGE 3
if run_stage3:
	# Run calwebb_spec3 and produce combined calibrated products
	if not params.tso_observation:
		params = S3_level3_process.run(params)
	else:
		params = S3_level3tso_process.run(params)

## Extract spectra
if run_extract:
	params = S3_extract.run(params)

	# For MIRI MRS, run special defringing routine, if needed
	if params.special_defringe and params.instrument == 'miri' and params.obs_type == 'ifu':
		params = S3_special_defringe.run(params)
			
## Combine spectra
if run_combine:
	params = S3_combine.run(params)

	# For MIRI MRS, correct for spectral leak
	if params.instrument == 'miri' and params.obs_type == 'ifu':
		params = S3_miri_spectral_leak.run(params)

print('## ============================================')
print('END PIPELINE PROCESSING')
print('## ============================================')
