from glob import glob
from jwst.pipeline.calwebb_detector1 import Detector1Pipeline
import numpy as np
import astropy.io.fits as fits
import os
from . import aux

current_dir = os.path.dirname(__file__)
cfg_file = os.path.join(current_dir, 'log.cfg')

def run(params):
	'''
	This function runs Stage 1 of the JWST pipeline on the uncalibrated data.

	Parameters
    ----------
    params : obj
    	Object containing all parameter settings for pipeline run.
	'''

	# Run for both science and (if needed) background observations
	obs = [params.obs_numb]
	if params.bkg_obs_numb is not None and params.bkg_obs_numb != params.obs_numb:
		obs.append(params.bkg_obs_numb)

	for oo in obs:
		# Get all uncal files
		input_files = np.array(sorted(glob(f'{params.data_dir}{params.prog_id}/Obs{oo}/uncal/*_uncal.fits')))

		# Exclude target acquisition observations, (for MIRI) imaging, and (for NIRSpec) NRS2 uncal files for PRISM and M grating observations
		input_files = aux.select_spec_files(input_files, params)
		nfiles = len(input_files)

		print(f'Stage 1: processing {nfiles} files from Obs{oo}:')
		print(f'{input_files}')

		# Process each one through jwst pipeline module calwebb_detector1
		outdir = f'{params.data_dir}{params.prog_id}/Obs{oo}/Stage1/'
		for i,fi in enumerate(input_files):
			print(f'Processing file {i+1} of {nfiles}...')
			Detector1Pipeline.call(fi, output_dir=outdir, save_results=True, save_calibrated_ramp=False, steps=params.stage1_rules, logcfg=cfg_file)

	print('Stage 1 detector processing complete!')

	return params