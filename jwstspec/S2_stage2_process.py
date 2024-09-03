import os
from glob import glob
import numpy as np
from jwst.pipeline.calwebb_spec2 import Spec2Pipeline
from . import aux

current_dir = os.path.dirname(__file__)
cfg_file = os.path.join(current_dir, 'log.cfg')

def run(params):
	'''
	This function runs Stage 2 of the JWST pipeline on the countrate data.

	Parameters
    ----------
    params : obj
    	Object containing all parameter settings for pipeline run.
	'''

	# Get rate(ints) files from Stage 1 processing
	if params.bkg_subtract == 'pixel' or params.bkg_subtract == None:
		input_files = aux.select_spec_files(np.array(sorted(glob(f'{params.data_dir}{params.prog_id}/Obs{params.obs_numb}/Stage1/*_rate{params.vers}.fits'))), params)
	elif params.bkg_subtract == 'asn':
		input_files = np.array(sorted(glob(f'{params.data_dir}{params.prog_id}/Obs{params.obs_numb}/Stage1/*_rate{params.vers}.json')))
	nfiles = len(input_files)

	# Process each one through jwst pipeline module calwebb_spec2
	print(f'Stage 2: Calibrating frames and building Level 2 spectrum products for {nfiles} files...')
	outdir = f'{params.data_dir}{params.prog_id}/Obs{params.obs_numb}/Stage2{params.stage2_suffix}/'
	os.makedirs(f'{outdir}', exist_ok=True)

	for i,fi in enumerate(input_files):
		print(f'Processing file {i+1} of {nfiles}...')

		# Custom step rules for various input types and specifications
		if params.obs_type == 'ifu':
			# Turn autocentroiding on
			try:
				params.stage2_rules['extract_1d'].update({'ifu_autocen' : True})
			except:
				params.stage2_rules['extract_1d'] = {'ifu_autocen' : True}

			# Set other parameters using pipeline settings
			if params.instrument == 'miri':
				params.stage2_rules['residual_fringe'] = {'skip' : False}				# Make sure residual_fringe runs on MIRI MRS datasets
				params.stage2_rules['extract_1d'].update({'ifu_rfcorr' : True}) 		# Turn on extra spectrum-level defringing step
			if hasattr(params, 'cube_align'):
				if params.cube_align == 'ifu':
					params.stage2_rules['cube_build'] = {'coord_system' : 'ifualign'}
				elif params.cube_align == 'internal':
					params.stage2_rules['cube_build'] = {'coord_system' : 'internal_cal'}
		if params.instrument == 'nirspec':
			params.stage2_rules['nsclean'] = {'skip' : True}			# Make sure to skip this (for now), since NSClean is currently handled in S1_noise_correct.py
		if params.bkg_subtract == 'asn':
			params.stage2_rules['bkg_subtract'] = {'skip' : False}		# Make sure pipeline runs the ASN background subtraction routine

		# Pipeline call
		Spec2Pipeline.call(fi, output_dir=outdir, save_results=True, steps=params.stage2_rules, logcfg=cfg_file)

	print('Stage 2 spectrum build complete!')

	return params