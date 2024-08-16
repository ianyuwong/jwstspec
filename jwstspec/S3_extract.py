from glob import glob
from . import aux
import numpy as np

def run(params):
	'''
	This function extracts spectra from either Stage2 and Stage3 spectral images or cubes.

	Parameters
    ----------
    params : obj
    	Object containing all parameter settings for pipeline run.
	'''

	# Get s2d/s3d files produced by standard JWST pipeline or custom processing
	if params.obs_type == 'ifu':
		affix = 's3d'
	elif params.obs_type == 'slit':
		affix = 's2d'

	if params.extract_stage == 'Stage2':			
		if len(params.vers) > 0:
			input_files = sorted(glob(f'{params.data_dir}{params.prog_id}/Obs{params.obs_numb}/Stage2{params.stage2_suffix}/*rate{params.vers}_{affix}.fits'))
		else:
			if params.instrument == 'nirspec':
				input_files = sorted(glob(f'{params.data_dir}{params.prog_id}/Obs{params.obs_numb}/Stage2{params.stage2_suffix}/*nrs?_{affix}.fits'))
			elif params.instrument == 'miri':
				input_files = sorted(glob(f'{params.data_dir}{params.prog_id}/Obs{params.obs_numb}/Stage2{params.stage2_suffix}/*[n,r,g][g,t,e]_{affix}.fits'))
	elif params.extract_stage == 'Stage3':
		input_files = sorted(glob(f'{params.data_dir}{params.prog_id}/Obs{params.obs_numb}/Stage3{params.stage3_suffix}/*_{affix}.fits'))
	nfiles = len(input_files)

	print(f'Stage 3: Extracting spectra from {nfiles} files...')

	# Read in data
	if params.obs_type == 'ifu':
		input_set = aux.read_cubes(input_files, suffix=params.extr_suffix)
	if params.obs_type == 'slit':
		# Indicate if pairwise nod subtraction was carried out
		if params.bkg_obs_numb == params.obs_numb:
			nod_subtract = True
		else:
			nod_subtract = False
		input_set = aux.read_spectra(input_files, suffix=params.extr_suffix, nod_subtract=nod_subtract)

	# cent_data = np.genfromtxt('centroids_bkgfromcharon2_miri.dat')

	for i,inp in enumerate(input_set):
		# Extract spectrum and save it
		print("================")
		print(f'Processing file {i+1} of {nfiles}...')
		# fix_centroid = [cent_data[i][0], cent_data[i][1]]
		inp.extract_spec(extr_aper_rad=params.extr_aper_rad, bkg_aper_in=params.bkg_aper_in, bkg_aper_out=params.bkg_aper_out, fix_centroid=params.fix_centroid, \
			bkg_sig_clip=params.bkg_sig_clip, pix_sig_clip=params.pix_sig_clip, extr_method=params.extr_method, window_width=params.window_width, save_cleaned=params.save_cleaned)

		# Create wavelength-stacked plot
		if params.extr_method == 'aperture':
			inp.master_plot()

	print('Stage 3 spectral extraction complete!')

	return params