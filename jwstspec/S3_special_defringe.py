from glob import glob
import numpy as np
from jwst.residual_fringe import utils as rfutils
import pickle
import copy

def run(params):
	'''
	This function runs the spectrum-level defringing routine on MIRI MRS spectra.

	Parameters
    ----------
    params : obj
    	Object containing all parameter settings for pipeline run.
	'''

	# Grab spectrum files
	if params.extr_method == 'aperture':
		inputdir = f'{params.data_dir}{params.prog_id}/Obs{params.obs_numb}/outputs/{params.extr_aper_rad}_{params.bkg_aper_in}_{params.bkg_aper_out}{params.extr_suffix}/'
	elif params.extr_method == 'PSF':
		inputdir = f'{params.data_dir}{params.prog_id}/Obs{params.obs_numb}/outputs/PSF_{params.extr_aper_rad}_{params.bkg_aper_in}{params.extr_suffix}/'
	elif params.extr_method == 'box':
		inputdir = f'{params.data_dir}{params.prog_id}/Obs{params.obs_numb}/outputs/box_{params.extr_aper_rad}_{params.bkg_aper_in}{params.extr_suffix}/'
	spec_files = np.array(sorted(glob(inputdir+'*spectrum.pickle')))

	for i, file in enumerate(spec_files):
		print(f'Defringing {file} ({i+1}/{len(spec_files)})')
		with open(file,'rb') as f:
			spectrum = pickle.load(f)
		wave = spectrum.wave
		if params.spec_bkg_sub:
			flux = spectrum.flux
		else:
			flux = spectrum.fluxbkg
		detector = spectrum.detector
		nonzero_inds = np.where(flux != 0)[0]

		# For Stage 2 products, extract the two channel segments within the spectrum
		if params.extract_stage == 'Stage2':
			split_ind = (nonzero_inds[1:] - nonzero_inds[:-1]).argmax()
			seg_endpoints = [[nonzero_inds[0], nonzero_inds[split_ind]+1], [nonzero_inds[split_ind+1], nonzero_inds[-1]+1]]
			if detector == 'MIRIFUSHORT':
				channels = ['1', '2']
			elif detector == 'MIRIFULONG':
				channels = ['3', '4']

		# For Stage 3, the inputs are already separated into segments
		elif params.extract_stage == 'Stage3':
			channels = [spectrum.grating[2]]
			seg_endpoints = [[nonzero_inds[0], nonzero_inds[-1]+1]]

		# Run each channel segment through the defringing routine
		corrflux = copy.deepcopy(flux)
		for j, ch in enumerate(channels):
			start, end = seg_endpoints[j]
			try:
				corrflux[start:end] = rfutils.fit_residual_fringes_1d(flux[start:end], wave[start:end], channel=ch, dichroic_only=False, max_amp=None)
			except:
				print(f'Defringing failed in channel {ch}')

		# Save defringed flux into new pickle file
		if params.spec_bkg_sub:
			spectrum.flux = corrflux
		else:
			spectrum.fluxbkg = corrflux
		new_file = file.split('.pickle')[0] + '-defringed.pickle'
		with open(new_file,'wb') as f:
			pickle.dump(spectrum, f)

	return params