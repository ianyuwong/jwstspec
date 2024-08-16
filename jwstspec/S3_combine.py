from glob import glob
from . import aux
import numpy as np

def run(params):
	'''
	This function cleans and combines spectra from one or more dither/nod positions.

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
	if params.special_defringe:
		spec_files = np.array(sorted(glob(inputdir+'*spectrum-defringed.pickle')))
	else:
		spec_files = np.array(sorted(glob(inputdir+'*spectrum.pickle')))

	resultsdir = inputdir

	# Parse gratings and detectors to group exposures together
	gratings = np.array([fi.split('_')[-4] for fi in spec_files])
	detectors = np.array([fi.split('_')[-3] for fi in spec_files])
	uniq_grat = np.unique(gratings)
	uniq_detectors = np.unique(detectors)
	print(f'Found gratings: {uniq_grat}')
	print(f'Found detectors: {uniq_detectors}')

	groups = []
	for grat in uniq_grat:
		for det in uniq_detectors:
			group = spec_files[np.where((gratings==grat) & (detectors==det))]
			if len(group) > 0:
				groups.append(group)

	# Combine spectra obtained in each filter
	for i,group in enumerate(groups):
		wave, flux, fluxerr, meta = aux.spec_combine(group, resultsdir, params.spec_bkg_sub, params.special_defringe, params.spec_sig_clip, params.spec_window_half)

		# Save spectrum
		aux.spec_save(wave, flux, fluxerr, resultsdir, meta)

	return params
