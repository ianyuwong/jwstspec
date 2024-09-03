from glob import glob
import numpy as np
from jwst.spectral_leak.spectral_leak_step import SpectralLeakStep
from stdatamodels.jwst import datamodels

def run(params):
	'''
	This function applies spectral leak correction on MIRI MRS spectra.

	Parameters
    ----------
    params : obj
    	Object containing all parameter settings for pipeline run.
	'''

	# Define directories and appending text
	if params.extr_method == 'aperture':
		inputdir = f'{params.data_dir}{params.prog_id}/Obs{params.obs_numb}/outputs/{params.extr_aper_rad}_{params.bkg_aper_in}_{params.bkg_aper_out}{params.extr_suffix}/'
	elif params.extr_method == 'PSF':
		inputdir = f'{params.data_dir}{params.prog_id}/Obs{params.obs_numb}/outputs/PSF_{params.extr_aper_rad}_{params.bkg_aper_in}{params.extr_suffix}/'
	elif params.extr_method == 'box':
		inputdir = f'{params.data_dir}{params.prog_id}/Obs{params.obs_numb}/outputs/box_{params.extr_aper_rad}_{params.bkg_aper_in}{params.extr_suffix}/'
	append = ''
	if not params.spec_bkg_sub:
		append += '-nobkgsub'
	if params.special_defringe:
		append += '-defringed'

	# Get files that contains 1B and 3A sub-bands
	try:
		ch1b_file = glob(inputdir+f'*ch1*MEDIUM_MIRIFUSHORT_spectrum{append}.dat')[0]
		ch3a_file = glob(inputdir+f'*ch3*SHORT_MIRIFULONG_spectrum{append}.dat')[0]
	except:
		print('Missing Ch1b and/or Ch3a data. Skipping spectral leak correction.')
		return

	# Extract spectral segments
	ch1b_start = 5.65
	ch1b_end = 6.65
	ch3a_start = 11.55
	ch3a_end = 13.46

	ch1b_data = np.genfromtxt(ch1b_file)
	wave1b, flux1b = ch1b_data[:, 0], ch1b_data[:, 1]
	ind1b = np.where((wave1b > ch1b_start) * (wave1b < ch1b_end))
	wave1b, flux1b = wave1b[ind1b], flux1b[ind1b]

	ch3a_data = np.genfromtxt(ch3a_file)
	wave3a, flux3a = ch3a_data[:, 0], ch3a_data[:, 1]
	ind3a = np.where((wave3a > ch3a_start) * (wave3a < ch3a_end))
	wave3a, flux3a = wave3a[ind3a], flux3a[ind3a]

	# Get reference file
	x1d_file = glob(f'{params.data_dir}{params.prog_id}/Obs{params.obs_numb}/*/*x1d.fits')[0] 	# any x1d file for the observation
	spectral_leak_step = SpectralLeakStep()
	ref_file = spectral_leak_step.get_reference_file(x1d_file, 'mrsptcorr')

	# Extract spectral leak correction
	leak_ref = datamodels.MirMrsPtCorrModel(ref_file)
	leak_wave = leak_ref.leakcor_table.wavelength
	leak_percent = leak_ref.leakcor_table.frac_leak

	# Calculate leak correction and apply to 3A spectrum
	leak = np.interp(wave3a, 2*wave1b, flux1b) * np.interp(wave3a, leak_wave, leak_percent)
	flux3a_corr = flux3a - leak

	# Recover missing long wavelengths
	w = np.where((wave3a > 13.2) & (np.isnan(leak)))
	flux3a_corr[w] = flux3a[w]

	# Save corrected 3A spectrum
	ch3a_data[ind3a, 1] = flux3a_corr
	new_file = ch3a_file.split('.dat')[0] + '-specleakcorr.dat'
	header = '\n'.join(np.genfromtxt(ch3a_file, max_rows=16, delimiter='xx',dtype='str', comments='$')).replace('# ','')
	np.savetxt(new_file, ch3a_data, delimiter='\t', fmt='%.6f', header=header)

	print('MIRI spectral leak correction complete!')

	return params