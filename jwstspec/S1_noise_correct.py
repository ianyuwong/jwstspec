from glob import glob
import astropy.io.fits as fits
from astropy.stats import sigma_clip
import jwst.pipeline.calwebb_spec2 as pipe
import scipy.ndimage as nd
import numpy as np
import copy

def run(params):
	'''
	This function removes the 1/f readnoise artifacts from NIRSpec data

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
		# Get rate files produced by standard JWST pipeline
		input_files = sorted(glob(f'{params.data_dir}{params.prog_id}/Obs{oo}/Stage1/*rate.fits'))
		nfiles = len(input_files)
		if nfiles == 0:
			raise Warning('ERROR: no rate files found in Stage1 directory!')

		print(f'Stage 1: Running 1/f detector readnoise correction on {nfiles} files  from Obs{oo} using {params.readnoise_correct} method...')

		# Separate affix for different destriping methods
		if params.readnoise_correct == 'nsclean':
			vers = '0'
		elif params.readnoise_correct == 'constant':
			vers = '1' 
		elif params.readnoise_correct == 'moving_median':
			vers = '2'

		for i,fi in enumerate(input_files):
			print(f'Destriping {fi.split("/")[-1]} ({i+1}/{nfiles})...')

			# For NSClean, run the files through the first three steps of the calwebb_spec2 pipeline module and save outputs
			if params.readnoise_correct == 'nsclean':
				import os
				current_dir = os.path.dirname(__file__)
				cfg_file = os.path.join(current_dir, 'log.cfg')

				outdir = f'{params.data_dir}{params.prog_id}/Obs{oo}/Stage1/'
				result = pipe.assign_wcs_step.AssignWcsStep.call(fi, logcfg=cfg_file)
				if params.obs_type == 'ifu':
					result = pipe.msaflagopen_step.MSAFlagOpenStep.call(result, logcfg=cfg_file)
				result = pipe.nsclean_step.NSCleanStep.call(result, skip=False, save_results=True, output_dir=outdir, suffix='ratecorr0', logcfg=cfg_file)

			# Otherwise, manually destripe each file
			else:
				# Open rate file and deepcopy it to safely handle I/O
				hdulist = fits.open(fi)
				hdul = copy.deepcopy(hdulist)
				hdulist.close()
				data = hdul['SCI', 1].data

				# Flag bad pixels
				dq = hdul['DQ', 1].data
				spec_dq = dq > 0

				if params.obs_type == 'ifu':
					# Use cal file to get on-sky area to exclude for 1/f noise correction
					cal_file = fi.replace('Stage1/','cal/').replace('_rate','_cal')
					spec_mask = ~np.isnan(fits.getdata(cal_file))

					# Create a few pixels' worth of buffer
					spec_mask = nd.binary_dilation(spec_mask, iterations=2)

					# Also mask out fixed slit region at center of detectors
					spec_mask[900:1100,:] = True
				else:
					# For slit spectra, select outer rectangular regions to keep for destriping (top 5 and bottom 5 rows)
					spec_mask = np.zeros(data.shape)
					spec_mask[5:-5,:] = 1

				full_mask = np.logical_or(spec_mask, spec_dq)
				
				# Mask out on-sky pixels and remove stripes
				img = np.ma.MaskedArray(data, mask=full_mask)
				clipped = sigma_clip(img, axis=0, sigma=3.0, maxiters=None)
				if params.readnoise_correct == 'constant':
					stripes = np.array([np.ma.median(clipped, axis=0),]*img.shape[0])
				elif params.readnoise_correct == 'moving_median':
					import numbamisc
					stripes = numbamisc.median_filter(clipped,np.ones([201,1]))

				# Remove noise model and save as new file
				hdul['SCI',1].data = data - stripes
				hdul['PRIMARY',1].header['HISTORY'] = f'Readnoise correction carried out with S1_noise_correct.py using {params.readnoise_correct} method'
				hdul.writeto(fi.replace(f'_rate',f'_ratecorr{vers}'), overwrite=True)
				hdul.close()

	print('Stage 1 readnoise correction complete!')

	return params