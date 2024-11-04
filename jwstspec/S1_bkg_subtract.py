import os
from glob import glob
import astropy.io.fits as fits
import numpy as np
from jwst.associations import asn_from_list, association_io
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase
from jwst.background.background_sub import background_sub
from jwst import datamodels
from . import aux

def run(params):
	'''
	This function handles background subtraction for all modes, either through direct pixel-by-pixel subtraction or ASN file creation.

	Parameters
    ----------
    params : obj
    	Object containing all parameter settings for pipeline run.
	'''

	# Define input science and background file directories
	sci_dir = f'{params.data_dir}{params.prog_id}/Obs{params.obs_numb}/Stage1/'
	bkg_dir = f'{params.data_dir}{params.prog_id}/Obs{params.bkg_obs_numb}/Stage1/'
	if params.obs_numb == params.bkg_obs_numb:	# Set nod = True if using pairwise dither subtraction from same observation
		nod = True 
		bkg_ind = [1,0,3,2]						# Index handling for pairwise subtraction
	else:
		nod = False

	# Get files
	vers = ''
	if params.tso_observation:
		vers += 'ints'
	if hasattr(params, 'readnoise_correct'):
		if params.readnoise_correct == 'nsclean':
			vers = 'corr0'
		elif params.readnoise_correct == 'constant':
			vers = 'corr1' 
		elif params.readnoise_correct == 'moving_median':
			vers = 'corr2'
	sci_files = aux.select_spec_files(np.array(sorted(glob(f'{sci_dir}*_rate{vers}.fits'))), params)
	bkg_files = aux.select_spec_files(np.array(sorted(glob(f'{bkg_dir}*_rate{vers}.fits'))), params)
	nscifiles = len(sci_files)
	nbkgfiles = len(bkg_files)
	if nscifiles == 0:
		raise Warning(f'ERROR: no science data files found in Obs{params.obs_numb} Stage1 directory!')
	if nbkgfiles == 0:
		raise Warning(f'ERROR: no background data files found in Obs{params.bkg_obs_numb} Stage1 directory!')

	# Sort by band/grating and detector
	if params.obs_type == 'ifu':
		if params.instrument == 'nirspec':
			key = 'GRATING'
		else:
			key = 'BAND'
		# Use grating and detector
		bands = np.array(['SHORT','SHORT','MEDIUM','MEDIUM','LONG','LONG','PRISM','G140M','G235M','G395M','G140H','G140H','G235H','G235H','G395H','G395H'])
		detectors = np.array(['MIRIFUSHORT','MIRIFULONG','MIRIFUSHORT','MIRIFULONG','MIRIFUSHORT','MIRIFULONG','NRS1','NRS1','NRS1','NRS1','NRS1','NRS2','NRS1','NRS2','NRS1','NRS2'])
	else:
		key = 'FILTER'
		# Use filter and detector
		bands = np.array(['CLEAR','P750L'])
		detectors = np.array(['NRS2','MIRIMAGE'])
	ncomb = len(bands)
	sci_groups = [np.zeros(0).astype('str'), ] * ncomb
	bkg_groups = [np.zeros(0).astype('str'), ] * ncomb

	# Populate groups and collect grating wheel assembly tilt values (for NIRSpec)
	sci_gwa_xtil = np.zeros(ncomb)		
	sci_gwa_ytil = np.zeros(ncomb)
	bkg_gwa_xtil = np.zeros(ncomb)		
	bkg_gwa_ytil = np.zeros(ncomb)
	for fi in sci_files:
		header_sci = fits.getheader(fi, 'PRIMARY')
		band, detector = header_sci[key], header_sci['DETECTOR']
		w = np.where((bands == band) & (detectors == detector))[0][0]
		sci_groups[w] = np.append(sci_groups[w],fi)
		if params.instrument == 'nirspec':
			sci_gwa_xtil[w] = header_sci['GWA_XTIL']
			sci_gwa_ytil[w] = header_sci['GWA_YTIL']
	for fi in bkg_files:
		header_bkg = fits.getheader(fi, 'PRIMARY')
		band, detector = header_bkg[key], header_bkg['DETECTOR']
		w = np.where((bands == band) & (detectors == detector))[0][0]
		bkg_groups[w] = np.append(bkg_groups[w],fi)
		if params.instrument == 'nirspec':
			bkg_gwa_xtil[w] = header_bkg['GWA_XTIL']
			bkg_gwa_ytil[w] = header_bkg['GWA_YTIL']

	# For NIRSpec, check to see if grating wheel tilt values allow for association-based background subtraction
	if params.instrument == 'nirspec' and params.bkg_subtract == 'asn':
		tolerance = 1.0e-8		# Same tolerance as JWST calibration pipeline
		diff_xtil = np.absolute(sci_gwa_xtil - bkg_gwa_xtil)
		diff_ytil = np.absolute(sci_gwa_ytil - bkg_gwa_ytil)

		# If tolerance not met, use manual pixel-by-pixel subtraction instead
		if np.max(diff_xtil) > tolerance:
			print('WARNING: Relative grating wheel assembly tilts for science and background observations not within tolerance!')
			print('Switching over to bkg_subtract == "pixel" method')
			params.bkg_subtract = 'pixel'
			params.vers.replace('asn', 'pixel')

	## Method 1: Subtract average background pixel-by-pixel for each science file
	if params.bkg_subtract == 'pixel':
		print(f'Stage 1: Subtracting background using Obs{params.obs_numb} (science) and Obs{params.bkg_obs_numb} (background)...')
		for i in range(ncomb):
			if len(sci_groups[i]) == 0:
				continue
			for j, sci_file in enumerate(sci_groups[i]):
				sci_model = datamodels.open(sci_file)
				if not nod:
					bkg_list = bkg_groups[i]
				elif nod:	
					ind = bkg_ind[j]
					bkg_list = bkg_groups[i][ind:ind+1]
				bkg_model, result = background_sub(sci_model, bkg_list, 3.0, None)		# Default step settings in Spec2Pipeline
				out_file = sci_file.replace(f'_rate{vers}',f'_rate{params.vers}')

				hdulist = fits.open(sci_file)
				hdulist['SCI',1].data = result.data
				hdulist['ERR',1].data = result.err
				hdulist['DQ',1].data = result.dq
				hdulist['PRIMARY',1].header['HISTORY'] = f'Pixel-by-pixel background subtraction carried out by S1_pixel_bkg_subtract.py using Obs{params.bkg_obs_numb}'
				hdulist.writeto(out_file, overwrite=True)
				hdulist.close()

		print('Stage 1 background subtraction complete!')

	## Method 2: Create associations for each pair of science-background files and save ASN files
	elif params.bkg_subtract == 'asn':
		print(f'Stage 1: Creating ASN files using Obs{params.obs_numb} (science) and Obs{params.bkg_obs_numb} (background)...')
		for i in range(ncomb):
			if len(sci_groups[i]) == 0:
				continue
			for j in range(len(sci_groups[i])):
				# This is a default Level 2 association
				asn_name = sci_groups[i][j].split(f'.fits')[0].split('/')[-1].replace(f'_rate{vers}',f'_rate{params.vers}')
				asn = asn_from_list.asn_from_list([sci_groups[i][j]], rule=DMSLevel2bBase, format='json')
				asn['asn_type'] = "spec2"
				asn['program'] = params.prog_id
				asn['asn_id'] = f'o{params.obs_numb}'

				# Add all appropriate background files (have to do this manually, due to directory differences)
				if not nod:
					for k in range(len(bkg_groups[i])):
						asn['products'][0]['members'].append({"expname": f"{bkg_groups[i][k]}","exptype": "background"})
				elif nod:	# There should just be a pair in the case of nod == True
					ind = int(not(bool(j)))
					asn['products'][0]['members'].append({"expname": f"{bkg_groups[i][ind]}","exptype": "background"})

				# Save ASN file
				serialized = association_io.json.dump(asn)[1]
				output_file = f'{sci_dir}{asn_name}.json'
				f = open(output_file, 'w')
				f.write(serialized)
				f.close()
				
		print('Stage 1 ASN file creation complete!')

	return params
