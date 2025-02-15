import astropy.io.fits as fits
from glob import glob
import numpy as np
from jwst.associations import asn_edit, asn_from_list, association_io
from jwst.associations.lib.rules_level3 import Asn_Lv3TSO
from jwst.pipeline.calwebb_tso3 import Tso3Pipeline
import os
from . import aux

current_dir = os.path.dirname(__file__)
cfg_file = os.path.join(current_dir, 'log.cfg')

def run(params):
	'''
	This function runs Stage 3 of the JWST pipeline for TSO observations.

	Parameters
    ----------
    params : obj
    	Object containing all parameter settings for pipeline run.
	'''

	# Get cal files produced by standard JWST pipeline or custom processing		
	if 'subbkg' in params.vers:		# If background subtraction was performed
		input_files = np.array(sorted(glob(f'{params.data_dir}{params.prog_id}/Obs{params.obs_numb}/Stage2{params.stage2_suffix}/*rate{params.vers}_calints.fits')))
	else:
		if params.instrument == 'nirspec':
			input_files = np.array(sorted(glob(f'{params.data_dir}{params.prog_id}/Obs{params.obs_numb}/Stage2{params.stage2_suffix}/*nrs?_calints.fits')))
		elif params.instrument == 'miri':
			input_files = np.array(sorted(glob(f'{params.data_dir}{params.prog_id}/Obs{params.obs_numb}/Stage2{params.stage2_suffix}/*[e]_calints.fits')))
	input_files = aux.select_spec_files(input_files, params)

	# Define output directory
	outdir = f'{params.data_dir}{params.prog_id}/Obs{params.obs_numb}/Stage3{params.stage3_suffix}/'
	os.makedirs(f'{outdir}', exist_ok=True)

	print(f'Stage 3: Running calwebb_tso3 pipeline...')
	# Make association
	filt = fits.getheader(input_files[0], 'PRIMARY')['FILTER']
	asn_name = f'jw{params.prog_id}-o{params.obs_numb}_t001_{params.instrument}_{filt.lower()}'
	asn = asn_from_list.asn_from_list(input_files, rule=Asn_Lv3TSO, format='json', product_name=asn_name)
	asn['asn_type'] = "tso-spec3"
	asn['program'] = params.prog_id
	asn['asn_id'] = f'o{params.obs_numb}'
	asn['target'] = fits.getval(input_files[0], 'TARGNAME')

	# Save ASN file	
	serialized = association_io.json.dump(asn)[1]
	output_file = f'{outdir}{asn_name}_asn.json'
	f = open(output_file, 'w')
	f.write(serialized)
	f.close()

	# Process ASN file through jwst pipeline module calwebb_tso3
	Tso3Pipeline.call(output_file, output_dir=outdir, save_results=True, steps=params.stage3_rules, logcfg=cfg_file)

	print('Stage 3 processing complete!')

	return params