import astropy.io.fits as fits
from glob import glob
import numpy as np
from jwst.associations import asn_edit, asn_from_list, association_io
from jwst.associations.lib.rules_level3 import Asn_Lv3NRSIFU, Asn_Lv3NRSFSS, Asn_Lv3MIRMRS, Asn_Lv3SpectralTarget
from jwst.pipeline.calwebb_spec3 import Spec3Pipeline
import os
from . import aux

current_dir = os.path.dirname(__file__)
cfg_file = os.path.join(current_dir, 'log.cfg')

def run(params):
	'''
	This function runs Stage 3 of the JWST pipeline.

	Parameters
    ----------
    params : obj
    	Object containing all parameter settings for pipeline run.
	'''

	# Get cal files produced by standard JWST pipeline or custom processing		
	if len(params.vers) > 0:
		input_files = np.array(sorted(glob(f'{params.data_dir}{params.prog_id}/Obs{params.obs_numb}/Stage2{params.stage2_suffix}/*rate{params.vers}_cal.fits')))
	else:
		if params.instrument == 'nirspec':
			input_files = np.array(sorted(glob(f'{params.data_dir}{params.prog_id}/Obs{params.obs_numb}/Stage2{params.stage2_suffix}/*nrs?_cal.fits')))
		elif params.instrument == 'miri':
			input_files = np.array(sorted(glob(f'{params.data_dir}{params.prog_id}/Obs{params.obs_numb}/Stage2{params.stage2_suffix}/*[n,r,g][g,t,e]_cal.fits')))
	input_files = aux.select_spec_files(input_files, params)

	# Sort by grating (NIRSpec only)
	if params.instrument == 'nirspec':
		gratings = np.array(['PRISM','G140M','G235M','G395M','G140H','G235H','G235H','G395H'])
		ngrat = len(gratings)
		sci_groups = [np.zeros(0).astype('str'), ] * ngrat
		for fi in input_files:
			header_sci = fits.getheader(fi, 'PRIMARY')
			grating = header_sci['GRATING']
			w = np.where((gratings == grating))[0][0]
			sci_groups[w] = np.append(sci_groups[w],fi)

		# Only keep groups that have files
		used_gratings = [gratings[i] for i in range(ngrat) if len(sci_groups[i]) > 0]
		sci_groups = [grp for grp in sci_groups if len(grp) > 0]
		ngroups = len(sci_groups)

	# For MIRI, group everything together
	elif params.instrument == 'miri':
		sci_groups = [input_files]
		if params.obs_type == 'ifu':
			used_gratings = ['COMBINED']
		else:
			used_gratings = ['P750L']
		ngroups = 1

	# Create associations for each group of dithers, save ASN files, and run calwebb_spec3
	if params.instrument == 'nirspec':
		if params.obs_type == 'ifu':
			rule = Asn_Lv3NRSIFU
		else:
			rule = Asn_Lv3NRSFSS
	elif params.instrument == 'miri':
		if params.obs_type == 'ifu':
			rule = Asn_Lv3MIRMRS
		else:
			rule = Asn_Lv3SpectralTarget
	outdir = f'{params.data_dir}{params.prog_id}/Obs{params.obs_numb}/Stage3{params.stage3_suffix}/'
	os.makedirs(f'{outdir}', exist_ok=True)

	print(f'Stage 3: Running calwebb_spec3 pipeline...')
	# Custom step rules for various input types and specifications
	if params.obs_type == 'ifu':
		# Turn autocentroiding on
		try:
			params.stage3_rules['extract_1d'].update({'ifu_autocen' : True})
		except:
			params.stage3_rules['extract_1d'] = {'ifu_autocen' : True}

		# Set other parameters using pipeline settings
		if params.instrument == 'miri':
			params.stage3_rules['extract_1d'].update({'ifu_rfcorr' : True}) 		# Turn on extra spectrum-level defringing step
			params.stage3_rules['cube_build'] = {'output_type' : 'band'}			# Create separate Level 3 product for each MIRI sub-band (due to rotational modulation)
			params.stage3_rules['spectral_leak'] = {'skip' : False}					# Turn on optional spectral leak correction step
		if hasattr(params, 'cube_align'):
			if params.cube_align == 'ifu':
				cube_align_rule = {'coord_system' : 'ifualign'}
			elif params.cube_align == 'internal':
				cube_align_rule = {'coord_system' : 'internal_cal'}
			try:
				params.stage3_rules['cube_build'].update(cube_align_rule)
			except:
				params.stage3_rules['cube_build'] = cube_align_rule

	for i in range(ngroups):
		print(f'Combining dithers for {params.instrument.upper()} {used_gratings[i]}...')
		
		# Make association
		asn_name = f'jw{params.prog_id}-o{params.obs_numb}_t001_{params.instrument}_{used_gratings[i].lower()}'
		asn = asn_from_list.asn_from_list(sci_groups[i], rule=rule, format='json', product_name=asn_name)
		asn['asn_type'] = "spec3"
		asn['program'] = params.prog_id
		asn['asn_id'] = f'o{params.obs_numb}'
		asn['target'] = fits.getval(input_files[0], 'TARGNAME')

		# Save ASN file	
		serialized = association_io.json.dump(asn)[1]
		output_file = f'{outdir}{asn_name}_asn.json'
		f = open(output_file, 'w')
		f.write(serialized)
		f.close()

		# Process ASN file through jwst pipeline module calwebb_spec3
		Spec3Pipeline.call(output_file, output_dir=outdir, save_results=True, steps=params.stage3_rules, logcfg=cfg_file)

	print('Stage 3 processing complete!')

	return params