from . import aux
import os
import subprocess

def run(params):
	'''
	This function calls jwst_mast_query from command line and downloads the necessary JWST data from MAST

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
		# Create default directory structure, if needed
		files_exist = os.path.exists(f'{params.data_dir}{params.prog_id}/Obs{oo}/uncal/download_log.txt')	# Download record
		aux.make_dirs(params.data_dir, params.prog_id, oo, params.obs_type, params.dwnld_all)

		# Download data from MAST, if needed
		if params.download or not files_exist:
			print('Stage 0: Downloading data from MAST...')
			cmd_list = ["jwst_download.py","--propID",f"{params.prog_id}","--obsnums",f"{oo}","--instrument",f"{params.instrument}",
						"--obsmode",f"{params.obs_type}","--lookbacktime","100000","--clobber","--filetypes","uncal"] 
			if params.dwnld_all:
				cmd_list.extend(["rate","rateints","cal","calints","x1d","x1dints","s2d","s3d"])
			subprocess.run(cmd_list, input='y', text=True)

			# Move relevant data products to corresponding directories within data_dir
			print('Moving files to default locations...')
			aux.move_files(params.dwnld_dir, params.data_dir, params.prog_id, oo, params.obs_type, params.dwnld_all)

			print('Stage 0 data ingest complete!')
		else:
			print('Skipping data download...')

	return params