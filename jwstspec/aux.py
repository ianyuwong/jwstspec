import numpy as np
import astropy.io.fits as fits
from astropy.stats import SigmaClip, sigma_clip
import matplotlib as mpl 
mpl.use('tkagg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import Circle
from matplotlib.widgets import RectangleSelector
from matplotlib.ticker import MultipleLocator
from scipy.optimize import curve_fit, leastsq
from scipy.interpolate import interp1d, splrep, BSpline
from astropy.visualization import imshow_norm, LogStretch
import pickle
import copy
import os
import datetime
import warnings
import photutils

#Suppress non-critical warnings
warnings.filterwarnings('ignore', category=UserWarning, append=True)


#============================================================================================================
#%%%%%%%%%% OBJECT CLASSES %%%%%%%%%%
#============================================================================================================


class ifu_cube(object):
	'''
	IFU data cube object
	'''

	def __init__(self, filename, suffix=''):

		# Populate the ifu_cube object with data and header information
		self.longname = filename
		self.filename = filename.split('/')[-1]
		hdulist = fits.open(filename)
		self.hdulist = copy.deepcopy(hdulist)
		hdulist.close()

		self.meta = self.hdulist['PRIMARY',1].header			# Primary header metadata
		self.inst = self.meta['INSTRUME']
		self.detector = self.meta['DETECTOR']
		self.header = self.hdulist['SCI',1].header	# Header of SCI array
		self.imgs = self.hdulist['SCI',1].data		# Flux values
		self.errs = self.hdulist['ERR',1].data		# Flux errors
		self.dqs = self.hdulist['DQ',1].data 		# Data quality flags
		self.nimgs = len(self.imgs)
		self.pxscale = self.header['CDELT1']	# Pixel scale (in degrees)
		self.nx = self.header['NAXIS1']
		self.ny = self.header['NAXIS2']
		self.program = self.meta['PROGRAM']
		self.obs = self.meta['OBSERVTN']
		self.fluxunit = self.header['BUNIT']
		if self.inst == 'NIRSPEC':
			self.grating = self.meta['GRATING']
		elif self.inst == 'MIRI':
			if self.meta['BAND'] == 'MULTIPLE':									# Stage3 processing combining sub-bands
				self.grating = f'ch{self.meta["CHANNEL"]}'						# Use channel number instead, to prevent duplicate names
			elif self.meta['BAND'] in ['SHORT','MEDIUM','LONG']:				# Stage3 processing using separate files for each sub-band
				self.grating = f'ch{self.meta["CHANNEL"]}{self.meta["BAND"]}'	# Use channel number instead and band name
			else:
				self.grating = self.meta['BAND']
		if self.fluxunit == 'MJy/sr':
			self.conv_factor = self.header['PIXAR_SR']	# Convert to MJy
			self.imgs *= self.conv_factor
			self.errs *= self.conv_factor
			self.header['PHOTMJSR'] = 1.0 			# Fix to 1.0 for all IFU observations
													#	Some observations had wrong values saved

		self.visit = self.meta['VISIT']
		self.exposure = self.meta['EXPOSURE'].zfill(2)
		self.targ = self.meta['TARGNAME']
		self.targ_ra = self.meta['TARG_RA']
		self.targ_dec = self.meta['TARG_DEC']
		self.suffix = suffix
		self.outputbase = f'{filename.split("Obs")[0]}Obs{self.obs}/outputs/'

		# Mask out nonzero data quality pixels and off-sky pixels
		self.cube_mask = (self.dqs > 0) | (np.isnan(self.imgs) == True)
				
		# Save masked cube
		self.mask_imgs = np.ma.array(self.imgs, mask=self.cube_mask)
		self.mask_errs = np.ma.array(self.errs, mask=self.cube_mask)

		# Determine wavelength array (in microns)
		try:
			# Use SCI extension WCS
			self.wave = self.header['CRVAL3'] + (np.arange(self.header['NAXIS3']) - self.header['CRPIX3'] + 1) * self.header['CDELT3']
		except:
			# Determine wavelength array from WCS-TABLE extension (MIRI MRS data)
			wcs_table = self.hdulist['WCS-TABLE',1].data
			self.wave = wcs_table.field('wavelength')[0][:,0]


	def extract_spec(self, extr_aper_rad, bkg_aper_in, bkg_aper_out, fix_centroid=None, bkg_sig_clip=5.0, pix_sig_clip=8.0, extr_method='PSF', window_width=10, save_cleaned=False):
		'''
		Determines object centroid location and extracts spectrum

		Parameters
	   	----------
	    extr_aper_rad : float or int
	    	Radius of circular extraction aperture, in pixels,
	      	centered on the computed centroid position.
	    bkg_aper_in : float or int
	      	Inner radius, in pixels, of background extraction region for aperture photometry.
	    bkg_aper_out : float or int
	      	Outer radius, in pixels, of background extraction region for aperture photometry.
	      	Should be at least a few pixels larger than bkg_aper_in.
	    fix_centroid : list; optional
	    	User-specified centroid position to be used for extraction.
	    	Default is None.
	   	bkg_sig_clip : float or int; optional
	      	Threshold of sigma clipping for flagging background outliers.
	      	Default is 5.0.
	   	pix_sig_clip : float or int; optional
	      	Threshold of sigma clipping for flagging outliers in individual pixel arrays and for PSF fitting.
	      	Default is 8.0.
	    extr_method : str; optional
		   	Method for extracting flux in IFU cubes.
		window_width : int; optional
			Half-width of sliding median window for template PSF fitting.
			Default is 10 slices.
		save_cleaned : bool; optional
			Save outlier-cleaned data cube as a separate s3d file.
			Default is False.
	   	'''
		
		print(f"Working on {self.filename}...")

		# Save extraction parameters in object
		self.extr_method = extr_method
		self.pix_sig_clip = pix_sig_clip
		self.extr_aper_rad = extr_aper_rad
		self.bkg_aper_in = bkg_aper_in
		self.bkg_aper_out = bkg_aper_out
		self.fix_centroid = fix_centroid
		self.bkg_sig_clip = bkg_sig_clip
		self.window_width = window_width
		self.save_cleaned = save_cleaned
			
		# Define output directory
		if self.extr_method == 'aperture':
			self.outputdir = f'{self.outputbase}aper_{self.extr_aper_rad}_{self.bkg_aper_in}_{self.bkg_aper_out}{self.suffix}/'
		elif self.extr_method == 'PSF':
			self.outputdir = f'{self.outputbase}PSF_{self.extr_aper_rad}_{self.bkg_aper_in}{self.suffix}/'
		elif self.extr_method == 'box':
			self.outputdir = f'{self.outputbase}box_{self.extr_aper_rad}_{self.bkg_aper_in}{self.suffix}/'
		os.makedirs(self.outputdir, exist_ok=True)

		# Collapse IFU cube along wavelength axis
		self.master_frame = np.ma.median(self.mask_imgs, axis=0)
		self.master_err = np.ma.median(self.mask_errs, axis=0)

		# Use specified centroid position, if supplied
		if self.fix_centroid is not None and (type(self.fix_centroid[0]) is not str):
			self.centroid = self.fix_centroid

		# Otherwise, compute global flux centroid iteratively, using maximum flux location as initial guess
		else:
			# Restrict to central region due to high thermal backgrounds (MIRI) or FOV edge artifacts (NIRSpec)
			trimmed_master = copy.deepcopy(self.master_frame)
			new_mask = np.ones(np.shape(trimmed_master)).astype('bool')
			new_mask[10:35,10:35] = False
			trimmed_master.mask = np.logical_or(trimmed_master.mask, new_mask)
			self.centr_guess = np.unravel_index(trimmed_master.argmax(),trimmed_master.shape)[::-1]
			position = self.centr_guess
			
			iter = 0
			test_radius = 3		# For MIRI, you need small aperture if source is near the edge...
			while iter < 3:
				aper = photutils.aperture.CircularAperture(position, r=test_radius)
				stats = photutils.aperture.ApertureStats(self.master_frame, aper)
				position = stats.centroid
				iter += 1

			# Use 2D Gaussian fit of local region for a more accurate centroid position
			box_range = int(test_radius)	# Half-length of subarray sides
			position = position.astype('int')[::-1]
			subarray = self.master_frame[position[0]-box_range:position[0]+box_range+1,position[1]-box_range:position[1]+box_range+1]
			self.centroid = photutils.centroids.centroid_2dg(subarray) + position[::-1] - box_range

		print("Global centroid: x=","{:.3f}".format(self.centroid[0]),", y=","{:.3f}".format(self.centroid[1]))

		# Save integer-rounded centroid x and y positions, for use in box and PSF extraction methods
		self.centroidx = round(self.centroid[0])
		self.centroidy = round(self.centroid[1])

		# Flag outliers in slices and save cleaned data cube
		self.outlier_flag()

		# Shift centroid if relative offsets given
		if self.fix_centroid is not None and type(self.fix_centroid[0]) is str:
			self.centroid = self.centroid + np.asarray(self.fix_centroid).astype('float')
			self.centroidx = round(self.centroidx + float(self.fix_centroid[0]))
			self.centroidy = round(self.centroidy + float(self.fix_centroid[1]))
			print("Offset centroid coordinates given. Extracting spectrum using x=","{:.3f}".format(self.centroid[0]),", y=","{:.3f}".format(self.centroid[1]))

		# Extract spectrum
		if self.extr_method == 'aperture':
			self.aper_phot()
		elif self.extr_method == 'PSF':
			self.PSF_phot()
		elif self.extr_method == 'box':
			self.box_phot()
		print("Spectrum extracted!")	


	def outlier_flag(self):
		'''
		Automatically masks individual outliers pixel by pixel and saves cleaned data cube
		'''

		### Step 1: Slice-by-slice outlier masking in background region

		# Get mask for region inside background region
		inner_mask = np.zeros((self.ny, self.nx)).astype('bool')
		if self.extr_method == 'aperture':
			X, Y = np.meshgrid(np.arange(self.nx), np.arange(self.ny))
			diff = np.sqrt((X - self.centroidx)**2 + (Y - self.centroidy)**2)
			inner_mask[diff < self.bkg_aper_in] = True
		else:
			inner_mask[max(0,self.centroidy-self.bkg_aper_in):self.centroidy+self.bkg_aper_in+1, max(0,self.centroidx-self.bkg_aper_in):self.centroidx+self.bkg_aper_in+1] = True

		# Go slice by slice and flag outliers
		for i in range(self.nimgs):
			mask = np.logical_or(self.cube_mask[i], inner_mask)
			img = np.ma.MaskedArray(copy.deepcopy(self.imgs[i]), mask=mask)

			# Use sigma-clipping to find outliers
			clipped = sigma_clip(img, sigma=self.bkg_sig_clip, maxiters=None)

			# Isolate identified outliers
			outlier_mask = clipped.mask != mask
			self.dqs[i][outlier_mask] = 1000

			# Merge outlier mask with original slice mask
			self.cube_mask[i] = np.logical_or(self.cube_mask[i], outlier_mask)

		### Step 2: Pixel-by-pixel outlier masking in central region

		# Get x and y indices of pixels that lie in inner region
		yarr, xarr = np.where(inner_mask)

		# Go through pixel arrays one by one and flag outliers
		for i in range(len(xarr)):
			px_array = copy.deepcopy(self.mask_imgs[:,yarr[i],xarr[i]])
			px_array_err = copy.deepcopy(self.mask_errs[:,yarr[i],xarr[i]])
			ind = np.arange(len(px_array))
			med_err = np.ma.median(px_array_err / px_array)

			# Clip masked values and outliers in relative fluxerr
			ww = np.where((px_array.mask == False) & (px_array_err / px_array < 10 * med_err))
			clip = px_array[ww]
			clip_err = px_array_err[ww]
			clip_ind = ind[ww]
			
			if len(clip_ind) < 10:
				# Skip if pixel is at/beyond edge of field of view
				continue
			else:
				# Smooth pixel array and calculate residuals
				tck = splrep(clip_ind, clip, w=1/clip_err, s=2*len(clip_ind))	# The smoothing threshold is set empirically to handle most datasets
				diff = px_array - BSpline(*tck)(ind)

				# Find and flag outliers in both flux and fluxerr, update masked cube
				outlier_mask = (abs(diff) > self.pix_sig_clip * px_array_err) | (px_array_err > 100 * med_err)
				new_mask = np.logical_or(px_array.mask, outlier_mask)
				self.cube_mask[:,yarr[i],xarr[i]] = new_mask

				# plt.plot(ind,px_array,'r.')
				# plt.plot(ind[~new_mask],px_array[~new_mask],'g.')
				# plt.plot(clip_ind,BSpline(*tck)(clip_ind),'b-')
				# plt.savefig(f'{xarr[i]}-{yarr[i]}.png')
				# plt.close()

		# Update masked image and error cubes
		self.mask_imgs = np.ma.array(self.imgs, mask=self.cube_mask)
		self.mask_errs = np.ma.array(self.errs, mask=self.cube_mask)

		# Save cleaned cube, if needed
		if self.save_cleaned:
			outflux = copy.deepcopy(self.imgs) / self.conv_factor
			outerr = copy.deepcopy(self.errs) / self.conv_factor
			outflux[self.cube_mask] = np.nan
			outerr[self.cube_mask] = np.nan

			hdulist = copy.deepcopy(self.hdulist)
			hdulist['SCI',1].data = outflux
			hdulist['ERR',1].data = outerr
			hdulist['DQ',1].data = self.dqs
			hdulist.writeto(self.longname.replace('_s3d','_s3d_cleaned'), overwrite=True)
			hdulist.close()


	def box_phot(self):
		'''
		Carries out flux extraction from a box with background subtraction.
		'''

		# Extract flux, background from box centered on rounded centroid pixel
		flux = np.zeros(self.nimgs)
		flux_bkg = np.zeros(self.nimgs)
		fluxerr = np.zeros(self.nimgs)
		bkg = np.zeros(self.nimgs)
		mask_aper = np.zeros(self.nimgs)		# Mask-in-aperture flag array
		mask_close = np.zeros(self.nimgs)		# Mask-close-to-center flag array
		for i in range(self.nimgs):
			# Skip slice if it is all zeros or masked values (NaNs)
			if np.ma.sum(self.mask_imgs[i]) == 0 or np.sum(self.mask_imgs[i].mask) == self.nx*self.ny:
				continue

			# Adjust flux scaling
			rescale = 1e12 						# Convert to uJy for more manageable values
			img = self.mask_imgs[i] * rescale
			err = self.mask_errs[i] * rescale

			# Calculate background and mask out-of-box pixels
			box_width = self.extr_aper_rad
			box_width2 = self.bkg_aper_in
			sub_cube = img
			sub_cube2 = copy.deepcopy(sub_cube)
			sub_cube2.mask[max(0,self.centroidy-box_width2):self.centroidy+box_width2+1, max(0,self.centroidx-box_width2):self.centroidx+box_width2+1] = True
			master_subbkg = (sub_cube - np.ma.median(sub_cube2)).data
			mask = np.zeros(master_subbkg.shape).astype('bool')
			mask[max(0,self.centroidy-box_width):self.centroidy+box_width+1, max(0,self.centroidx-box_width):self.centroidx+box_width+1] = True
			master_subbkg[mask == False] = 0

			# Populate arrays
			bkg[i] = np.ma.median(sub_cube2)
			flux[i] = np.sum(master_subbkg)
			flux_bkg[i] = np.ma.sum(sub_cube[mask == True])
			fluxerr[i] = np.sqrt(np.ma.sum(err[mask == True] ** 2))

			# Flag wavelength slice if there is at least one masked point within fitting aperture
			aper_mask_sum = np.sum(img.mask[self.centroidy-box_width:self.centroidy+box_width+1, self.centroidx-box_width:self.centroidx+box_width+1])
			if aper_mask_sum > 0:
				mask_aper[i] = 1

			# Flag wavelength slice if there is at least one masked point within 7x7 box around centroid
			close_mask_sum = np.sum(img.mask[self.centroidy-3:self.centroidy+4, self.centroidx-3:self.centroidx+4])
			if close_mask_sum > 0:
				mask_close[i] = 1

		# Create spectrum object and save it
		spec = spectrum(self, flux_bkg/rescale, flux/rescale, fluxerr/rescale, None, None, None, bkg/rescale, None, mask_aper, mask_close)
		spec.save_spec()
			

	def PSF_phot(self):
		'''
		Carries out photometric extraction using empirical PSF fitting.
		'''

		print('Starting PSF fitting...')

		# Initiate data arrays
		bkg_per_pixel = np.zeros(self.nimgs)	# Background value per pixel
		bkg_aper = np.zeros(self.nimgs)			# Total background in aperture
		flux_aper_bkg = np.zeros(self.nimgs)	# Flux with no background subtraction
		flux_aper = np.zeros(self.nimgs)		# Background-subtracted flux
		fluxerr_aper = np.zeros(self.nimgs)		# Flux error
		mask_aper = np.zeros(self.nimgs)		# Mask-in-aperture flag array
		mask_close = np.zeros(self.nimgs)		# Mask-close-to-center flag array

		# Carry out iterative PSF fitting with outlier flagging for every wavelength slice one by one
		for i in range(self.window_width,self.nimgs-self.window_width-1):
			# Skip slice if it is all zeros or masked values (NaNs)
			if np.ma.sum(self.mask_imgs[i]) == 0 or np.sum(self.mask_imgs[i].mask) == self.nx*self.ny:
				continue

			# Adjust flux scaling
			rescale = 1e12 						# Convert to uJy for more manageable values
			img = copy.deepcopy(self.mask_imgs[i]) * rescale
			err = copy.deepcopy(self.mask_errs[i]) * rescale

			# Create template PSF
			box_width = int(np.round(self.extr_aper_rad))
			box_width2 = self.bkg_aper_in
			sub_cube = np.ma.median(self.mask_imgs[i-self.window_width:i+self.window_width+1,:,:], axis=0)
			sub_cube2 = copy.deepcopy(sub_cube)
			sub_cube2.mask[max(0,self.centroidy-box_width2):self.centroidy+box_width2+1, max(0,self.centroidx-box_width2):self.centroidx+box_width2+1] = True
			master_subbkg = (sub_cube - np.ma.median(sub_cube2)).data
			mask = np.zeros(master_subbkg.shape).astype('bool')
			mask[max(0,self.centroidy-box_width):self.centroidy+box_width+1, max(0,self.centroidx-box_width):self.centroidx+box_width+1] = True
			master_subbkg[mask == False] = 0
			self.psf_int = master_subbkg / np.ma.sum(master_subbkg)			

			# Define fitting region mask
			fit_mask = np.zeros(img.shape).astype('bool')
			fit_mask[max(0,self.centroidy-box_width2):self.centroidy+box_width2+1, max(0,self.centroidx-box_width2):self.centroidx+box_width2+1] = True
			fit_mask[max(0,self.centroidy-box_width):self.centroidy+box_width+1, max(0,self.centroidx-box_width):self.centroidx+box_width+1] = False

			# Iterative least-squares fitting
			iter = 0
			while iter < 3:
				try:
					res = leastsq(self.resid_eval, [np.ma.max(img),np.ma.median(img)], args=(img, err, fit_mask), full_output=True)
					model, model_bkg = self.model_eval(res[0], fit_mask)
					ratio = (img - model_bkg) / err

					abs_ratio = np.absolute(ratio)
					img.mask[abs_ratio > self.pix_sig_clip] = True 	# Mask out outliers
					iter += 1
				except:
					res = [None, None]
					iter = 3
				
			# Skip slice if optimization fails
			if res[1] is None:	
				print(f'Slice {i+1} of {self.nimgs}, {"{:.4f}".format(self.wave[i])} nm')	
				print(f'Skipping...PSF fitting did not converge!')
				continue

			# flux_aper[i] = np.sum(model)												# Integrated flux of best-fit PSF model
			# flux_aper_bkg[i] = np.sum(model_bkg)										# Integrated PSF flux with background level added (not really used)	
			# bkg_per_pixel[i] = res[0][1]												# Best-fit background level

			# Convert background level to counts for shot noise calculation
			bkg_per_pixel[i] = res[0][1]												# Best-fit background level
			uJy_per_cts	= self.header['PHOTMJSR'] * self.header['PIXAR_SR'] / self.meta['EFFEXPTM'] * rescale
			bkg_shot_noise = np.sqrt(abs(bkg_per_pixel[i]) * uJy_per_cts)				

			# Calculate flux uncertainty using weighted error array with background contribution
			# aper_mask = np.ones(img.shape).astype('bool')
			# aper_mask[self.centroidy-box_width:self.centroidy+box_width+1, self.centroidx-box_width:self.centroidx+box_width+1] = False
			# weight_err = np.ma.MaskedArray(np.sqrt((model/np.max(model) * err)**2 + bkg_shot_noise**2), mask=aper_mask)
			weight_err = np.ma.MaskedArray(np.sqrt((model/np.max(model) * err)**2 + bkg_shot_noise**2))

			# Define extraction aperture and carry out circular aperture extraction on the best-fit PSF model
			extr_aper = photutils.aperture.CircularAperture(self.centroid, r=self.extr_aper_rad)
			aper_stats = photutils.aperture.ApertureStats(model, extr_aper, error=weight_err)
			aper_stats_bkg = photutils.aperture.ApertureStats(model_bkg, extr_aper, error=weight_err)
			flux_aper[i] = aper_stats.sum
			flux_aper_bkg[i] = aper_stats_bkg.sum
			bkg_aper[i] = bkg_per_pixel[i] * aper_stats.sum_aper_area.value

			# Between the optimization flux uncertainty estimate and the total weighted error, choose the larger one
			flux_choices = [np.sqrt(res[1][0][0]), aper_stats.sum_err]	
			fluxerr_aper[i] = max(flux_choices)										# Total flux error

			print(f'Slice {i+1} of {self.nimgs}, {"{:.4f}".format(self.wave[i])} nm')
			print(f'Flux = {"{:.4f}".format(flux_aper[i])}')
			print(f'Fluxerr = {"{:.4f}".format(fluxerr_aper[i])}')
			print(f'Bkg = {"{:.4f}".format(bkg_aper[i])}')

			# Flag wavelength slice if there is at least one masked point within fitting aperture
			aper_mask_sum = np.sum(img.mask[self.centroidy-box_width:self.centroidy+box_width+1, self.centroidx-box_width:self.centroidx+box_width+1])
			if aper_mask_sum > 0:
				mask_aper[i] = 1

			# Flag wavelength slice if there is at least one masked point within 7x7 box around centroid
			close_mask_sum = np.sum(img.mask[self.centroidy-3:self.centroidy+4, self.centroidx-3:self.centroidx+4])
			if close_mask_sum > 0:
				mask_close[i] = 1

		# Create spectrum object and save it
		spec = spectrum(self, flux_aper_bkg/rescale, flux_aper/rescale, fluxerr_aper/rescale, None, None, None, bkg_per_pixel/rescale, bkg_aper/rescale, mask_aper, mask_close)
		spec.save_spec()


	def aper_phot(self):
		'''
		Carries out aperture extraction and background subtraction 
		'''

		# Initiate data arrays
		cent_x = np.zeros(self.nimgs)			# Centroid x-position
		cent_y = np.zeros(self.nimgs)			# Centroid y-position
		fwhm = np.zeros(self.nimgs)				# Source FWHM (in pixels)
		bkg_per_pixel = np.zeros(self.nimgs)	# Background value per pixel
		bkg_aper = np.zeros(self.nimgs)			# Total background in aperture
		flux_aper_bkg = np.zeros(self.nimgs)	# Flux with no background subtraction
		flux_aper = np.zeros(self.nimgs)		# Background-subtracted flux
		fluxerr_aper = np.zeros(self.nimgs)		# Flux error
		mask_aper = np.zeros(self.nimgs)		# Mask-in-aperture flag array
		mask_close = np.zeros(self.nimgs)		# Mask-close-to-center flag array

		# Remove NaNs from flux arrays
		self.imgs[np.isnan(self.imgs)] = 0.0
		self.errs[np.isnan(self.errs)] = 0.0

		# Cycle through wavelength slices to calculate background, and extract flux
		for i in range(self.nimgs):

			# Skip slice if it is all zeros or masked values (NaNs)
			if np.ma.sum(self.mask_imgs[i]) == 0 or np.sum(self.mask_imgs[i].mask) == self.nx*self.ny:
				continue

			img = self.imgs[i]
			err = self.errs[i]

			# Mask flagged values
			mask = self.cube_mask[i]

			# Define apertures and extract aperture information
			sigclip = SigmaClip(sigma=self.bkg_sig_clip, maxiters=None)
			extr_aper = photutils.aperture.CircularAperture(self.centroid, r=self.extr_aper_rad)
			close_aper = photutils.aperture.CircularAperture(self.centroid, r=3)
			back_aper = photutils.aperture.CircularAnnulus(self.centroid, r_in=self.bkg_aper_in, r_out=self.bkg_aper_out)
			aper_stats = photutils.aperture.ApertureStats(img, extr_aper, error=err, mask=mask, sigma_clip=None)
			bkg_stats = photutils.aperture.ApertureStats(img, back_aper, mask=mask, sigma_clip=sigclip)
			mask_stats = photutils.aperture.ApertureStats(mask, extr_aper)
			mask_stats_close = photutils.aperture.ApertureStats(mask, close_aper)	# Check for masked pixels near center of PSF

			# Compute background
			bkg_per_pixel[i] = bkg_stats.median	
			area = aper_stats.sum_aper_area.value	# Total extraction aperture area
			bkg_aper[i] = bkg_per_pixel[i] * area	

			# Extract background-subtracted flux
			flux_aper_bkg[i] = aper_stats.sum
			flux_aper[i] = flux_aper_bkg[i] - bkg_aper[i]
			fluxerr_aper[i] = aper_stats.sum_err

			# Save per-slice centroid and FWHM information (skip if aperture stats do not work)
			img_subbkg = img - bkg_stats.median
			aper_stats_subbkg = photutils.aperture.ApertureStats(img_subbkg, extr_aper, error=err, mask=mask, sigma_clip=None)
			position = aper_stats_subbkg.centroid
			if position[0] < 0 or position[0] > 50 or position[1] < 0 or position[1] > 50:
				continue
			else:
				cent_x[i], cent_y[i] = position[0], position[1]
				fwhm[i] = aper_stats_subbkg.fwhm.value
		
			# Flag wavelength slice if there is at least one masked point within aperture
			if mask_stats.sum > 0:
				mask_aper[i] = 1

			# Flag wavelength slice if there is at least one masked point within 3 pixels of centroid
			if mask_stats_close.sum > 0:
				mask_close[i] = 1

		# Create spectrum object and save it
		os.makedirs(self.outputdir, exist_ok=True)
		spec = spectrum(self, flux_aper_bkg, flux_aper, fluxerr_aper, cent_x, cent_y, fwhm, bkg_per_pixel, bkg_aper, mask_aper, mask_close)
		spec.save_spec()


	def resid_eval(self, arr, img, err, mask):
		'''
		Function that returns flattened error-normalized residual array
		'''

		model = arr[0] * self.psf_int
		modelsum_bkg = np.ma.MaskedArray(model + arr[1], mask=mask)
		resid = (img - modelsum_bkg) / err

		return resid.compressed()


	def model_eval(self, arr, mask):
		'''
		Function that evalutes the PSF model and returns the flux with and without background
		'''

		model = np.ma.MaskedArray(arr[0] * self.psf_int, mask=mask)
		modelsum_bkg = model + arr[1]

		return model, modelsum_bkg


	def master_plot(self):
		'''
		Plots the master frame with annotated centroid and apertures for circular aperture extraction
		'''

		outf = f'{self.outputdir}Visit{self.visit}_{self.grating}_{self.detector}_exp{self.exposure}_masterframe.pdf'

		# Flatten low outliers
		maskimg = copy.deepcopy(self.master_frame)
		threshold = np.nanpercentile(maskimg.filled(np.nan), 5)
		maskimg[np.ma.where(maskimg<threshold)] = threshold

		fig = plt.figure(figsize=(5,5))
		ax = fig.add_subplot(111)
		ax.imshow(maskimg,norm=colors.PowerNorm(0.2),cmap='viridis',origin='lower')
		ax.plot(self.centroid[0],self.centroid[1],'k.',ms=6)
		circ1 = Circle(xy=(self.centroid[0],self.centroid[1]),radius=self.bkg_aper_in,facecolor='none',edgecolor='red',linestyle='-')
		ax.add_artist(circ1)
		circ2 = Circle(xy=(self.centroid[0],self.centroid[1]),radius=self.bkg_aper_out,facecolor='none',edgecolor='red',linestyle='-')
		ax.add_artist(circ2)
		circ3 = Circle(xy=(self.centroid[0],self.centroid[1]),radius=self.extr_aper_rad,facecolor='none',edgecolor='black',linestyle='-')
		ax.add_artist(circ3)
		ax.xaxis.set_major_locator(MultipleLocator(10))
		ax.yaxis.set_major_locator(MultipleLocator(10))
		ax.xaxis.set_minor_locator(MultipleLocator(5))
		ax.yaxis.set_minor_locator(MultipleLocator(5))
		ax.tick_params(labelsize=12)
		ax.set_xlabel('x [px]',fontsize=14)
		ax.set_ylabel('y [px]',fontsize=14)

		plt.savefig(outf)
		plt.close()

#============================================================================================================

class slit_spectrum(object):
	'''
	Fixed-slit spectrum object
	'''

	def __init__(self, filename, suffix='', nod_subtract=False):

		# Populate the slit_spectrum object with data and header information
		self.filename = filename.split('/')[-1]
		hdulist = fits.open(filename)
		self.meta = hdulist['PRIMARY',1].header			# Primary header metadata
		self.inst = self.meta['INSTRUME']
		self.detector = self.meta['DETECTOR']
		self.header = hdulist['SCI',1].header	# Header of SCI array
		self.img = hdulist['SCI',1].data		# Flux values
		self.err = hdulist['ERR',1].data		# Flux errors
		self.program = self.meta['PROGRAM']
		self.obs = self.meta['OBSERVTN']
		self.fluxunit = self.header['BUNIT']
		if self.inst == 'NIRSPEC':
			self.grating = self.meta['GRATING']
		elif self.inst == 'MIRI':
			self.grating = self.meta['FILTER']
		if self.fluxunit == 'MJy/sr':
			self.img *= self.header['PIXAR_SR']	# Convert to MJy
			self.err *= self.header['PIXAR_SR']
		self.visit = self.meta['VISIT']
		self.exposure = self.meta['EXPOSURE'].zfill(2)
		self.targ = self.meta['TARGNAME']
		self.nx = self.header['NAXIS1']
		self.ny = self.header['NAXIS2']
		self.suffix = suffix
		self.outputbase = f'{filename.split("Obs")[0]}Obs{self.obs}/outputs/'

		# Establish orientation for instrument
		if self.inst == 'NIRSPEC':
			self.slit_direction = 'y'
			self.slit_length = self.ny
			self.spec_length = self.nx
		elif self.inst == 'MIRI':
			self.slit_direction = 'x'
			self.slit_length = self.nx
			self.spec_length = self.ny

		# Mask out bad pixels
		self.image_mask = (self.img == 0)|(np.isnan(self.img))

		# Exchange axes of MIRI LRS for ease in processing
		if self.inst == 'MIRI':
			self.img = self.img.T
			self.err = self.err.T

		# Mask negative trace in case of pairwise nod subtraction
		if nod_subtract:
			collapsed = np.ma.median(np.ma.array(self.img, mask=self.image_mask), axis=1)
			min_loc = collapsed.argmin()
			self.custom_mask = np.zeros(self.img.shape).astype('bool')
			self.custom_mask[min_loc-8:min_loc+7,:] = True
			if self.inst == 'MIRI':
				self.custom_mask[40:,:] = True
			self.image_mask = np.logical_or(self.image_mask, self.custom_mask)

		# Save masked image
		self.mask_img = np.ma.array(self.img, mask=self.image_mask)
		self.mask_err = np.ma.array(self.err, mask=self.image_mask)

		# Determine wavelength array (in microns) from corresponding x1d file (maybe there's an easier way?...)
		x1d_file = filename.replace('s2d','x1d')
		self.wave = fits.getdata(x1d_file, 'EXTRACT1D')['WAVELENGTH']


	def extract_spec(self, extr_aper_rad, bkg_aper_in, bkg_aper_out, fix_centroid=None, bkg_sig_clip=3.0, pix_sig_clip=8.0, extr_method='PSF', window_width=20, save_cleaned=False):
		'''
		Determines spectrum centroid location and extracts spectrum

		Parameters
		----------
		extr_aper_rad : float or int; optional
			Radius of circular extraction aperture, in pixels,
			centered on the computed centroid position.
		bkg_aper_in : float or int; optional
	    	Inner radius, in pixels, of background extraction region for aperture photometry.
		bkg_aper_out : float or int; optional
			Outer radius, in pixels, of background extraction region for aperture photometry.
			Should be at least a few pixels larger than bkg_aper_in.
		fix_centroid : float or int; optional
			User-specified centroid position to be used for extraction.
			Default is None.
		bkg_sig_clip : float or int; optional
			Threshold of sigma clipping for determining background level.
			Default is 3.0.
		pix_sig_clip : float or int; optional
		  	Threshold of sigma clipping for flagging outliers in individual pixel arrays and for PSF fitting.
		  	Default is 8.0.
		extr_method : str; optional
		   	Method for extracting flux in s2d files: 'PSF' or 'aperture'.
		window_width : int; optional
			Half-width of sliding median window for template PSF fitting.
			Default is 20 rows/columns in the dispersion direction.
		save_cleaned : bool; optional
			This option is not used for slit spectra.
		'''

		print(f"Working on {self.filename}...")

		# Save extraction parameters in object
		self.extr_aper_rad = extr_aper_rad
		self.bkg_aper_in = bkg_aper_in
		self.bkg_aper_out = bkg_aper_out
		self.fix_centroid = fix_centroid
		self.bkg_sig_clip = bkg_sig_clip
		self.pix_sig_clip = pix_sig_clip
		self.window_width = window_width

		# Define output directory
		if extr_method == 'aperture':
			self.outputdir = f'{self.outputbase}{self.extr_aper_rad}_{self.bkg_aper_in}_{self.bkg_aper_out}{self.suffix}/'
		elif extr_method == 'PSF':
			self.outputdir = f'{self.outputbase}PSF_{self.extr_aper_rad}_{self.bkg_aper_in}{self.suffix}/'
		os.makedirs(self.outputdir, exist_ok=True)

		# Collapse image along wavelength axis
		self.master_array = np.ma.median(self.mask_img, axis=1)

		# Use specified centroid position, if supplied
		if self.fix_centroid is not None:
			self.centroid = self.fix_centroid

		# Otherwise, compute global flux centroid
		else:
			def gauss(x, back, peak, centr, sigma):
				return back + peak * np.exp(-(x - centr) ** 2 / (2 * sigma ** 2))

			res = curve_fit(gauss, np.arange(self.slit_length), self.master_array, 
	    		p0=[np.min(self.master_array), np.max(self.master_array), self.master_array.argmax(), 3])[0]
			self.psf_model = gauss(np.arange(self.slit_length), res[0], res[1], res[2], res[3])
			self.centroid = res[2]
			# self.centroid = self.master_array.argmax()
		print(f"Global centroid: {self.slit_direction} = {'{:.3f}'.format(self.centroid)}")

		# Extract spectrum
		if extr_method == 'aperture':
			self.aper_phot()
		elif extr_method == 'PSF':
			self.PSF_phot()
		else:
			raise ValueError(f'Photometric extraction method {extr_method} not recognized!')
		print("Spectrum extracted!")


	def PSF_phot(self):
		'''
		Carries out photometric extraction using empirical PSF fitting.
		'''

		print('Starting PSF fitting...')

		# Carry out iterative PSF fitting with outlier flagging for every column one by one
		flux = np.zeros(self.spec_length)
		flux_bkg = np.zeros(self.spec_length)
		fluxerr = np.zeros(self.spec_length)
		bkg = np.zeros(self.spec_length)
		mask_aper = np.zeros(self.spec_length)		# Mask-in-aperture flag array
		mask_close = np.zeros(self.spec_length)		# Mask-close-to-center flag array
		for i in range(self.window_width,self.spec_length-self.window_width-1):
			# Skip column if it is all zeros or masked values (NaNs)
			if np.ma.sum(self.mask_img[:,i]) == 0 or np.sum(self.mask_img.mask[:,i]) == self.ny:
				continue

			# Adjust flux scaling
			rescale = 1e12 						# Convert to uJy for more manageable values
			col = self.mask_img[:,i] * rescale
			err = self.mask_err[:,i] * rescale

			# Create template PSF
			box_width = self.extr_aper_rad
			box_width2 = self.bkg_aper_in
			sub_col = np.ma.median(self.mask_img[:,i-self.window_width:i+self.window_width+1], axis=1)
			sub_col2 = copy.deepcopy(sub_col)
			sub_col2.mask[max(0,round(self.centroid)-box_width2):round(self.centroid)+box_width2+1] = True
			master_subbkg = (sub_col - np.ma.median(sub_col2)).data
			mask = np.zeros(master_subbkg.shape).astype('bool')
			mask[max(0,round(self.centroid)-box_width):round(self.centroid)+box_width+1] = True
			master_subbkg[mask == False] = 0
			self.psf_int = master_subbkg / np.sum(master_subbkg)

			# Define fitting region mask
			fit_mask = np.zeros(self.psf_int.shape).astype('bool')
			fit_mask[max(0,round(self.centroid)-box_width2):round(self.centroid)+box_width2+1] = True
			fit_mask[max(0,round(self.centroid)-box_width):round(self.centroid)+box_width+1] = False

			# Iterative least-squares fitting
			iter = 0
			while iter < 3:
				try:
					res = leastsq(self.resid_eval, [np.ma.max(col),np.ma.median(col)], args=(col, err, fit_mask), full_output=True)		
					model, model_bkg = self.model_eval(res[0], fit_mask)
					ratio = (col - model_bkg) / err

					abs_ratio = np.absolute(ratio)
					col.mask[abs_ratio > self.pix_sig_clip] = True 	# Mask out outliers
					iter += 1
				except:
					res = [None, None]
					iter = 3

			# Skip slice if optimization fails
			if res[1] is None:	
				print(f'Column {i+1} of {self.spec_length}, {"{:.4f}".format(self.wave[i])} nm')	
				print(f'Skipping...PSF fitting did not converge!')
				continue

			flux[i] = np.sum(model)												# Integrated flux of best-fit PSF model
			flux_bkg[i] = np.sum(model_bkg)										# Integrated PSF flux with background level added (not really used)
			bkg[i] = res[0][1]													# Best-fit background level
			
			# Convert background level to counts for shot noise calculation
			uJy_per_cts	= self.header['PHOTMJSR'] * self.header['PIXAR_SR'] / self.meta['EFFEXPTM'] * rescale
			bkg_shot_noise = np.sqrt(abs(bkg[i]) * uJy_per_cts)

			# Calculate flux uncertainty using weighted error array with background contribution
			aper_mask = np.ones(len(col)).astype('bool')
			aper_mask[max(0,round(self.centroid)-box_width):round(self.centroid)+box_width+1] = False
			weight_err = np.ma.MaskedArray(np.sqrt((model/np.max(model) * err)**2 + bkg_shot_noise**2), mask=aper_mask)

			# Between the optimization flux uncertainty estimate and the weighted error, choose the larger one
			flux_choices = [np.sqrt(res[1][0][0]), np.sqrt(np.ma.sum(weight_err**2))]	
			fluxerr[i] = max(flux_choices)										# Total flux error

			weight_err = np.sqrt((model/np.max(model) * err)**2 + bkg_shot_noise**2)		# Weighted error array with background contribution
			fluxerr[i] = np.sqrt(np.ma.sum(weight_err**2))						# Total flux error

			print(f'Column {i+1} of {self.spec_length}, {"{:.4f}".format(self.wave[i])} nm')
			print(f'Flux = {"{:.4f}".format(flux[i])}')
			print(f'Fluxerr = {"{:.4f}".format(fluxerr[i])}')
			print(f'Bkg = {"{:.4f}".format(bkg[i])}')

			# Flag wavelength slice if there is at least one masked point within fitting aperture
			aper_mask_sum = np.sum(col.mask[int(self.centroid)-box_width:int(self.centroid)+box_width+1])
			if aper_mask_sum > 0:
				mask_aper[i] = 1

			# Flag wavelength slice if there is at least one masked point within 3 pixels of centroid
			close_mask_sum = np.sum(col.mask[int(self.centroid)-3:int(self.centroid)+4])
			if close_mask_sum > 0:
				mask_close[i] = 1

		# Create spectrum object and save it
		spec = spectrum(self, flux_bkg/rescale, flux/rescale, fluxerr/rescale, None, None, None, bkg/rescale, None, mask_aper, mask_close)
		spec.save_spec()


	def aper_phot(self):
		'''
		Carries out aperture extraction and background subtraction 
	    '''

		# Establish extraction regions
		centr_int = round(self.centroid)
		self.aper_low = max(centr_int - self.extr_aper_rad, 0)
		self.aper_high = min(centr_int + self.extr_aper_rad + 1, self.slit_length)
		self.back_a1 = max(centr_int - self.bkg_aper_out, 0)
		self.back_a2 = max(centr_int - self.bkg_aper_in + 1, 0)
		self.back_b1 = min(centr_int + self.bkg_aper_in, self.slit_length)
		self.back_b2 = min(centr_int + self.bkg_aper_out + 1, self.slit_length)
		extr_region = self.mask_img[self.aper_low:self.aper_high]
		extr_region_err = self.mask_err[self.aper_low:self.aper_high]
		bkg_region = np.concatenate([self.mask_img[self.back_a1:self.back_a2], self.mask_img[self.back_b1:self.back_b2]])

		# Determine background
		clipped = sigma_clip(bkg_region, axis=0, sigma=self.bkg_sig_clip, maxiters=None)
		bkg_per_pixel = np.ma.median(clipped, axis=0)
		bkg_aper = bkg_per_pixel * (self.extr_aper_rad * 2 + 1)

		# Extract flux
		flux_aper_bkg = np.ma.sum(extr_region, axis=0)
		flux_aper = flux_aper_bkg - bkg_aper
		fluxerr_aper = np.ma.sqrt(np.ma.sum(extr_region_err**2, axis=0))

		# Create spectrum object and save it
		spec = spectrum(self, flux_aper_bkg, flux_aper, fluxerr_aper, None, None, None, bkg_per_pixel, bkg_aper, None, None)
		spec.save_spec()


	def resid_eval(self, arr, col, err, mask):
		'''
		Function that returns flattened error-normalized residual array
		'''

		model = arr[0] * self.psf_int
		modelsum_bkg = np.ma.MaskedArray(model + arr[1], mask=mask)
		resid = (col - modelsum_bkg) / err

		return resid


	def model_eval(self, arr, mask):
		'''
		Function that evalutes the PSF model and returns the flux with and without background
		'''

		model = arr[0] * self.psf_int
		modelsum_bkg = np.ma.MaskedArray(model + arr[1], mask=mask)

		return model, modelsum_bkg


	def master_plot(self):
		'''
		Plots the master array with annotated centroid and apertures
		'''

		outf = f'{self.outputdir}Visit{self.visit}_{self.grating}_{self.detector}_exp{self.exposure}_masterarray.pdf'

		fig = plt.figure(figsize=(5,5))
		ax = fig.add_subplot(111)
		ax.plot(np.arange(self.slit_length),self.master_array,'r-')
		if self.fix_centroid is None:
			ax.plot(np.arange(self.slit_length),self.psf_model,'b-.')
		ax.axvline(self.centroid,color='black',linestyle='--')
		ax.axvline(self.aper_low,color='black',linestyle='-')
		ax.axvline(self.aper_high-1,color='black',linestyle='-')
		ax.axvline(self.back_a1,color='black',linestyle=':')
		ax.axvline(self.back_a2-1,color='black',linestyle=':')
		ax.axvline(self.back_b1,color='black',linestyle=':')
		ax.axvline(self.back_b2-1,color='black',linestyle=':')
		ax.tick_params(labelsize=12)
		ax.set_xlim([0,self.slit_length])
		ax.set_xlabel(f'{self.slit_direction} [px]',fontsize=14)
		ax.set_ylabel('Flux [MJy]',fontsize=14)

		plt.savefig(outf)
		plt.close()

#============================================================================================================

class spectrum(object):
	'''
	Extracted spectrum object
	'''

	def __init__(self, input_obj, flux_aper_bkg, flux_aper, fluxerr_aper, cent_x, cent_y, fwhm, bkg_per_pixel, bkg_aper, mask_aper, mask_close):

		# Extract relevant information from input_obj and inputs
		self.filename = input_obj.filename
		self.wave = input_obj.wave
		self.fluxbkg = flux_aper_bkg
		self.flux = flux_aper 
		self.fluxerr = fluxerr_aper
		self.centroid = input_obj.centroid
		self.cent_x = cent_x
		self.cent_y = cent_y
		self.fwhm = fwhm
		self.bkg_per_pixel = bkg_per_pixel
		self.bkg = bkg_aper
		self.mask_aper = mask_aper
		self.mask_close = mask_close
		self.extr_aper_rad = input_obj.extr_aper_rad
		self.pix_sig_clip = input_obj.pix_sig_clip
		self.bkg_aper_in = input_obj.bkg_aper_in
		self.bkg_aper_out = input_obj.bkg_aper_out
		self.fix_centroid = input_obj.fix_centroid
		self.bkg_sig_clip = input_obj.bkg_sig_clip
		self.window_width = input_obj.window_width
		self.program = input_obj.program
		self.obs = input_obj.obs
		self.inst = input_obj.inst
		self.grating = input_obj.grating
		self.detector = input_obj.detector
		self.visit = input_obj.visit
		self.meta = input_obj.meta
		self.exposure = input_obj.exposure
		self.targ = input_obj.targ
		self.outputdir = input_obj.outputdir

	def save_spec(self):
		'''
		Saves spectrum object into pickle file
		'''

		outf = f'{self.outputdir}Visit{self.visit}_{self.grating}_{self.detector}_exp{self.exposure}_spectrum.pickle'
		f = open(outf,'wb')
		pickle.dump(self,f)
		f.close()

#============================================================================================================

class params(object):
	'''
	Pipeline run parameters object
	'''

	def __init__(self):

		self.stage1_rules = {}
		self.stage2_rules = {}
		self.stage3_rules = {}

	def add_params(self, prog_id, obs_numb, bkg_obs_numb, instrument, obs_type, tso_observation, dwnld_dir, data_dir, download, dwnld_all,
					readnoise_correct, bkg_subtract, cube_align, stage1_suffix, stage2_suffix, stage3_suffix, extract_stage, extr_method, extr_suffix, extr_aper_rad,
					bkg_aper_in, bkg_aper_out, window_width, pix_sig_clip, bkg_sig_clip, fix_centroid, save_cleaned, spec_bkg_sub, spec_sig_clip, spec_window_half, special_defringe):
		'''
		Add parameters. See run_jwstspec.py file for parameter definitions.
		'''

		self.prog_id = prog_id
		self.obs_numb = obs_numb
		self.bkg_obs_numb = bkg_obs_numb
		self.instrument = instrument
		self.obs_type = obs_type
		self.tso_observation = tso_observation
		self.dwnld_dir = dwnld_dir
		self.data_dir = data_dir
		self.download = download
		self.dwnld_all = dwnld_all
		self.readnoise_correct = readnoise_correct
		self.bkg_subtract = bkg_subtract
		self.cube_align = cube_align
		self.stage1_suffix = stage1_suffix
		self.stage2_suffix = stage2_suffix
		self.stage3_suffix = stage3_suffix
		self.extract_stage = extract_stage
		self.extr_method = extr_method
		self.extr_suffix = extr_suffix
		self.extr_aper_rad = extr_aper_rad
		self.bkg_aper_in = bkg_aper_in
		self.bkg_aper_out = bkg_aper_out
		self.window_width = window_width
		self.pix_sig_clip = pix_sig_clip
		self.bkg_sig_clip = bkg_sig_clip
		self.fix_centroid = fix_centroid
		self.save_cleaned = save_cleaned
		self.spec_bkg_sub = spec_bkg_sub
		self.spec_sig_clip = spec_sig_clip
		self.spec_window_half = spec_window_half
		self.special_defringe = special_defringe

		# Run sanity check on parameters
		self.sanity_check()


	def sanity_check(self):
		'''
		Do a quick check to make sure parameter settings aren't conflicting or unreasonable
		'''

		# Define affix for different data types and custom destriping methods
		vers = ''
		if self.tso_observation:
			vers += 'ints'
		if self.instrument == 'nirspec' and hasattr(self, 'readnoise_correct'):
			if self.readnoise_correct == 'nsclean':
				vers += 'corr0'
			elif self.readnoise_correct == 'constant':
				vers += 'corr1' 
			elif self.readnoise_correct == 'moving_median':
				vers += 'corr2'
		if self.bkg_subtract in ['pixel', 'asn']:
			vers += f'-subbkg_{self.bkg_subtract}'
		self.vers = vers

		# Check instrument and observation type settings
		if self.instrument not in ['nirspec', 'miri'] or self.obs_type not in ['ifu', 'slit', 'slitless']:
			raise ValueError(f'Instrument configuration {self.instrument.upper()} {self.obs_type.upper()} not supported.')

		# Check destriping method
		if self.instrument == 'nirspec':
			if self.readnoise_correct not in [None, 'nsclean', 'constant', 'moving_median']:
				raise ValueError(f'Destriping method {self.destriping} not recognized.')

		# Check bkg_subtract method
		if self.bkg_subtract is not None:
			if self.bkg_subtract not in ['asn', 'pixel']:
				raise ValueError(f'Background subtraction method {self.bkg_subtract} not recognized.')
			if self.bkg_obs_numb is None:
				raise ValueError(f'Must provide background observation number if bkg_subtract is specified.')

		if hasattr(self, 'extract_stage'):	
			# Check extract_stage
			if self.extract_stage not in ['Stage2', 'Stage3']:
				raise ValueError(f'params.extract_stage must be Stage2 or Stage3.')

			# For PSF fitting, warn user if window width is larger than recommended
			if self.extr_method == 'PSF' and self.obs_type == 'ifu' and ((self.instrument == 'nirspec' and self.window_width > 10) or (self.instrument == 'miri' and self.window_width > 5)):
				print(f'WARNING! window_width = {self.window_width} is larger than recommended settings (<=10 for NIRSpec, <=5 for MIRI)')
				go_ahead = input('Continue? [y/n]')
				if go_ahead == 'y':
					pass
				else:
					stop

			# For box and IFU extraction methods and all slit spectrum extractions, enforce integer aperture size requirement
			if self.extr_method == 'PSF':
				if not isinstance(self.bkg_aper_in, int):
					raise TypeError(f'Background aperture size for PSF extraction must be an integer.')
			if self.extr_method == 'box':
				if not isinstance(self.extr_aper_rad, int) or not isinstance(self.bkg_aper_in, int):
					raise TypeError(f'Aperture sizes for box extraction must be integers.')
			if self.obs_type == 'slit':
				if not isinstance(self.extr_aper_rad, int) or not isinstance(self.bkg_aper_in, int) or not isinstance(self.bkg_aper_out, int):
					raise TypeError(f'Aperture sizes for slit observations must be integers.')


#============================================================================================================
#%%%%%%%%%%%% OTHER FUNCTIONS %%%%%%%%%%%%%
#============================================================================================================


def make_dirs(data_dir, prog_id, obs_numb, obs_type, dwnld_all):
	'''
	Creates default directory structure in preparation for JWST data download from MAST

	Parameters
	----------
	data_dir : str
		Path to top-level data directory.
	prog_id : str
	   	Five-digit program ID.
	obs_numb : str
		Three-digit observation number.
	obs_type : str
	   	'ifu' or 'slit' or 'slitless'.
	dwnld_all : bool
		Choice to download all data types. If False, only uncal files are downloaded.
	'''

	if obs_type not in ['ifu','slit', 'slitless']:
		raise ValueError(f'Observation type {obs_type} not recognized!')

	os.makedirs(f'{data_dir}/{prog_id}/Obs{obs_numb}/uncal/', exist_ok=True)			# uncal files from MAST
	if dwnld_all:
		os.makedirs(f'{data_dir}/{prog_id}/Obs{obs_numb}/rate/', exist_ok=True)			# rate and rateints files from MAST (ramp-fitted frames)
		os.makedirs(f'{data_dir}/{prog_id}/Obs{obs_numb}/cal/', exist_ok=True)			# cal or calints files from MAST
		if obs_type == 'slit':
			os.makedirs(f'{data_dir}/{prog_id}/Obs{obs_numb}/s2d/', exist_ok=True)		# s2d files from MAST (slit spectroscopy)
		elif obs_type == 'ifu':
			os.makedirs(f'{data_dir}/{prog_id}/Obs{obs_numb}/s3d/', exist_ok=True)		# s3d files from MAST (IFU cubes)
		os.makedirs(f'{data_dir}/{prog_id}/Obs{obs_numb}/x1d/', exist_ok=True)			# x1d or x1dints files from MAST
	os.makedirs(f'{data_dir}/{prog_id}/Obs{obs_numb}/Stage1/', exist_ok=True)			# Output directory for jwst Stage 1 pipeline
	os.makedirs(f'{data_dir}/{prog_id}/Obs{obs_numb}/Stage2/', exist_ok=True)			# Output directory for jwst Stage 2 pipeline
	os.makedirs(f'{data_dir}/{prog_id}/Obs{obs_numb}/Stage3/', exist_ok=True)			# Output directory for jwst Stage 3 pipeline
	os.makedirs(f'{data_dir}/{prog_id}/Obs{obs_numb}/outputs/', exist_ok=True)			# Default results directory
	print('Data directories created.')


def move_files(dwnld_dir, data_dir, prog_id, obs_numb, obs_type, dwnld_all):
	'''
	Moves relevant data downloaded from MAST to appropriate directories within data_dir

	Parameters
	----------
	dwnld_dir : str
   		Path to jwst_mast_query download directory.
	data_dir : str
		Path to top-level data directory.
	prog_id : str
	   	Five-digit program ID.
	obs_numb : str
		Three-digit observation number.
	obs_type : str
	   	'ifu' or 'slit' or 'slitless'.
	dwnld_all : bool
		Choice to download all data types. If False, only uncal files are downloaded.
	'''

	# Account for obs_numb naming convention inconsistency in jwst_mast_query
	if obs_numb[0] == '0':
		limit = -2
	else:
		limit = -3
	prefix = f'{dwnld_dir}{prog_id}/obsnum{obs_numb[limit:]}/'

	# Moving uncal files (only science images, excluding target acquisition or direct images)
	cmd = f'mv {prefix}jw{prog_id}{obs_numb}???*uncal.fits {data_dir}{prog_id}/Obs{obs_numb}/uncal/'
	os.system(cmd)

	if dwnld_all:
		# Remove MIRI MRS image uncals, if they exist
		if obs_type == 'ifu':
			cmd = f'rm {data_dir}{prog_id}/Obs{obs_numb}/uncal/*mirimage_uncal.fits'
			os.system(cmd)
		
		# Moving rate and rateints files
		cmd = f'mv {prefix}jw{prog_id}{obs_numb}???*rate*.fits {data_dir}{prog_id}/Obs{obs_numb}/rate/'
		os.system(cmd)

		# Moving cal or calints files
		cmd = f'mv {prefix}jw{prog_id}{obs_numb}???*_cal*.fits {data_dir}{prog_id}/Obs{obs_numb}/cal/'
		os.system(cmd)

		# Moving x1d or x1dints files
		cmd = f'mv {prefix}jw{prog_id}{obs_numb}???*_x1d*.fits {data_dir}{prog_id}/Obs{obs_numb}/x1d/'
		os.system(cmd)

		if obs_type == 'slit':
			# Moving s2d files
			cmd = f'mv {prefix}jw{prog_id}{obs_numb}???*s2d.fits {data_dir}{prog_id}/Obs{obs_numb}/s2d/'
			os.system(cmd)
		elif obs_type == 'ifu':
			# Moving s3d files
			cmd = f'mv {prefix}jw{prog_id}{obs_numb}???*s3d.fits {data_dir}{prog_id}/Obs{obs_numb}/s3d/'
			os.system(cmd)

	# Write download log file
	today = datetime.date.today().strftime("%b-%d-%Y")
	clock = datetime.datetime.now().strftime("%H:%M:%S")
	logstr = f'Data downloaded from MAST at {clock} (local time) on {today}.'
	with open(f'{data_dir}{prog_id}/Obs{obs_numb}/uncal/download_log.txt','w') as f:
		f.write(logstr)


def select_spec_files(input_files, params):
	'''
	Selects for only the spectroscopic data files with science data from a list of input files

	Parameters
	----------
	input_files : numpy.ndarray
	   	List of input files.
	params : obj
		Pipeline run parameters object.
 
	Returns
	-------
	select_files : numpy.ndarray
   		Selected files.
	'''

	select = np.ones(len(input_files)).astype('bool')
	for i,fi in enumerate(input_files):
		header = fits.getheader(fi, 'PRIMARY')

		# Exclude target acquisition and imaging data files
		if 'TA' in header['EXP_TYPE'] or 'IMAGE' in header['EXP_TYPE']:
			select[i] = False

		# Exclude NRS2 data files for PRISM and M-grating NIRSpec observations
		if params.instrument == 'nirspec' and header['GRATING'][-1] != 'H' and header['DETECTOR'] == 'NRS2':
			select[i] = False
			
	select_files = input_files[select]

	return select_files


def read_cubes(input_files, suffix=''):
	'''
	Reads in s3d.fits input files and creates a list of IFU data cube objects

	Parameters
	----------
	input_files : list
	   	List of s3d files.
 
	Returns
	-------
	cube_set : list
   		List of ifu_cube data objects, populated with data frames 
   		and relevant information.
   	'''
	
	print(f"Importing {len(input_files)} files:")
	print(input_files)
	cube_set = [ifu_cube(fi, suffix=suffix) for fi in input_files]
	return cube_set


def read_spectra(input_files, suffix='', nod_subtract=False):
	'''
	Reads in s2d.fits input files and creates a list of slit spectrum objects

	Parameters
	----------
	input_files : list
		List of s2d files.
 
	Returns
	-------
	spectrum_set : list
   		List of slit_spectrum data objects, populated with data frames 
   		and relevant information.
	'''
	
	print(f"Importing {len(input_files)} files:")
	print(input_files)
	spectrum_set = [slit_spectrum(fi, suffix=suffix, nod_subtract=nod_subtract) for fi in input_files]
	return spectrum_set


def spec_combine(group, resultsdir, spec_bkg_sub=True, special_defringe=False, spec_sig_clip=5.0, window_half=10):
	'''
	Takes a group of spectrum files and returns a single outlier-cleaned, combined spectrum

	Parameters
	----------
	group : list
    	List of spectrum pickle files.
	resultsdir : string
		Directory for saving spectra and plots (usually same location as extracted spectra).
	spec_bkg_sub : bool; optional
		Choice for whether to utilize background-subtracted flux (True) or not (False).
	special_defringe : bool; optional
		Whether or not special defringing was carried out on spectrum.
	spec_sig_clip : float; optional
		Threshold for sigma-clipping of each spectrum prior to combining.
	window_half : int; optional
		Half-width of the moving median filter window. Default is 10.
	'''

	# Initialize lists for data
	nspec = len(group)
	wave_grp = []
	flux_grp = []
	fluxerr_grp = []
	mask_grp = []

	# Populate lists
	for i,file in enumerate(group):
		with open(file,'rb') as f:
			spectrum = pickle.load(f)

		wave = spectrum.wave
		if spec_bkg_sub:
			flux = spectrum.flux
		else:
			flux = spectrum.fluxbkg
		fluxerr = spectrum.fluxerr
		mask_close = spectrum.mask_close

		# Create meta data dictionary
		if i == 0:
			meta = dict()
			meta['program'] = spectrum.program
			meta['obs'] = spectrum.obs
			meta['visit'] = spectrum.visit
			meta['grating'] = spectrum.grating
			meta['detector'] = spectrum.detector
			meta['target'] = spectrum.targ
			meta['bkg_sub'] = spec_bkg_sub
			meta['defringe'] = special_defringe
			meta['instrument'] = spectrum.inst
			meta['jwst_ver'] = spectrum.meta['CAL_VER']
			meta['crds_context'] = spectrum.meta['CRDS_CTX']
			meta['date'] = spectrum.meta['DATE']
			meta['extr_aper_rad'] = spectrum.extr_aper_rad
			print(f"Working on {nspec} spectra in {meta['grating']} on {meta['detector']}...")

		# Mask bad entries
		if hasattr(flux, 'mask'):
			mask = flux.mask
			flux = flux.data
			fluxerr = fluxerr.data
		else: 
			mask = np.zeros(len(wave))
		mask[np.isnan(flux)] = 1
		mask[flux == 0] = 1
		mask[flux == 1e-8] = 1 					# points that were not handled in S3_special_defringe
		mask[np.isnan(fluxerr)] = 1
		mederr = np.nanmedian(fluxerr[flux != 0]/flux[flux != 0])
		mask[fluxerr > 10*mederr*flux] = 1 		# anomalously large relative flux uncertainties
		flux = np.ma.array(flux, mask=mask)
		fluxerr = np.ma.array(fluxerr, mask=mask)

		# Iteratively mask outliers
		print(f'Cleaning outliers in spectrum {i+1} of {nspec}...')
		iter = 0
		while iter < 10:
			wave, flux, fluxerr, mask = spec_clean(wave, flux, fluxerr, mask, spec_sig_clip, window_half)
			iter += 1

		wave_grp.append(wave)
		flux_grp.append(flux)
		fluxerr_grp.append(fluxerr)
		mask_grp.append(mask)

	# Check that the spectra have the same length
	lengths = np.array([len(flux) for flux in flux_grp])
	if len(np.unique(lengths)) > 1:
		raise ValueError(f"Spectra are not of the same length! {lengths}")

	# Check that wavelength arrays are (nearly) identical
	if nspec > 1:
		diff = np.absolute(np.diff(np.asarray(wave_grp),axis=0))
		if np.max(diff) > 1e-4:	# Allowance of 0.1 nm
			raise ValueError(f"MISMATCH: Wavelength arrays differ as much as {max(diff)}!")

	# Do manual masking of remaining outliers
	wave_grp, flux_grp, fluxerr_grp, mask_grp = spec_mask_manual(wave_grp, flux_grp, fluxerr_grp, mask_grp)

	# Mean averaging (weighted means are unreliable for small numbers of dithers!)
	master_mask = np.asarray(mask_grp)
	flux_grp = np.ma.array(flux_grp, mask=master_mask)
	fluxerr_grp = np.ma.array(fluxerr_grp, mask=master_mask)
	
	wave_combine = wave_grp[0]
	flux_combine = np.ma.mean(flux_grp, axis=0)
	fluxerr_combine = np.sqrt(np.ma.sum(fluxerr_grp**2, axis=0)) / np.sum(~master_mask.astype('bool'), axis=0)
	mask_combine = flux_combine._mask
	print(f'Spectra combined.')

	# Save individual dither spectra
	append = ''
	if not meta['bkg_sub']:
		append += '-nobkgsub'
	if meta['defringe']:
		append += '-defringed'
	ditherdir = f'{resultsdir}dither/'
	os.makedirs(ditherdir, exist_ok=True)
	for i in range(nspec):
		outf = f"{ditherdir}{meta['target'].replace(' ','').replace('/','')}_{meta['grating']}_{meta['detector']}_spectrum_Dit{i+1}{append}.dat"
		np.savetxt(outf, np.array([wave_grp[i], flux_grp[i]*1e12, fluxerr_grp[i]*1e12]).T, delimiter='\t', fmt='%.5f')

	# Clean combined spectrum again, since there are still some notable outliers...
	wave_combine, flux_combine, fluxerr_combine, mask = spec_clean(wave_combine, flux_combine, fluxerr_combine, mask_combine, spec_sig_clip, window_half)
	
	# Plot individual and combined spectra
	print(f'Plotting spectra...')
	spec_plot(wave_combine, flux_combine, fluxerr_combine, resultsdir, meta, [wave_grp, flux_grp, fluxerr_grp])

	# Return combined spectrum
	return wave_combine, flux_combine, fluxerr_combine, meta


def spec_clean(x, y, yerr, mask, sig, window_half):

	# Mask outliers using moving median filter
	for i in range(window_half,len(x)):
		if y.mask[i]:
			continue
		else:
			seg = np.ma.concatenate([y[i-window_half:i],y[i+1:i+window_half+1]])
			med,std = np.ma.median(seg),np.ma.std(seg)
			if abs(y[i]-med) > sig*std:
				mask[i] = 1

	mask_y = np.ma.array(y, mask=mask)
	mask_yerr = np.ma.array(yerr, mask=mask)
	
	return x, mask_y, mask_yerr, mask


def spec_clip(x1, y1, yerr1, x2, y2, yerr2):

	nspec = len(x1)
	if nspec != len(x2):
		raise RuntimeError('Lengths of input spectrum lists are not the same!')

	# Run through each pair of spectra one-by-one
	for i in range(nspec):
		# Trim object spectrum wavelength range
		select = np.where(np.isnan(y2[i]) == False)
		wavmin, wavmax = np.min(x2[i][select]), np.max(x2[i][select])
		clip = np.where((x1[i] > wavmin) & (x1[i] < wavmax))
		x1[i] = x1[i][clip]
		y1[i] = y1[i][clip]
		yerr1[i] = yerr1[i][clip]

		# Interpolate the stellar spectrum to the wavelength grid of the object spectrum
		f = interp1d(x2[i][select], y2[i][select], kind='cubic')
		ferr = interp1d(x2[i][select], yerr2[i][select], kind='cubic')
		x2[i] = x1[i]
		y2[i] = f(x1[i])
		yerr2[i] = ferr(x1[i])
		
	return x1, y1, yerr1, x2, y2, yerr2


def spec_mask_manual(xgrp, ygrp, yerrgrp, maskgrp):

	# Scaling factor to make text file more friendly
	rescale = 1e12

	print("Click+release to select masked points.")
	print("Press 'q' to mask.")

	# Define interactive rectangular selection function and widget
	def line_select_callback(eclick, erelease):
	    'eclick and erelease are the press and release events'
	    x1, y1 = eclick.xdata, eclick.ydata
	    x2, y2 = erelease.xdata, erelease.ydata
	    print("(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))
	    
	def toggle_selector(event):
	    if event.key in ['q'] and toggle_selector.RS.active: 
	        toggle_selector.RS.set_active(False)

	# Plot individual spectra and manually mask outliers using rectangle widget	
	colors = ['b','g','y','r','c','m']
	symbols = ['.','^','s','o']
	
	done = False
	while not done:
		fig, current_ax = plt.subplots(figsize=(15,15))
		xmin = []
		xmax = []
		ymin = []
		ymax = []
		for i in range(len(xgrp)):
			fmt = colors[i % 6] + symbols[int(i // 6)]
			unmasked = np.where(ygrp[i].mask==False)
			plt.plot(xgrp[i][unmasked], ygrp[i][unmasked]*rescale, fmt, label=f'Dither {i+1}')
			xmin.append(min(xgrp[i][unmasked]))
			xmax.append(max(xgrp[i][unmasked]))
			ymin.append(min(ygrp[i][unmasked]*rescale))
			ymax.append(max(ygrp[i][unmasked]*rescale))
		toggle_selector.RS = RectangleSelector(current_ax, line_select_callback,
                                  useblit=True, button=[1],
                                  minspanx=5, minspany=5, spancoords='pixels',interactive=True)
		
		# Get plotting range
		global_min = min(ymin)
		global_max = max(ymax)
		global_range = global_max - global_min
		current_ax.set_xlim(min(xmin) - 0.1, max(xmax) + 0.1)
		current_ax.set_ylim(global_min - 0.1*global_range, global_max + 0.1*global_range)

		plt.connect('key_press_event', toggle_selector)
		plt.xlabel('Wavelength [um]', fontsize=14)
		plt.ylabel(f'Flux [uJy]',fontsize=14)
		plt.legend(loc='upper right')
		plt.show()
		
		# Mask all points within selected box
		[x1,x2,y1,y2] = toggle_selector.RS.extents
		for i in range(len(xgrp)):
			ww = np.where((xgrp[i] > x1) & (xgrp[i] < x2) & (ygrp[i]*rescale > y1) & (ygrp[i]*rescale < y2))
			ygrp[i].mask[ww] = True 
			yerrgrp[i].mask[ww] = True
			maskgrp[i][ww] = 1

		# Continue until done
		inp = input("Done? [y/n]")
		if inp == 'y':
			done = True

	return xgrp, ygrp, yerrgrp, maskgrp


def spec_plot(x, y, yerr, resultsdir, meta, dither_grp):

	# Default output file for plot
	append = ''
	if not meta['bkg_sub']:
		append += '-nobkgsub'
	if meta['defringe']:
		append += '-defringed'
	outf = f"{resultsdir}Visit{meta['visit']}_{meta['grating']}_{meta['detector']}_spectra{append}.pdf"
	colors = ['b','g','y','r','c','m']
	symbols = ['.','^','s','o']

	fig = plt.figure(figsize=(15,5))
	ax = fig.add_subplot(111)
	unmasked = np.where(y.mask==False)
	ax.errorbar(x[unmasked],y[unmasked],yerr=yerr[unmasked],fmt='k.',capsize=0,label='Combined',zorder=5)
	for i in range(len(dither_grp[0])):
		fmt = colors[i % 6] + symbols[int(i // 6)]
		unmasked = np.where(dither_grp[1][i].mask==False)
		ax.errorbar(dither_grp[0][i][unmasked],dither_grp[1][i][unmasked],yerr=dither_grp[2][i][unmasked],fmt=fmt,capsize=0,label=f'Dither {i+1}')
	if meta['grating'] == 'PRISM':
		ax.xaxis.set_major_locator(MultipleLocator(0.5))
		ax.xaxis.set_minor_locator(MultipleLocator(0.1))
	elif meta['instrument'] == 'MIRI':
		ax.xaxis.set_major_locator(MultipleLocator(1.0))
		ax.xaxis.set_minor_locator(MultipleLocator(0.2))
	else:
		ax.xaxis.set_major_locator(MultipleLocator(0.2))
		ax.xaxis.set_minor_locator(MultipleLocator(0.05))
	ax.tick_params(labelsize=12)
	ax.set_xlabel('Wavelength [um]',fontsize=14)
	ax.set_ylabel(f'Flux [MJy]',fontsize=14)
	plt.legend(loc='best')
	plt.savefig(outf)
	plt.close()


def spec_save(wave, flux, fluxerr, resultsdir, meta):

	# Scaling factor to make text file more friendly
	rescale = 1e12

	print(f'Spectrum saved!')
	append = ''
	if not meta['bkg_sub']:
		append += '-nobkgsub'
	if meta['defringe']:
		append += '-defringed'
	outf = f"{resultsdir}{meta['target'].replace(' ','').replace('/','')}_Visit{meta['visit']}_{meta['grating']}_{meta['detector']}_spectrum{append}.dat"

	# Fill masked flux/fluxerr values with NaNs
	flux_out = flux.filled(np.nan) * rescale
	fluxerr_out = fluxerr.filled(np.nan) * rescale

	# Populate header with metadata
	header = f'''=====
Object:   {meta['target']}
Program:  {meta['program']}
ObsNum:   {meta['obs']}
Instr:    {meta['instrument']}
Visit:    {meta['visit']}
Grating:  {meta['grating']} 
Detector: {meta['detector']}
Bkg_sub:  {meta['bkg_sub']}
jwst_ver: {meta['jwst_ver']}
CRDS_cxt: {meta['crds_context']}
Proc_date:{meta['date']}
Extr_aper:{meta['extr_aper_rad']}
=====
Wavelength (um)\tFlux (uJy)\t Flux Error (uJy)
====='''

	# Save formatted arrays
	np.savetxt(outf, np.array([wave, flux_out, fluxerr_out]).T, delimiter='\t', fmt='%.6f', header=header)
