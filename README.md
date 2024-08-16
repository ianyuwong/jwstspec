# jwstspec
A generalized pipeline for customized end-to-end processing and spectral extraction of JWST spectroscopic data 

(Under construction...)

There are two scripts that can be used to run the pipeline:
1. **run_jwstspec.py** - This is a generalized run function that allows for full functionality, including spectral extraction of both NIRSpec and MIRI spectroscopic data.
2. **run_jwstspec_simple.py** - This is an abridged run function exclusively for MIRI data that contains only the data download and processing steps.

Users should make copies of those files to edit and run.

### Installation
We recommend that users first create and activate a new conda environment:
```
conda create -n <env_name> python=3.11
conda activate <env_name>
```
If needed, install required packages (see below for list of packages):
```
pip install <pkg_name>
```
Finally, go to the directory where `jwstspec` will be installed and clone the git repository:
```
cd <repo_location>
git clone https://github.com/ianyuwong/jwstspec.git
```
No further installation is necessary. The sample run scripts are provided in the top-level directory.



### Pipeline workflow
The pipeline consists of four stages:
- **Stage 0:** This stage automatically downloads the data products from MAST using a command line call to `jwst_mast_query` and creates the pipeline-specific directory architecture that will be used in subsequent stages. For information about installing and configuring jwst_mast_query, consult the instructions on the [package GitHub page](https://github.com/spacetelescope/jwst_mast_query).
- **Stage 1:** This stage runs the uncalibrated data through Stage 1 of the [official JWST pipeline](https://github.com/spacetelescope/jwst) `jwst`. This stage also handles background subtraction: If the user supplies a background observation number, the pipeline processes both the science and background data. If `bkg_subtract='asn'` is selected, an association file is generated that is passed to the subsequent stage for background subtraction. If `bkg_subtract = 'pixel'`, manual image-minus-image subtraction is carried out. For NIRSpec data, the user can choose to run various defringing routines to correct for 1/f readnoise. The 'nsclean' option is recommended, which uses the corresponding step contained in the official JWST pipeline. 
- **Stage 2:** This stage runs Stage 2 of the official JWST pipeline on the countrate products generated in the previous stage to produce calibrated products. For MIRI MRS data, the *residual_fringe* and *ifu_rfcorr* steps are run by default.
- **Stage 3:** This stage runs Stage 3 of the official JWST pipeline to produce dither/nod-combined calibrated products. If desired, the pipeline then proceeds to spectral extraction, which can be carried out on either the Stage 2 or the Stage 3 data products. Various methods for spectral extraction are available, including standard aperture extraction and a special empirical PSF fitting routine. Simultaneous background estimation and subtraction can be performed. For MIRI MRS data, the user can choose to run the *ifu_rfcorr* step on the individual sub-band spectra. The final step cleans and averages individual spectra from one or more dithers to produce the final combined spectra. For MIRI MRS data, the spectral leak correction is run by default to produce corrected spectra.

For users who do not wish to use jwst_mast_query or would like to run the pipeline on a server, Stage 0 can be run with `download = False`, which will just set up the directory architecture without downloading data. The necessary data products can then be moved into the respective directories for subsequent stages of the pipeline.


### Required packges
- jwst
- jwst_mast_query
- numpy
- scipy
- astropy
- matplotlib
- photutils   (*only if running spectral extraction*)
- numbamisc   (*only if running `defringing = 'moving_median'` method for NIRSpec defringing*)
