# VertexWiseR v1.4.0 (TBC)

## NEW FEATURES
* New function getting_started(): allows users to run through the package's tutorials and to download all the demo data at once. Article vignettes now propose downloading the demo data in bulk rather than running the line every time.
* New function plot_overlay_surf(): allows users to plot two surfaces on top of each other, with varying opacity, color maps and value ranges (for cortical surfaces). A new tutorial/vignette called 'Overlaying plots and transparent thresholding' was also created to showcase it.
* New function DTSERIESvextract(): allows users to extract vertex-wise surface-based CIFTI dense time-series data from an individual dtseries .nii file from HCP, fMRIprep or XCP-D preprocessed directories.
* New options for plot_surf: alpha setting to reduce transparency of non-NaN vertices, and option to not plot the NaN vertices.

## FIXES
* When VWRfirstrun() was set to promptless, it checked the wrong directory path for brainstat's fslr32k data and considered it missing.
* Successful check message stopped showing except when promptless was set as FALSE, this was a mistake.
* Fix to ensure the right numpy version gets installed when installing miniconda, in case brainstat does not upon installation
* Fix to still run UV if RETICULATE_PYTHON was set to NA (a suggestion given in our article for users whose systems already have a python installed but want to dedicate VWR to a virtual environment)

# VertexWiseR v1.3.2

## FIXES
* Function argument: FSLRvextract()'s argument "dscaler" changed to "dscalar"

# VertexWiseR v1.3.1

## NEW FEATURES
* Slight modification of the cluster building process, which speeds up TFCE analyses.
* plot_surf() now includes the transparent_bg argument which uses a boolean to make the plot's background transparent (default is FALSE, white background).
* fslr32k can now be used in surf_to_vol() and decode_surf_data()
* BrainStat's Yeo parcellation data are no longer required by VWRfirstrun(). They were fetched by default in BrainStat's SLM script, when VertexWiseR used it with hippocampal surfaces. 
* Slightly clearer message in VWRfirstrun()'s check of the BrainStat data path: the actual user's home path is printed instead of "$HOME_DIR/".

## FIXES
* Documentation mistake: surf_to_vol() does NOT work with fsaverage6.

# VertexWiseR v1.3.0

## NEW FEATURES
* TFCE computation is now optimized for speed: previously, the parallel steps in the permutation loop were building the linear model again every time; now the model will be saved as a file in a temporary directory (tempdir(), which will automatically be cleaned up), and loaded within the loops instead of rebuilt.
 
## FIXES
* TFCE_vertex_analysis_mixed(): The random variable was not being  processed identically in BrainStat across the permuted and unpermuted models, because of a setting factorizing the variable in only one case. Now random variables are factorized in both. This error affected the coefficients and t-stat estimations, and by extension the results of Example 2, which no longer shows negative clusters and now displays a smaller positive cluster for TFCE. This makes the outcome more consistent with the RFT results (TFCE being more conservative instead of the opposite). We apologise for overlooking this inconsistency in the code. Please refer to Example 2 as presented in its vignette/[website page](https://cogbrainhealthlab.github.io/VertexWiseR/articles/VertexWiseR_Example_2.html) for the most up-to-date and accurate results.

# VertexWiseR v1.2.1

## NEW FEATURES
* Reticulate's [last update](https://posit.co/blog/reticulate-1-41/) allows users to install ephemeral Python environments with UV instead of requiring a stable Python/Miniconda installation. If users create their own with py_require() before running VertexWiseR, such environment will be selected automatically. If no Python environment is found, VertexWiseR now gives the choice to either install an ephemeral environment with UV, or to install Miniconda or Python via the classic ways.

## FIXES
* SURFvextract() now gives a proper error message if subjects' surface measure files could not be found. It also will get rid of the sublist.txt which was outputted automatically if subj_ID was set to FALSE. 
* Fix for messages that were silenced by mistake during Miniconda's installations process.
* Fix for the script using pip/pip3 to install properly vtk (9.3.1) when choosing the classic Python installation.

# VertexWiseR v1.2.0

## NEW FEATURES

* Update for RFT_vertex_analysis: now outputs also the unthresholded tstat map; Rdoc fixed accordingly (setting p=1 as previously advised caused errors, it doesn't simply output the unthresholded tstat map)
* SURFvextract() now has the optional argument fshomepath. This ensures FreeSurfer's environment is accessed by R and set up again if the function is used from RStudio and the FREESURFER_HOME variable has been defined before running Rstudio. This is needed because RStudio does not inherit the system variables set before opening it from a terminal. 
* If the surf_data argument (in modelling or smoothing functions) is a list object with the subject list along with the surface matrix (as outputted by extraction functions with subj_ID=TRUE), the code automatically detects the matrix in that list, named "surf_obj", instead of forcing users to specify the matrix. 
* Modelling functions now accept a string with the path to the .rds file outputted by extraction functions, instead of only the matrix itself as the surf_data argument
* New function CAT12vextract() which allows surface data resampled to 32k meshes in CAT12 to be extracted and converted to a surface .rds object. This works with any measure applicable for the 32k resampled meshes ('thickness', 'depth', 'fractaldimension', 'gyrification', and 'toroGI20mm'). 
* New vignette/article gives advice on how to solve various Python-related issues
* Numpy's version check is no longer present as pip replaces it upon BrainStat's installation with the compatible version

## FIXES
* Python package 'vtk' causes issues in latest 9.4.0 versions. The correct 9.3.1 version is now installed when installing Miniconda or reticulate's Python environment via VWRfirstrun(). 
* The cmap argument in plot_surf() is now converted to class color if not in that class by default
* Surface extracters fix: Working directory is restored before saving the RDS instead of on exit for HIPvextract() and FSLRvextract(). This ensures the filename is not interpreted as relative to the sdirpath (subjects directory).
* Fixed a parameter that was skipping VWR_check if smooth_surf() or TFCE_threshold() were nested to avoid VWR_check's repetition. This caused issues for non-interactive sessions/other nesting situations. Now, instead: For TFCE_threshold(), users calling it may manually set the argument VWR_check=FALSE if they want it to be skipped. For smooth_surf(), now VWR_check will be skipped if the parent function is identified as "model_check" as it will mean the smooth function is run from a function that already have a VWR_check (RFT_vertex_analysis, TFCE_vertex_analysis and TFCE_vertex_analysis_mixed functions). 
* Fix for non-miniconda Python installation: numpy, vtk and brainstat were installed via system('') but the optional use of pip or pip3 did not work properly in Windows due to a bash syntax error. The pip3 installation is now only triggered if an error occurs with the first system('pip ...') call. Furthermore, reticulate's environment is loaded immediately after the Python installation to make sure the pip function can be used subsequently to install the required packages. Other minor improvements in the messages printed during the VWRfirstrun() installation process.
* VWRfirstrun(): One message was still showing even if promptless=TRUE, this is now fixed
* FSLRvextract(): Default filename fixed (simply 'fslr32k.rds')
* To avoid instability, the reticulate Miniconda installation is now set at v.24.9.2 instead of the latest

# VertexWiseR v1.1.0

## NEW FEATURES
* New plotsurf_3d() function: allows surface data to be plotted in a 3D viewer via the plotly interface
* The [extraction tutorial](https://cogbrainhealthlab.github.io/VertexWiseR/articles/VertexWiseR_surface_extraction.html) now also gives FreeSurfer preprocessed data for demonstration	
* plot_surf() now accepts atlas ROI values from atlases supported by atlas_to_surf()
* atlas_to_surf() now works with hippocampal surface data
* smooth_surf() allows to enter a path to a surface .rds object instead of inputting the object itself

## FIXES

* For miniconda/python libraries installations via VWRfirstrun(), pip3 added as an alternative used if pip install fails
* Fixed Surface extraction tutorial (download.file had the wrong url/method, untar added as alternative as unzip had possible failures for the demo surface data)
* Fixes for MacOS (better custom path management, SURFvextract() compatibility issue fixed) 
