# VertexWiseR v1.2.0

## NEW FEATURES

* Update for RFT_vertex_analysis: now outputs also the unthresholded tstat map; Rdoc fixed accordingly (setting p=1 as previously advised caused errors, it doesn't simply output the unthresholded tstat map)
* SURFvextract() now has the optional argument fshomepath. This ensures FreeSurfer's environment is accessed by R and set up again if the function is used from RStudio. This is needed because RStudio does not inherit the system variables set before opening it from a terminal. 
* If the surf_data argument (in modelling or smoothing functions) is a list object with the subject list along with the surface matrix, the code automatically detects the matrix (now named "surf_obj" by the extracter functions) within that list, instead of forcing users to specify the matrix. 
* Modelling functions now accept a string with the path to the extracted surface .rds file instead of the matrix itself as the surf_data argument

## FIXES
* Python package 'vtk' causes issues in latest 9.4.0 versions. The correct 9.3.1 version is now installed when installing Miniconda or reticulate's Python environment via VWRfirstrun(). 
* The cmap argument in plot_surf() is now converted to class color if not in that class by default
* Surface extracters fix: Working directory is restored before saving the RDS instead of on exit for HIPvextract() and FSLRvetract(). This ensures the filename is not interpreted as relative to the sdirpath (subjects directory).

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
