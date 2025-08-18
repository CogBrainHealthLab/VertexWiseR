#' @title DTSERIESvextract
#'
#' @description Extracts vertex-wise surface-based CIFTI dense time-series data from an individual dtseries .nii file from HCP, fMRIprep or XCP-D preprocessed directories, and stores it as a single .RDS file.
#' @details The function extracts the data from the dtseries.nii file provided, and organizes the left and right hemisphere vertex data for each subject as rows in a N x 64984 data matrix within a .rds object. 
#' @param dtseries A string object containing the full path to the dtseries files to extract from.
#' @param filename A string object containing the desired name of the output RDS file. Default is 'fslr32k.rds' in the R temporary directory (tempdir()).
#' @param silent A logical object to determine whether messages will be silenced. Set to 'FALSE' by default
#' @param VWR_check A boolean object specifying whether to check and validate system requirements. Default is TRUE.
#' 
#' @returns A .RDSfile containing a surface data matrix object, with N time-point x M vertices dimensions and can be readily used by VertexWiseR statistical analysis functions. Each row corresponds to a time point in order and contains the left to right hemispheres' vertex-wise values.
#' @examples
#' #demo cifti dtseries from openneuro
#' #(ds005012, sub-18_ses-1_task-mid, run-01, 
#' #reduced to 50 time points)
#' download.file(paste0("https://github.com/CogBrainHealthLab",
#' "/VertexWiseR/blob/main/inst/demo_data/",
#' "demo_91k_bold.dtseries.nii?raw=TRUE"),
#' destfile=paste0(tempdir(),
#' "/demo_91k_bold.dtseries.nii"), 
#' mode = "wb")
#'              
#' sub_dtseries=DTSERIESvextract(
#' dtseries=paste0(tempdir(),
#'              "/demo_91k_bold.dtseries.nii"), 
#' silent=FALSE)
#' 
#' ##visualizing e.g. the first 4 frames of the fMRI volume
#' #plot_surf(sub_dtseries[c(1,10,20,40),], 
#' #            file="4frames.png")
#'
#' @export

DTSERIESvextract=function(dtseries, filename, silent=FALSE, VWR_check = TRUE)
{
  ##loading python library and functions
  #Check required python dependencies. If files missing:
  #Will prompt the user to get them in interactive session 
  #Will stop if it's a non-interactive session
  if (VWR_check == TRUE){
    if(silent==FALSE)
    {message("Checking for VertexWiseR system requirements ... ")
      check = VWRfirstrun('python/conda only')}
    else
    {check = VWRfirstrun('python/conda only', promptless = TRUE)}
    if (!is.null(check)) {return(check)}
  } else if(interactive()==FALSE) { return(message('Non-interactive sessions need requirement checks'))}
  
  #dtseries exists?
  if(!file.exists(dtseries))
  {stop('The dtseries file given does not exist at the set path.')}
  
  ## filename check
  if (missing("filename") & silent==FALSE) {
    warning(paste0('No filename argument was given. The matrix object "dtseries" will be saved in R temporary directory (tempdir(): ', tempdir(), ').\n'))
    filename=paste0(tempdir(),'/dtseries.rds')
  }

  ##reading and extracting data from the cifti file
  reticulate::py_require("nibabel")
  reticulate::source_python(paste0(system.file(package='VertexWiseR'),'/python/read_cifti.py'))
  cifti.timeseries.dat=t(read_cifti(dtseries))
  
  if(silent==FALSE) {message(paste0("Saving output as ",filename))}
  
  saveRDS(cifti.timeseries.dat, file=filename)
  
  if(silent==FALSE) {message("done!")}
  
  return(cifti.timeseries.dat)
}

