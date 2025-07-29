#' @title DTSERIESvextract
#'
#' @description Extracts vertex-wise surface-based CIFTI dense time-series data from an individual dtseries .nii file from HCP, fMRIprep or XCP-D preprocessed directories, and stores it as a single .RDS file.
#' @details The function extracts the data from the dtseries.nii file provided, and organizes the left and right hemisphere vertex data for each subject as rows in a N x 64984 data matrix within a .rds object. 
#' @param dtseries A string object containing the full path to the dtseries files to extract from.
#' @param wb_path The filepath to the workbench folder that you have previously downloaded and unzipped
#' @param filename A string object containing the desired name of the output RDS file. Default is 'fslr32k.rds' in the R temporary directory (tempdir()).
#' @param silent A logical object to determine whether messages will be silenced. Set to 'FALSE' by default
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
#' wb_path="/path/to/workbench",
#' filename="demo_91k_bold.dtseries.rds", 
#' silent=FALSE)
#' 
#' 
#' ##visualizing e.g. the first 4 frames of the fMRI volume
#' #plot_surf(sub_dtseries[c(1,10,20,40),], 
#' #            file="4frames.png")
#'
#' @importFrom ciftiTools ciftiTools.setOption read_xifti
#' @export

DTSERIESvextract=function(dtseries, wb_path, filename, silent=FALSE)
{
  #dtseries exists?
  if(!file.exists(dtseries))
  {stop('The dtseries file given does not exist at the set path.')}
  #wb_path exists? (exits without error for Rdoc example)
  if(!dir.exists(wb_path))
  {return(message('The path to the workbench directory could be not found in the set wb_path.'))}
  
  ## filename check
  if (missing("filename") & silent==FALSE) {
    warning(paste0('No filename argument was given. The matrix object "dtseries" will be saved in R temporary directory (tempdir(): ', tempdir(), ').\n'))
    filename=paste0(tempdir(),'/dtseries.rds')
  }
  
  ##load and configure ciftitools
  ciftiTools::ciftiTools.setOption('wb_path', wb_path)

  ##reading and extracting data from the cifti file
  cifti.obj=ciftiTools::read_xifti(dtseries,
                       brainstructures = c("left","right"))
  n_frames=ncol(cifti.obj$data$cortex_left)
  cifti.timeseries.dat=matrix(0, nrow=n_frames,ncol=64984)
  
  fslr32k_dat=matrix(NA, nrow=n_frames, ncol=64984)
  LH.idx=which(cifti.obj$meta$cortex$medial_wall_mask$left==TRUE)
  RH.idx=which(cifti.obj$meta$cortex$medial_wall_mask$right==TRUE)+32492
  cifti.timeseries.dat[,c(LH.idx,RH.idx)]=t(rbind(cifti.obj$data$cortex_left,cifti.obj$data$cortex_right))
  
  if(silent==FALSE) {message(paste0("Saving output as ",filename))}
  
  saveRDS(cifti.timeseries.dat, file=filename)
  
  if(silent==FALSE) {message("done!")}
  
  return(cifti.timeseries.dat)
}

