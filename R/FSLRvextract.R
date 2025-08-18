#' @title FSLRvextract
#'
#' @description Extracts vertex-wise surface-based measures for each subject from HCP, fMRIprep or XCP-D preprocessed directories, and stores it as a single .RDS file.
#' @details The function searches for the HCP preprocessed directory by listing out files with certain suffixes, extracts the data from these files, and organizes the left and right hemisphere vertex data for each subject as rows in a N x 64984 data matrix within a .rds object. 
#'
#' @param sdirpath A string object containing the path to the HCP or fMRIprep preprocessed directory. Default is the current working directory ("./").
#' @param filename A string object containing the desired name of the output RDS file. Default is 'fslr32k.rds' in the R temporary directory (tempdir()).
#' @param dscalar A string object containing the filename suffix of the dscalar file. These dscalar files are named differently depending on the preprocessing pipeline used. Examples of filename suffixes are shown below
#' \itemize{
#'  \item `.thickness_MSMAll.32k_fs_LR.dscalar.nii` (HCP MSMAll pipeline)
#'  \item `.sulc_MSMAll.32k_fs_LR.dscalar.nii` (HCP MSMAll pipeline)
#'  \item `.thickness.32k_fs_LR.dscalar.nii` (HCP legacy pipeline)
#'  \item `.sulc.32k_fs_LR.dscalar.nii` (HCP legacy pipeline)
#'  \item `_space-fsLR_den-91k_thickness.dscalar.nii` (fMRIprep; using the `--cifti-output 91k` flag)
#'  \item `_space-fsLR_den-91k_curv.dscalar.nii` (fMRIprep; using the `--cifti-output 91k` flag)
#'  \item `_space-fsLR_den-91k_sulc.dscalar.nii` (fMRIprep; using the `--cifti-output 91k` flag)}
#' @param subj_ID A logical object to determine whether to return a list object containing both subject ID and data matrix.
#' @param silent A logical object to determine whether messages will be silenced. Set to 'FALSE' by default
#' @param VWR_check A boolean object specifying whether to check and validate system requirements. Default is TRUE.
#' 
#' @returns A .RDSfile with a list containing 1. the list of subject IDs (first element) and 2. a surface data matrix object (second element), or only a data matrix object. The matrix has N subjects x M vertices dimensions and can be readily used by VertexWiseR statistical analysis functions. Each row corresponds to a subject (in the order they are listed in the folder) and contains the left to right hemispheres' vertex-wise values.
#' @examples
#' dat_fslr32k=FSLRvextract(sdirpath="./", 
#' filename="dat_fslr32k.rds",
#' dscalar=".thickness_MSMAll.32k_fs_LR.dscalar.nii", 
#' subj_ID = TRUE, silent=FALSE)
#' @export

FSLRvextract=function(sdirpath="./", filename, dscalar, subj_ID = TRUE, silent=FALSE, VWR_check=TRUE)
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
  
  oldwd <- getwd() 
  
  if (!file.exists(sdirpath)) { stop('The path indicated in sdirpath could not be found.')}
  setwd(sdirpath)
  
  ## filename check
  if (missing("filename") & silent==FALSE) {
    warning(paste0('No filename argument was given. The matrix object "fslr32k.rds" will be saved in R temporary directory (tempdir()).\n'))
    filename=paste0(tempdir(),'/fslr32k.rds')
  }
  
  ## get filelists and subject lists
  filelist=list.files(pattern=dscalar, recursive=TRUE)
  sublist=gsub(dscalar, "",basename(filelist))
  
  ##Function stops if files not found
  if (length(filelist) ==0)
  {return(message(paste0('No *',dscalar, ' files could be found in the set sdirpath')))} 
  
  ##read data and save data for each subject as rows in a data matrix
  fslr32k_dat=matrix(NA, nrow=NROW(sublist), ncol=64984)
    
  reticulate::py_require("nibabel")
  reticulate::source_python(paste0(system.file(package='VertexWiseR'),'/python/read_cifti.py'))
  
  for (sub in 1:NROW(sublist))
  {
    if(silent==FALSE) {message(paste0("Processing (",sub,"/",NROW(sublist),") ",filelist[sub]," ..."))}
    fslr32k_dat[sub,]=read_cifti(filelist[sub])
  }
  fslr32k_dat=fslr32k_dat[order(sublist),]
  
  ##recode NAs to 0 because NA values will result in warnings for some of the package's functions
  fslr32k_dat[is.na(fslr32k_dat)]=0
  
  if(silent==FALSE) {message(paste0("Saving output as ",filename))}
  ##output file depending on subj_ID==T
  
  if(subj_ID==TRUE) {fslr32k_dat=list(sub_list=sublist, 
                                      surf_obj=fslr32k_dat)} 
  
  setwd(oldwd) #will restore user's working directory path on function break
  saveRDS(fslr32k_dat, file=filename)
  
  if(silent==FALSE) {message("done!")}
  
  return(fslr32k_dat)
}