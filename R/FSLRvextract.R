#' @title FSLRvextract
#'
#' @description Extracts vertex-wise surface-based measures for each subject from HCP, fMRIprep or XCP-D preprocessed directories, and stores it as a single .RDS file.
#' @details The function searches for the HCP preprocessed directory by listing out files with certain suffixes, extracts the data from these files, and organizes the left and right hemisphere vertex data for each subject as rows in a N x 64984 data matrix within a .rds object. 
#'
#' @param sdirpath A string object containing the path to the HCP or fMRIprep preprocessed directory. Default is the current working directory ("./").
#' @param wb_path The filepath to the workbench folder that you have previously downloaded and unzipped
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
#' 
#' @returns A .RDSfile with a list containing 1. the list of subject IDs (first element) and 2. a surface data matrix object (second element), or only a data matrix object. The matrix has N subjects x M vertices dimensions and can be readily used by VertexWiseR statistical analysis functions. Each row corresponds to a subject (in the order they are listed in the folder) and contains the left to right hemispheres' vertex-wise values.
#' @examples
#' dat_fslr32k=FSLRvextract(sdirpath="./", 
#' wb_path="/path/to/workbench",
#' filename="dat_fslr32k.rds",
#' dscalar=".thickness_MSMAll.32k_fs_LR.dscalar.nii", 
#' subj_ID = TRUE, silent=FALSE)
#' 
#' @importFrom ciftiTools ciftiTools.setOption readcii
#' @export

FSLRvextract=function(sdirpath="./", wb_path,filename, dscalar, subj_ID = TRUE, silent=FALSE)
{
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
  
  ##load and configure ciftitools
  ciftiTools::ciftiTools.setOption('wb_path', wb_path)
  
  ##read data and save data for each subject as rows in a data matrix
  fslr32k_dat=matrix(NA, nrow=NROW(sublist), ncol=64984)
  
  for (sub in 1:NROW(sublist))
  {
    if(silent==FALSE) {message(paste0("Processing (",sub,"/",NROW(sublist),") ",filelist[sub]," ..."))}
    
    dat.temp=ciftiTools::readcii(filelist[sub],brainstructures = c("left","right"))
    LH.idx=which(dat.temp$meta$cortex$medial_wall_mask$left==TRUE)
    RH.idx=which(dat.temp$meta$cortex$medial_wall_mask$right==TRUE)+32492
    fslr32k_dat[sub,c(LH.idx,RH.idx)]=c(dat.temp$data$cortex_left,dat.temp$data$cortex_right)
    remove(dat.temp)
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
