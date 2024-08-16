#' @title FSLRvextract
#'
#' @description Extracts vertex-wise surface-based measures for each subject from HCP preprocessed directory, and stores it as a single .RDS file.
#' @details The function searches for the HCP preprocessed directory by listing out files with certain suffixes, extract the data from these files, and organize the left and right hemisphere vertex data for each subject as rows in a N x 64984 data matrix within a .rds object. 
#'
#' @param sdirpath A string object containing the path to the HCP preprocessed directory. Default is the current working directory ("./").
#' @param wb_path The filepath to the workbench folder that you have previously downloaded and unzipped
#' @param filename A string object containing the desired name of the output RDS file. Default is 'fslr32k_measure.rds' in the R temporary directory (tempdir()).
#' @param measure A string object containing the name of the measure of interest. Options are 'thickness' or 'sulc'. Default is 'thickness'.
#' @param MSMAll A logical object to determine whether the filename of the dscalar file includes the "MSMAll", which will be the case if the data is preprocessed using the more recent HCP MSMAll pipeline. Set to 'TRUE' by default
#' @param subj_ID A logical object to determine whether to return a list object containing both subject ID and data matrix.
#' @param silent A logical object to determine whether messages will be silenced. Set to 'FALSE' by default
#' 
#' @returns A .RDSfile with a list containing 1. the list of subject IDs (first element) and 2. a surface data matrix object (second element), or only a data matrix object. The matrix has N subjects x M vertices dimensions and can be readily used by VertexWiseR statistical analysis functions. Each row corresponds to a subject (in the order they are listed in the folder) and contains the left to right hemispheres' vertex-wise values.
#' @examples
#' FSLRvextract(sdirpath = "./", filename = paste0(tempdir(),"/fslr32k.RDS"), measure = "thickness") 
#' 
#' @importFrom ciftiTools ciftiTools.setOption readcii
#' @export

FSLRvextract=function(sdirpath="./", wb_path,MSMAll=TRUE, filename, measure="thickness", subj_ID = TRUE, silent=FALSE)
{
  oldwd <- getwd()
  on.exit(setwd(oldwd)) #will restore user's working directory path on function break
  
  setwd(sdirpath)
  
  ## check if an acceptable measure is entered
  if (is.na(match(measure,c("thickness","sulc"))))
  {
  stop("Invalid measure entered. Please enter one of the following in the measure argument—'thickness_MSMAll','sulc_MSMAll','thickness','sulc'")
  }
  
  if(MSMAll==TRUE) {measure=paste0(measure,"_MSMAll")}
  
  ## filename check
  if (missing("filename")) {
    warning(paste0('No filename argument was given. The matrix object "fslr32k_', measure,'.rds will be saved in R temporary directory (tempdir()).\n'))
    filename=paste0(tempdir(),'/fslr32k_',measure,'.rds')
  }
  
  ## get filelists and subject lists
  filelist=list.files(pattern=paste0(measure,".32k_fs_LR.dscalar.nii"), recursive=TRUE)
  sublist=gsub(paste0(".",measure,".32k_fs_LR.dscalar.nii"), 
               "",
               basename(filelist))
  
  ##Function stops if files not found
  if (length(filelist) ==0)
  {return(message(paste0('No *',paste0(measure,".32k_fs_LR.dscalar.nii"), 'files could be found in the set sdirpath')))} 
  
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
    fslr32k_dat[,c(LH.idx,RH.idx)]=c(dat.temp$data$cortex_left,dat.temp$data$cortex_right)
    remove(dat.temp)
  }
  
  if(silent==FALSE) {message(paste0("Saving output as ",filename))}
  ##output file depending on subj_ID==T
  fslr32k_dat=fslr32k_dat[order(sublist),]
  
  if(subj_ID==TRUE) {saveRDS(list(sublist,fslr32k_dat), file=filename)} 
  else  {saveRDS(fslr32k_dat, file=filename)}
  
  if(silent==FALSE) {message("done!")}
  
  return(fslr32k_dat)
}
