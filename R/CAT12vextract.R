#' @title CAT12vextract
#'
#' @description Extracts vertex-wise surface-based measures for each subject from a \href{https://neuro-jena.github.io/cat12-help/}{CAT12} preprocessed directory, resampled to a 32k mesh, and stores them all as a single .RDS file.
#' @details The function searches inside the CAT12 preprocessed for 32k meshes (.gii) with the user-selected measure, extracts the data from these files, and organizes the left and right hemisphere vertex data for each subject as rows in a N x 64984 data matrix within a .rds object. Python and reticulate are required as the \href{https://nipy.org/nibabel/}{NiBabel} package is used to import .gii files outputted by CAT12. 
#'
#' @param sdirpath A string object containing the path to the CAT12 subjects preprocessed directory. Default is the current working directory ("./").
#' @param filename A string object containing the desired name of the output RDS file. Default is 'CAT12_measure.rds' in the R temporary directory (tempdir()).
#' @param measure A string object containing the name of the measure of interest. Options are 'thickness', 'depth', 'fractaldimension', 'gyrification', and 'toroGI20mm'. Default is 'thickness.'
#' @param subj_ID A logical object to determine whether to return a list object containing both subject ID and data matrix. Subject IDs are assumed to be the top directory names in the sdirpath.
#' @param silent A logical object to determine whether messages will be silenced. Set to 'FALSE' by default
#' @param VWR_check A boolean object specifying whether to check and validate system requirements. Default is TRUE.
#' 
#' @returns A .RDSfile with a list containing 1. the list of subject IDs (first element) and 2. a surface data matrix object (second element), or only a data matrix object. The matrix has N subjects x M vertices dimensions and can be readily used by VertexWiseR statistical analysis functions. Each row corresponds to a subject (in the order they are listed in the folder) and contains the left to right hemispheres' vertex-wise values.
#' @examples
#' CAT12vextract(sdirpath="./", 
#' filename='thickness.rds', 
#' measure='thickness', 
#' subj_ID = TRUE, 
#' VWR_check=FALSE)
#' 
#' @importFrom stringr str_extract
#' @export

CAT12vextract=function(sdirpath="./", filename, measure='thickness', subj_ID = TRUE, silent=FALSE, VWR_check=TRUE)
{
  
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
  
  
  if (!file.exists(sdirpath)) { stop('The path indicated in sdirpath could not be found.')}

  ## filename check
  if (missing("filename")) {
    warning(paste0('No filename argument was given. The matrix object "CAT12_', measure,'.rds" will be saved in R temporary directory (tempdir()).\n'))
    filename=paste0(tempdir(),'/CAT12_', measure,'.rds')
  }
  
  ## Identify mesh files
    filelist=list.files(path = sdirpath, pattern="mesh.*\\.gii$", 
                      recursive=TRUE)
    #only mesh files with the requested measure
    measlist<-stringr::str_extract(filelist, 
                                    "(?<=mesh\\.).*?(?=\\.resampled)")
    filelist=filelist[which(measlist==measure)]
    
    #Export subject list
    sublist<-stringr::str_extract(filelist, "(?<=k\\.).*?(?=_T)") 
    
  
  ##Function stops if files not found
  if (length(filelist)==0)
  {return(message(paste('No mesh .gii files for', measure,'could be found in the set sdirpath')))} 
  
  #import nibabel package
  nibabel=reticulate::import("nibabel")
  
  ##read data and save data for each subject as rows in a data matrix
  CAT12_dat=matrix(NA, nrow=NROW(filelist), ncol=64984)
  
  for (sub in 1:NROW(filelist))
  {
    giiobj=nibabel$load(filelist[sub]);
    CAT12_dat[sub,]=giiobj$darrays[[1]]$data;
  }
  CAT12_dat=CAT12_dat[order(filelist),]
  
  ##recode NAs to 0 because NA values will result in warnings for some of the package's functions
  CAT12_dat[is.na(CAT12_dat)]=0;
  
  if(silent==FALSE) {message(paste0("Saving output as ", filename))}
  ##output file depending on subj_ID==T

  if(subj_ID==TRUE) {CAT12_dat=list(sub_list=sublist, 
                                       surf_obj=CAT12_dat)} 
  

  saveRDS(CAT12_dat, file=filename)

  if(silent==FALSE) {message("done!")}
  
  return(CAT12_dat)
}
