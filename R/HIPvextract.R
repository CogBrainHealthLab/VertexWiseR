#' @title HIPvextract
#'
#' @description Extracts hippocampal vertex-wise surface-based measures for each subject in the 'HippUnfold' subjects directory, and stores it as a single .RDS file.
#' @details The function searches for the hippocampal surface data by listing out files with certain suffixes, extract the data from these files, and organize the left and right hippocampal vertex data for each subject as rows in a N x 14524 data matrix within a .rds object. 
#'
#' @param sdirpath A string object containing the path to the 'HippUnfold' subjects directory. Default is the current working directory ("./").
#' @param filename A string object containing the desired name of the output RDS file. Default is 'hip_measure.rds' in the R temporary directory (tempdir()).
#' @param measure A string object containing the name of the measure of interest. Options are 'thickness','curvature','gyrification' and 'surfarea' (For more information see \href{https://hippunfold.readthedocs.io/en/latest/outputs/output_files.html#surface-metrics}{the 'HippUnfold' documentation}). Default is thickness.
#' @param subj_ID A logical object stating whether to return a list object containing both subject ID and data matrix.
#'
#' @returns A .RDSfile with a list containing 1. the list of subject IDs (first element) and 2. a surface data matrix object (second element), or only a data matrix object. The matrix has N subjects x M vertices dimensions and can be readily used by VertexWiseR statistical analysis functions. Each row corresponds to a subject (in the order they are listed in the folder) and contains the left to right hemispheres' hippocampal vertex-wise values.
#' @examples
#' HIPvextract(sdirpath = "./", filename = paste0(tempdir(),"/hip_data.RDS"), measure = "thickness") 
#' 
#' @importFrom gifti readgii
#' @export

HIPvextract=function(sdirpath="./", filename, measure="thickness", subj_ID = TRUE)
{
  oldwd <- getwd()
  on.exit(setwd(oldwd)) #will restore user's working directory path on function break
  
  setwd(sdirpath)
  
  if (missing("filename")) {
    warning(paste0('No filename argument was given. The matrix object "hip_', measure,'.rds will be saved in R temporary directory (tempdir()).\n'))
    filename=paste0(tempdir(),'/hip_',measure,'.rds')
  }
  
  ## get filelists and subject lists
  lh.filelist=list.files(pattern=paste("_hemi-L_space-T1w_den-0p5mm_label-hipp_",measure,".shape.gii",sep=""), recursive=TRUE)
  rh.filelist=gsub(paste("_hemi-L_space-T1w_den-0p5mm_label-hipp_",measure,".shape.gii",sep=""),
                   paste("_hemi-R_space-T1w_den-0p5mm_label-hipp_",measure,".shape.gii",sep=""),
                   lh.filelist)
  sublist=gsub(paste("_hemi-L_space-T1w_den-0p5mm_label-hipp_",measure,".shape.gii",sep=""), 
               "",
               basename(lh.filelist))
  
  ##Function stops if files not found
  if (length(lh.filelist) == 0 | length(rh.filelist) == 0)
  {return(message('HippUnfold data could not be found in the set sdirpath'))} 

  ##read data and save data for each subject as rows in a data matrix
  hip_dat=matrix(NA, nrow=NROW(sublist), ncol=14524)
  
  for (sub in 1:NROW(sublist))
  {
    lh=gifti::readgii(lh.filelist[sub])
    rh=gifti::readgii(rh.filelist[sub])
    hip_dat[sub,]=c(lh$data$normal,rh$data$normal)
  }

  ##output file depending on subj_ID==T
  hip_dat=hip_dat[order(sublist),]
    if(subj_ID==TRUE)
    {
    saveRDS(list(sublist,hip_dat), file=filename)
    } else
    {
    saveRDS(hip_dat, file=filename)
    }
}

