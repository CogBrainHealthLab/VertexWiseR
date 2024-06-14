#' @title SURFvextract
#'
#' @description Extracts whole-brain vertex-wise surface-based measures for each subject in a 'FreeSurfer' output subjects directory, resamples the data to a common surface template, and stores it as a .rds file. This function requires the 'FreeSurfer' environment to be preset.
#' @details The function runs system shell commands that will produce in the set subjects directory: 1) a sorted list of subjects "sublist.txt"; 2) a link file to the selected surface fsaverage template. 3) left and right hemisphere .mgh maps outputted by 'FreeSurfer' 's mris_preproc. 
#' This function was currently not tested on a MacOS system.
#'
#' @param sdirpath A string object containing the path to the 'FreeSurfer' subjects directory. Default is the current working directory ("./").
#' @param filename A string object containing the desired name of the output RDS file. Default is 'brain_measure.rds' in the R temporary directory (tempdir()).
#' @param template A string object containing the name of surface template (available: 'fsaverage5', 'fsaverage6'). Default is fsaverage5.
#' @param measure A string object containing the name of the measure of interest. Options are thickness, curv, sulc, area, and volume (for freesurfer 7.4.1 or later). Default is thickness.
#' @param subj_ID A logical object stating whether to include subject IDs (folder names in the subjects directory) as a first column to the output matrix. Default is TRUE.
#'
#' @returns A .RDSfile with a list containing 1. the list of subject IDs (first element) and 2. a surface data matrix object (second element), or a data matrix object. The matrix can be used readily by VertexWiseR statistical analysis functions. Each row corresponds to a subject (in the same order as 1) and contains the left to right hemispheres' vertex-wise values.
#' @examples
#' SURFvextract(sdirpath = "freesurfer_subjdir", 
#' filename=paste0(tempdir(), "/CTv.rds"), template="fsaverage5", 
#' measure="curv") 
#' @importFrom freesurferformats read.fs.mgh
#' @importFrom utils read.delim
#' @export

SURFvextract=function(sdirpath="./", filename, template='fsaverage5', measure = 'thickness', subj_ID = TRUE) 
{ 
  
  if (missing("filename")) {
    warning(paste0('No filename argument was given. The matrix object "brain_', measure,'.rds will be saved in R temporary directory (tempdir()).\n'))
    filename=paste0(tempdir(),'/brain_',measure,'.rds')
  }

  
#check if sdirpath contains freesurfer surf folders 
dircount = dir(path=sdirpath, recursive=TRUE, pattern="surf$", include.dirs = TRUE)
if (length(dircount)==0) { return(message('FreeSurfer surface folders could not be found in the set sdirpath')) }
  
#list subject folders and excludes template folders
alldirs=dir(path=sdirpath, pattern="surf$", recursive=TRUE, include.dirs=TRUE) 
alldirs=alldirs[-grep("fsaverage5|fsaverage6|fsaverage", alldirs)]
#checks subject with specific surf measure data (rh.measure file) 
sublist=list.files(alldirs, pattern=paste0("rh.",measure), recursive=TRUE, full.names=TRUE)
#flags subjects with no appropriate surf measure data inside surf/
missinglist=dirname(alldirs[! alldirs %in% dirname(sublist)])
if (length(missinglist) != 0) { warning(paste("Some of the subject directories do not contain the required surface files:", toString(missinglist)))}

#Stores subjects IDs with available data in sublist.txt
sublist = dirname(dirname(sublist)) #remove path to files
write.table(sublist, file = paste0(sdirpath,'/sublist.txt'), row.names = FALSE, col.names = FALSE, quote = FALSE)

#Calls Freesurfer to extract vertex-wise thickness data from the sample and resample it to the fsaverage5 common-space surface; and concatenate it into mgh files
Sys.setenv(SUBJECTS_DIR=sdirpath)
system(paste0("ln -s $FREESURFER_HOME/subjects/", template, " -t $SUBJECTS_DIR \n
       mris_preproc --f $SUBJECTS_DIR/sublist.txt --target ", template, " --hemi lh --meas ", measure, " --out $SUBJECTS_DIR/lh.mgh \n 
       mris_preproc --f $SUBJECTS_DIR/sublist.txt --target ", template, " --hemi rh --meas ", measure, " --out $SUBJECTS_DIR/rh.mgh"));

#Reads mgh files to stores and assign the thickness values to each subject in a matrix object usable by VertexWiseR. Appends a column with the subject IDs if required by the user.
if (subj_ID == TRUE) 
{
SURFdata= t(rbind(drop(freesurferformats::read.fs.mgh(paste0(sdirpath,"lh.mgh"))),drop(freesurferformats::read.fs.mgh(paste0(sdirpath,"rh.mgh")))));
SURFdata=list(sublist,SURFdata); 
} else {
SURFdata=t(rbind(drop(freesurferformats::read.fs.mgh(paste0(sdirpath,"lh.mgh"))),drop(freesurferformats::read.fs.mgh(paste0(sdirpath,"rh.mgh")))));
}

saveRDS(SURFdata, file=filename)
}