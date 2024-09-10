#' @title SURFvextract
#'
#' @description Extracts whole-brain vertex-wise surface-based measures for each subject in a 'FreeSurfer' output subjects directory, resamples the data to a common surface template, and stores it as a .rds file. This function requires the 'FreeSurfer' environment to be preset and a 'FreeSurfer' license key.
#' @details The function runs system shell commands that will produce in the set subjects directory: 1) a sorted list of subjects "sublist.txt"; 2) a link file to the selected surface fsaverage template. 3) left and right hemisphere .mgh maps outputted by 'FreeSurfer' 's mris_preproc. 
#'
#' @param sdirpath A string object containing the path to the 'FreeSurfer' subjects directory. Default is the current working directory ("./").
#' @param filename A string object containing the desired name of the output RDS file. Default is 'brain_measure.rds' in the R temporary directory (tempdir()).
#' @param template A string object containing the name of surface template (available: 'fsaverage5', 'fsaverage6'). Default is fsaverage5.
#' @param measure A string object containing the name of the measure of interest. Options are thickness, curv, sulc, area, and volume (for freesurfer 7.4.1 or later). Default is thickness.
#' @param subj_ID A logical object stating whether to include subject IDs (folder names in the subjects directory) as a first column to the output matrix. Default is TRUE.
#'
#' @returns A .RDSfile with a list containing 1. the list of subject IDs (first element) and 2. a surface data matrix object (second element), or only a data matrix object. The matrix has N subjects x M vertices dimensions and can be used readily by VertexWiseR statistical analysis functions. Each row corresponds to a subject (in the order they are listed in the folder) and contains the left to right hemispheres' vertex-wise values.
#' @examples
#' SURFvextract(sdirpath = "freesurfer_subjdir", 
#' filename=paste0(tempdir(), "/CTv.rds"), template="fsaverage5",
#' measure="curv") 
#' @importFrom freesurferformats read.fs.mgh
#' @importFrom utils read.delim write.table
#' @importFrom stringr str_remove
#' @export

SURFvextract=function(sdirpath="./", filename, template='fsaverage5', measure = 'thickness', subj_ID = TRUE) 
{ 
  
  #check if FREESURFER_HOME has been set  
  if(Sys.getenv('FREESURFER_HOME')=="") {return(message('No FREESURFER_HOME variable has been set. SURFvextract() will not be able to work without FreeSurfer.'))}
  #check if FREESURFER_HOME/license.txt exists  
  if(!file.exists(paste0(Sys.getenv('FREESURFER_HOME'),'/license.txt')) )
  {return(message('No FREESURFER_HOME/license.txt file was found. FreeSurfer requires a license to be fully operational. See https://surfer.nmr.mgh.harvard.edu/fswiki/License'))}
  #check if FreeSurfer has been setup (test freeview command)
  testFS <- suppressWarnings(try(system("which freesurfer", intern = TRUE, ignore.stderr = TRUE),silent=TRUE))
  if  (!length(testFS) > 0) {return(message('It seems that FreeSurfer has not been set up in your source environment as \'which freesurfer\' results in an error. SURFvextract() will not be able to work without FreeSurfer.'))}
  
  
  if (missing("filename")) {
    warning(paste0('No filename argument was given. The matrix object "brain_', measure,'.rds will be saved in R temporary directory (tempdir()).\n'))
    filename=paste0(tempdir(),'/brain_',measure,'.rds')
  }

#check if sdirpath contains last slash (will fail in unix if not)
if (substr(sdirpath, nchar(sdirpath), nchar(sdirpath)) != '/')
  {sdirpath=paste0(sdirpath,'/')}
  
#check if sdirpath contains freesurfer surf folders 
dircount = dir(path=sdirpath, recursive=TRUE, pattern="surf$", include.dirs = TRUE)
if (length(dircount)==0) { return(message('FreeSurfer surface folders could not be found in the set sdirpath')) }
  
#list subject folders and excludes template folders
alldirs=dir(path=sdirpath, pattern="surf$", recursive=TRUE, include.dirs=TRUE) 
#remove paths that lead to the fsaverage template, if they are included
onlysubjsurf=alldirs[-grep("fsaverage5|fsaverage6|fsaverage", alldirs)]
if(length(onlysubjsurf) > 0) {alldirs=onlysubjsurf}
#checks subject with specific surf measure data (rh.measure file) 
sublist=list.files(paste0(sdirpath, alldirs), pattern=paste0("rh.",measure), recursive=TRUE, full.names=TRUE)
#flags subjects with no appropriate surf measure data inside surf/
missinglist=dirname(alldirs[! paste0(sdirpath, alldirs) 
                            %in% dirname(sublist)])
if (length(missinglist) != 0) { warning(paste("Some of the subject directories do not contain the required surface files:", toString(missinglist)))}

#Extract subjects IDs from folders with available data in sublist.txt
sublist = dirname(dirname(sublist)) #remove subfolders
sublist = stringr::str_remove(sublist, 
                                 pattern=sdirpath) #remove full path
write.table(sublist, file = paste0(sdirpath,'/sublist.txt'), row.names = FALSE, col.names = FALSE, quote = FALSE)

#Calls Freesurfer to extract vertex-wise thickness data from the sample and resample it to the fsaverage5 common-space surface; and concatenate it into mgh files
Sys.setenv(SUBJECTS_DIR=sdirpath)
system(paste0("ln -s $FREESURFER_HOME/subjects/", template, " -t $SUBJECTS_DIR"), ignore.stderr = TRUE)
system(paste0("mris_preproc --f $SUBJECTS_DIR/sublist.txt --target ", template, " --hemi lh --meas ", measure, " --surfreg lh.sphere --out $SUBJECTS_DIR/lh.mgh \n 
       mris_preproc --f $SUBJECTS_DIR/sublist.txt --target ", template, " --hemi rh --meas ", measure, "--surfreg rh.sphere --out $SUBJECTS_DIR/rh.mgh"));

#Reads mgh files to stores and assign the thickness values to each subject in a matrix object usable by VertexWiseR. Appends a column with the subject IDs if required by the user.
if (subj_ID == TRUE) 
{
SURFdata= t(rbind(drop(freesurferformats::read.fs.mgh(paste0(sdirpath,"lh.mgh"))),drop(freesurferformats::read.fs.mgh(paste0(sdirpath,"rh.mgh")))));
SURFdata=list(sublist,SURFdata); 
} else {
SURFdata=t(rbind(drop(freesurferformats::read.fs.mgh(paste0(sdirpath,"lh.mgh"))),drop(freesurferformats::read.fs.mgh(paste0(sdirpath,"rh.mgh")))));
}

saveRDS(SURFdata, file=filename)
return(SURFdata)
}