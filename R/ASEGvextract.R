#' @title ASEGvextract
#'
#' @description Extracts surface-based measures from ASeg subcortical segmentations in each subject. These must be specifically computed from a FreeSurfer output subjects directory using the python toolbox "SubCortexMesh". ASEGvextract() reads through the metrics .vtk files outputted by SubCortexMesh (surface_metrics/ directory), extracts the vertex-wise values, and stores them as a .rds file, merging left and right hemispheres when applicable, for each (selected) subcortical region. 
#' @details The function looks for standard subject folder names (starting "sub-") and surface file names as defined by SubCortexMesh (starting "left-","right-", or "brain-stem"). ASEGvextract() requires reticulate and python to read the .vtk files. 
#' 
#' Number of vertices in the (bilateral) matrix for each region-of-interest:
#'  - 2044: accumbens area
#'  - 3430: amygdala
#'  - 6940: caudate nuclei
#'  - 39214: cerebellum
#'  - 8132: hippocampus
#'  - 3200: pallidum
#'  - 8394: putamen
#'  - 7768: thalamus
#'  - 7144: ventral diencephalon
#'  - 9452: brain stem
#'  - 95718: all aseg ROIs
#' @param sdirpath A string object containing the path to the 'SubCortexMesh' surface metrics directory, containing surface-based metrics (.vtk files). Default is the current working directory ("./").
#' @param outputdir A string object containing the path of the directory where all ROI RDS files will be stored. Default is 'subcortices' in the R temporary directory (tempdir()).
#' @param measure A string object containing the name of the measure of interest. Options include thickness, and surfarea. Default is thickness.
#' @param roilabel A string object or vector of string objects containing the name(s) of the subcortical region(s) to extract (both hemispheres automatically included when applicable). E.g. "thalamus" or c("thalamus", "caudate"). Default is all regions in sdirpath.
#' @param subj_ID A logical object to determine whether to return a list object containing both subject ID and data matrix.
#' @param silent A logical object to determine whether messages will be silenced. Set to 'FALSE' by default
#' @param VWR_check A boolean object specifying whether to check and validate system requirements. Default is TRUE.
#'
#' @returns A directory containing - for each bilateral subcortical region separately - a .RDS files, each with a list containing 1. the list of subject IDs (first element) and 2. a surface data matrix object (second element), or only a data matrix object. Each matrix has N subjects x M vertices dimensions and can be used readily by VertexWiseR statistical analysis functions. Each row corresponds to a subject (in the order they are listed in the folder) and contains the left to right hemispheres' vertex-wise values (if applicable).
#' @examples
#' ASEGvextract(sdirpath = "subcortexmesh_output_metrics", 
#' outputdir=paste0(tempdir(), "\\subcortices"), measure="surfarea") 
#' @importFrom reticulate import
#' @importFrom stringr str_extract
#' @export

ASEGvextract=function(sdirpath="./", outputdir, measure = 'thickness', roilabel, subj_ID = TRUE, silent=FALSE, VWR_check=TRUE) 
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
  
  #Output directory
  if (missing("outputdir")) {
    warning(paste0('No outputdir argument was given. The matrix objects will be saved in a directory named "subcortices" inside the R temporary directory (tempdir()).\n'))
    outputdir=paste0(tempdir(),'\\subcortices')
  } else if (!dir.exists(outputdir)) {
    dir.create(outputdir)
  }
  
  #list all Left/Right files
  lh.filelist=list.files(path = sdirpath,
                         pattern=paste0("left.+",measure),
                         recursive=T)
  rh.filelist=list.files(path = sdirpath,
                         pattern=paste0("right.+",measure),
                         recursive=T)
  
  #remove mismatched ROIs if one of the two bilateral missing
  mismatch=lh.filelist[which(! lh.filelist %in% gsub('right','left',rh.filelist))]
  mismatch=c(mismatch, rh.filelist[which(! rh.filelist %in% gsub('left','right',lh.filelist))])
  if(silent==FALSE & length(mismatch)>0){
    message('Given the bilateral counterpart of the surface is missing, the following will be dropped:'); 
    message(paste(mismatch, collapse = "\n"))
  }
  if (length(which(lh.filelist %in% mismatch))!=0)
  {lh.filelist=lh.filelist[-which(lh.filelist %in% mismatch)]}
  if (length(which(rh.filelist %in% mismatch))!=0)
  {rh.filelist=rh.filelist[-which(rh.filelist %in% mismatch)]}
  
  #Brainstem, if computed, added separately as no hemispheric split
  bstem.filelist=list.files(path = sdirpath, pattern=paste0("brain-stem.+",measure),recursive=T)
  #All aseg, similarly already completely merged
  allaseg.filelist=list.files(path = sdirpath, pattern=paste0("allaseg.+",measure),recursive=T)
  
  #subcortical ROIs list
  subcortical_list=unique(stringr::str_extract(lh.filelist, paste0("(?<=left-).+?(?=_", measure, ")")))
  if(length(bstem.filelist)>0) {subcortical_list=append(subcortical_list, 'brain-stem')}
  if(length(allaseg.filelist)>0) {subcortical_list=append(subcortical_list, 'allaseg')}
  
  #reduced to the roilabel if specified
  if (!missing(roilabel))
  {
    roilabel=gsub(pattern = 'left-|right-', '', roilabel) #take out hemisphere if left in
    subcortical_list=subcortical_list[which(subcortical_list %in% roilabel)]
  }  
    
  ##Function stops if files not found
  if (length(subcortical_list) == 0) 
  {return(message('Subcortical meshes could not be found in the set sdirpath.'))}
  
  #subject list
  sublist=unique(stringr::str_extract(lh.filelist, "sub-[^/]+"))
  
  #prepare VTK surface reader
  vtk=reticulate::import('vtk')
  #vtk to array converter
  vtk_getscalar=function(vtk_path)
  {
    reader=vtk$vtkPolyDataReader()
    reader$SetFileName(vtk_path); reader$Update();
    scalars=reader$GetOutput()
    scalars=scalars$GetPointData()$GetScalars()
    n <- scalars$GetNumberOfTuples()
    values <- numeric(n)
    for (i in 0:(n-1)) {values[i+1] <- scalars$GetValue(i)}
    return(values)
  }
  
  #prepare matrices for each subcortices separately
  #precalculated from fsaverage's aseg converted to surfaces in subcortexmesh
  aseg_matrices=list(
    "accumbens-area"=matrix(NA, nrow=NROW(sublist), ncol=2044),
    "amygdala"=matrix(NA, nrow=NROW(sublist), ncol=3430),
    "caudate"=matrix(NA, nrow=NROW(sublist), ncol=6940),
    "cerebellum-cortex"=matrix(NA, nrow=NROW(sublist), ncol=39214),
    "hippocampus"=matrix(NA, nrow=NROW(sublist), ncol=8132),
    "pallidum"=matrix(NA, nrow=NROW(sublist), ncol=3200),
    "putamen"=matrix(NA, nrow=NROW(sublist), ncol=8394),
    "thalamus"=matrix(NA, nrow=NROW(sublist), ncol=7768),
    "ventraldc"=matrix(NA, nrow=NROW(sublist), ncol=7144),
    "brain-stem"=matrix(NA, nrow=NROW(sublist), ncol=9452),
    "allaseg"=matrix(NA, nrow=NROW(sublist), ncol=95718)
    )
  
  #extract vtk scalars for each subject and subcortical region
  for (sub in 1:NROW(sublist))
  {
    if (silent==FALSE){message(paste0("Extracting ", sublist[sub], "\'s data... [",sub,"/", NROW(sublist),"]"))}
    #narrow down lists to particular subject
    lh.filelist.sub=lh.filelist[grep(sublist[sub],lh.filelist)]
    rh.filelist.sub=rh.filelist[grep(sublist[sub],rh.filelist)]
    bstem.filelist.sub=bstem.filelist[grep(sublist[sub],bstem.filelist)]
    allaseg.filelist.sub=allaseg.filelist[grep(sublist[sub],allaseg.filelist)]
      
    for (vol in 1:NROW(subcortical_list))
    {
      vol_label=subcortical_list[vol]
      if (silent==FALSE){message(paste0('=> ', vol_label))}
      
      if (vol_label!='brain-stem' & vol_label!='allaseg')
      {
        if (length(grep(vol_label,lh.filelist.sub)) != 0 ){
        lh_vtk_path=paste0(sdirpath, lh.filelist.sub[grep(vol_label,lh.filelist.sub)])
        lh=vtk_getscalar(lh_vtk_path)
        rh_vtk_path=paste0(sdirpath, rh.filelist.sub[grep(vol_label,rh.filelist.sub)])
        rh=vtk_getscalar(rh_vtk_path)
        aseg_matrices[[vol_label]][sub,]=c(lh,rh)
        }
      } else if (vol_label=='brain-stem')
      {
        if (length(bstem.filelist.sub) != 0 ){
        bstem_path=paste0(sdirpath,bstem.filelist.sub[grep(vol_label,bstem.filelist.sub)])
        aseg_matrices[[vol_label]][sub,]=vtk_getscalar(bstem_path)}
      } else if (vol_label=='allaseg')
      {  
        if (length(allaseg.filelist.sub) != 0 ){
        allaseg_path=paste0(sdirpath,allaseg.filelist.sub[grep(vol_label,allaseg.filelist.sub)])
        aseg_matrices[[vol_label]][sub,]=vtk_getscalar(allaseg_path)}
      }
    }   
  }
  
  if (silent==FALSE){message(paste0("Saving data to ", outputdir,"..."))}
  
  #output each ROI in its separate RDS matrix
  for (surf_dat in 1:length(aseg_matrices))
  {
    if (!all(is.na(aseg_matrices[[surf_dat]])))
    {
      surf_obj=aseg_matrices[[surf_dat]]
      
      if (subj_ID == TRUE) 
      {
        row.names(surf_obj)=sublist
        aseg_matrix=list(surf_obj=surf_obj, sub_list=sublist)
        #remove subjects with no data
        allna=which(apply(aseg_matrix[[1]], 1, 
                          function(x) all(is.na(x))))
        if (length(allna) > 0){
          aseg_matrix[[1]]=aseg_matrix[[1]][-allna,]
          aseg_matrix[[2]]=aseg_matrix[[2]][-allna]
        }
        #update object in the grand list to be returned 
        aseg_matrices[[surf_dat]]=aseg_matrix
      } else
      {
        allna=which(apply(surf_obj, 1, function(x) all(is.na(x))))
        if (length(allna) > 0) {aseg_matrix=surf_obj[-allna,]}
      }
      
      #save individual object separately
      saveRDS(aseg_matrix, file=paste0(outputdir,"/",tolower(names(aseg_matrices[surf_dat])),"_", measure,".rds"))
    }
  }
  
  if (silent==FALSE){message("done!")}
  
  #drop empty matrices
  missing_data=c()
  for (i in 1:length(aseg_matrices))
  {if (all(is.na(aseg_matrices[[i]]))) {missing_data=c(missing_data,i)}}
  if(!is.null(missing_data)){aseg_matrices=aseg_matrices[-missing_data]}
  
  return(aseg_matrices)
}