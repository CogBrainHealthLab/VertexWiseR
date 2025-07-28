#' @title Smooth surface
#'
#' @description Smooths surface data at defined full width at half maximum (FWHM) as per the corresponding template of surface data
#'
#' @param surf_data A N x V matrix object containing the surface data (N row for each subject, V for each vertex), in fsaverage5 (20484 vertices), fsaverage6 (81924 vertices), fslr32k (64984 vertices) or hippocampal (14524 vertices) space. See also Hipvextract(), SURFvextract() or FSLRvextract output formats. Alternatively, a string object containing the path to the surface object (.rds file) outputted by extraction functions may be given.
#' @param FWHM A numeric vector object containing the desired smoothing width in mm 
#' @param VWR_check A boolean object specifying whether to check and validate system requirements. Default is TRUE.
#'
#' @returns A matrix object with smoothed vertex-wise values
#' @examples
#' surf_data = readRDS(file = url(paste0("https://github.com",
#'"/CogBrainHealthLab/VertexWiseR/blob/main/inst/demo_data/",
#'"FINK_Tv_ses13.rds?raw=TRUE")))[1:3,]
#'
#' surf_data_smoothed=smooth_surf(surf_data, 10, VWR_check=FALSE);
#' @importFrom reticulate source_python
#' @export

## smooth surface data 
## FWHM input is measured in mm, which is subsequently converted into mesh units
smooth_surf=function(surf_data, FWHM, VWR_check=TRUE)
{
  #gets surface matrix if is surf_data is a list or path
  surf_data=get_surf_obj(surf_data)
  
  #Check required python dependencies. Will skip if smoothing function called as part of modelling smooth_FWHM argument
  #If files missing:
  #Will prompt the user to get them in interactive session 
  #Will stop if it's a non-interactive session 
  if (VWR_check == TRUE & deparse(sys.call(-1)[[1]]) != 'model_check'){
    message("Checking for VertexWiseR system requirements ... ")
    #brainstat must be installed for non hippocampal surf to run get_edgelist, check can be skipped for hippocampus (so 14524 vertices)
    if(max(dim(t(surf_data)))==14524) {
      check = VWRfirstrun('python/conda only')}
    else {check = VWRfirstrun(n_vert=max(dim(t(surf_data))))}
    if (!is.null(check)) {return(check)} 
  } else if(VWR_check == FALSE & interactive()==FALSE) { return(message('Non-interactive sessions need requirement checks'))}
  
  #Solves the "no visible binding for global variable" issue
  . <- mesh_smooth <- NULL 
  internalenv <- new.env()
  assign("mesh_smooth", mesh_smooth, envir = internalenv)
  
  #mesh_smooth() fails if surf_data is not a matrix object
  if (class(surf_data)[1] != 'matrix') {
    surf_data = as.matrix(surf_data) }
  
  ##source python function
  reticulate::source_python(paste0(system.file(package='VertexWiseR'),'/python/smooth.py'))
  
  n_vert=ncol(surf_data)
  ##select template, set its FWHM parameter and load its edgelist file
  
  if(n_vert==20484) 
  {
    edgelist <- get_edgelist('fsaverage5') 
    FWHM=FWHM/3.5 #converting mm to mesh units
  } else if(n_vert==81924) 
  {
    edgelist <- get_edgelist('fsaverage6') 
    FWHM=FWHM/1.4 #converting mm to mesh units
  } else if(n_vert==64984) 
  {
    edgelist <- get_edgelist('fslr32k') 
    FWHM=FWHM/2 #converting mm to mesh units
  } else if(n_vert==14524) 
  {
    edgelist_hip <- get('edgelist_hip') 
    edgelist <- edgelist_hip@data
    FWHM=FWHM/0.5 #converting m to mesh units
  } else {stop("surf_data vector should only contain 20484 (fsaverage5), 81924 (fsaverage6), 64984 (fslr32k) or 14524 (hippocampal vertices) columns")}
  
  #to mask out the 0-value vertices (e.g., medial wall), so as to prevent the border regions from being significantly diluted by the 0-value vertices	  
  idx0=which(colSums(data.matrix(surf_data))==0)
  if(length(idx0)>0)
  {	  
    edgelist=edgelist[-which(!is.na(match(edgelist[,1],idx0))),]
    edgelist=edgelist[-which(!is.na(match(edgelist[,2],idx0))),]
  }
  
  smoothed=mesh_smooth(surf_data,edgelist, FWHM)
  smoothed[is.na(smoothed)]=0
  return(smoothed)
}