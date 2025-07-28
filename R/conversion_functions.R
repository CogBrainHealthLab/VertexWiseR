#Secondary functions for conversions of surface data:
#- surf_to_atlas()
#  Surface to atlas
#- atlas_to_surf()
#  Atlas to surface
#- fs5_to_fs6()
#  fsaverage5 to fsaverage6
#- fs6_to_fs5()
#  fsaverage6 to fsaverage5
#- surf_to_vol()
#  Surface to volume

############################################################################################################################
############################################################################################################################


#' @title Surface to atlas
#'
#' @description Returns the mean or sum of vertex-wise surface data for each ROI of a selected atlas
#' @details The function currently works with the aparc/Desikan-Killiany-70, Destrieux-148, Glasser-360, Schaefer-100, Schaefer-200, Schaefer-400 atlases. ROI to vertex mapping data were obtained from the \href{https://github.com/MICA-MNI/ENIGMA/tree/master/enigmatoolbox/datasets/parcellations}{'ENIGMA toolbox'} ; data for Destrieux came from \href{https://github.com/nilearn/nilearn/blob/a366d22e426b07166e6f8ce1b7ac6eb732c88155/nilearn/datasets/atlas.py}{ 'Nilearn' 's nilearn.datasets.fetch_atlas_surf_destrieux}
#' 
#' For hippocampal data, the function currently works with the "bigbrain" 10-parcels atlas integrated in 'HippUnfold.' See also \doi{doi:10.1016/j.neuroimage.2019.116328}.
#'
#' @param surf_data A N x V matrix object containing the surface data (N row for each subject, V for each vertex), in fsaverage5 (20484 vertices), fsaverage6 (81924 vertices), fslr32k (64984 vertices) or hippocampal (14524 vertices) space. See also Hipvextract(), SURFvextract() or FSLRvextract output formats.
#' @param atlas A numeric integer object corresponding to the atlas of interest.  1=Desikan, 2=Destrieux-148, 3=Glasser-360, 4=Schaefer-100, 5=Schaefer-200, 6=Schaefer-400. Set to `1` by default. This argument is ignored for hippocampal surfaces.
#' @param mode A string indicating whether to extract the sum ('sum') or the average ('mean') of the ROI vertices values. Default is 'mean'.
#'
#' @returns A matrix object with ROI as column and corresponding average vertex-wise values as row
#' @seealso \code{\link{atlas_to_surf}}
#' @examples
#' CTv = runif(20484,min=0, max=100)
#' parcel_data = surf_to_atlas(CTv, 1)
#' @export

surf_to_atlas=function(surf_data,atlas,mode='mean') 
{  
  #check length of vector or ncol of matrix
  if(max(dim(t(surf_data)))!=20484 & max(dim(t(surf_data)))!=81924 
     & max(dim(t(surf_data)))!=14524 & max(dim(t(surf_data)))!=64984) 
  {stop("Length of surf_data is neither 20484, 81924, 14524: the object is not compatible with the function")}
  
  #atlas argument needed if not hippocampal data
  if(missing("atlas") & max(dim(t(surf_data)))!=14524) {stop("Please specify an atlas number among the following: 1=aparc/Desikan-Killiany-70, 2=Destrieux-148, 3=Glasser-360, 4=Schaefer-100, 5=Schaefer-200, 6=Schaefer-400")}
  
  ###
  #mapping fsaverage5 space vertice to atlas (Nx20484 vertices)
  if(max(dim(t(surf_data)))==20484) 
  {
    #load atlas mapping surf_data
    ROImap_fs5 <- get('ROImap_fs5')
    ROImap <- list(ROImap_fs5@data,ROImap_fs5@atlases)
    #init variables
    nregions=max(ROImap[[1]][,atlas])
    #set NAs to 0
    surf_data[is.na(surf_data)]=0 
    
    if (is.vector(surf_data)==TRUE) {surf_data=rbind(matrix(surf_data,ncol=20484,nrow=1),NA); isavector=TRUE} #if vector, converts to matrix, and adds empty NA row to make object 2 dims
    
    ROI=matrix(NA, nrow=NROW(surf_data), ncol=nregions)
    
    if(mode=='mean') {
      for (region in 1:nregions)  {ROI[,region]=rowMeans(surf_data[,which(ROImap[[1]][,atlas]==region)])}
      if (exists("isavector"))  {ROI=ROI[1,]}
    } 
    else if (mode=='sum') 
    {
      for (region in 1:nregions)  {ROI[,region]=rowSums(surf_data[,which(ROImap[[1]][,atlas]==region)])}
      if (exists("isavector")) {ROI=ROI[1,]} #removes empty row if it was vector
    }
    else 
    {
      stop('\nPlease indicate a mode: only "sum" or "mean" are available.')
    }
    return(ROI)
  }
  
  ###
  #mapping fsaverage6 space vertice to atlas (Nx81924 vertices)
  if(max(dim(t(surf_data)))==81924) 
  {
    #load atlas mapping surf_data
    ROImap_fs6 <- get('ROImap_fs6')
    ROImap <- list(ROImap_fs6@data,ROImap_fs6@atlases)
    #init variables
    nregions=max(ROImap[[1]][,atlas])
    #set NAs to 0
    surf_data[is.na(surf_data)]=0 
    
    if (is.vector(surf_data)==TRUE) {surf_data=rbind(matrix(surf_data,ncol=81924,nrow=1),NA); isavector=TRUE} #if vector, converts to matrix, and adds empty NA row to make object 2 dims
    
    ROI=matrix(NA, nrow=NROW(surf_data), ncol=nregions)
    
    if(mode=='mean') {
      for (region in 1:nregions)  {ROI[,region]=rowMeans(surf_data[,which(ROImap[[1]][,atlas]==region)])}
      if (exists("isavector"))  {ROI=ROI[1,]}
    } 
    else if (mode=='sum') 
    {
      for (region in 1:nregions)  {ROI[,region]=rowSums(surf_data[,which(ROImap[[1]][,atlas]==region)])}
      if (exists("isavector")) {ROI=ROI[1,]} #removes empty row if it was vector
    }
    else 
    {
      stop('\nPlease indicate a mode: only "sum" or "mean" are available.')
    }
    return(ROI)
  }
  
  ###
  #mapping fslr32k space vertice to atlas (Nx64984 vertices)
  if(max(dim(t(surf_data)))==64984) 
  {
    #load atlas mapping surf_data
    ROImap_fslr32k <- get('ROImap_fslr32k')
    ROImap <- list(ROImap_fslr32k@data,ROImap_fslr32k@atlases)
    #init variables
    nregions=max(ROImap[[1]][,atlas])
    #set NAs to 0
    surf_data[is.na(surf_data)]=0 
    
    if (is.vector(surf_data)==TRUE) {surf_data=rbind(matrix(surf_data,ncol=64984,nrow=1),NA); isavector=TRUE} #if vector, converts to matrix, and adds empty NA row to make object 2 dims
    
    ROI=matrix(NA, nrow=NROW(surf_data), ncol=nregions)
    
    if(mode=='mean') {
      for (region in 1:nregions)  {ROI[,region]=rowMeans(surf_data[,which(ROImap[[1]][,atlas]==region)])}
      if (exists("isavector"))  {ROI=ROI[1,]}
    } 
    else if (mode=='sum') 
    {
      for (region in 1:nregions)  {ROI[,region]=rowSums(surf_data[,which(ROImap[[1]][,atlas]==region)])}
      if (exists("isavector")) {ROI=ROI[1,]} #removes empty row if it was vector
    }
    else 
    {
      stop('\nPlease indicate a mode: only "sum" or "mean" are available.')
    }
    return(ROI)
  }
  
  ###
  #mapping hippocampal space vertice to atlas (Nx14524 vertices)
  if(max(dim(t(surf_data)))==14524) 
  {
    #load atlas mapping surf_data
    ROImap_hip <- get('ROImap_hip')
    ROImap <- list(ROImap_hip@data,ROImap_hip@atlases)
    #init variables
    nregions=max(ROImap[[1]][,1])
    #set NAs to 0
    surf_data[is.na(surf_data)]=0 
    
    if (is.vector(surf_data)==TRUE) {surf_data=rbind(matrix(surf_data,ncol=14524,nrow=1),NA); isavector=TRUE} #if vector, converts to matrix, and adds empty NA row to make object 2 dims
    
    ROI=matrix(NA, nrow=NROW(surf_data), ncol=nregions)
    
    if(mode=='mean') {
      for (region in 1:nregions)  {ROI[,region]=rowMeans(surf_data[,which(ROImap[[1]][,1]==region)])}
      if (exists("isavector"))  {ROI=ROI[1,]}
    } 
    else if (mode=='sum') 
    {
      for (region in 1:nregions)  {ROI[,region]=rowSums(surf_data[,which(ROImap[[1]][,1]==region)])}
      if (exists("isavector")) {ROI=ROI[1,]} #removes empty row if it was vector
    }
    else 
    {
      stop('\nPlease indicate a mode: only "sum" or "mean" are available.')
    }
    return(ROI)
  }
  
}

############################################################################################################################
############################################################################################################################

#' @title Atlas to surface
#'
#' @description Maps average parcellation surface values (e.g. produced with the surf_to_atlas() function) to the fsaverage5, fsaverage6 or fslr32k space
#' @details The function currently supports the Desikan-Killiany-70, Schaefer-100, Schaefer-200, Schaefer-400, Glasser-360, or Destrieux-148 atlases for cortical surfaces, and the 'bigbrain' 10-parcels atlas for hippocampal surfaces. ROI to vertex mapping data for 1 to 4 were obtained from the \href{https://github.com/MICA-MNI/ENIGMA/tree/master/enigmatoolbox/datasets/parcellations}{'ENIGMA toolbox'} ; and data for 5 from \href{https://github.com/nilearn/nilearn/blob/a366d22e426b07166e6f8ce1b7ac6eb732c88155/nilearn/datasets/atlas.py}{'Nilearn' 's nilearn.datasets.fetch_atlas_surf_destrieux} . atlas_to_surf() will automatically detect the atlas based on the number of columns.
#'
#' @param parcel_data A matrix or vector object containing average surface measures for each region of interest, see the surf_to_atlas() output format. 
#' @param template A string object stating the surface space on which to map the data ('fsaverage5', 'fsaverage6', 'fslr32k', 'CIT168' (hippocampal)).
#'
#' @returns A matrix or vector object containing vertex-wise surface data mapped in fsaverage5, fsaverage6, fslr32k, or CIT168 space
#' @seealso \code{\link{surf_to_atlas}}
#' @examples
#' parcel_data = t(runif(100,min=0, max=100));
#' surf_data = atlas_to_surf(parcel_data, template='fsaverage5');
#' @export

atlas_to_surf=function(parcel_data, template) 
{
  #load atlas mapping surface data
  if (template=='fsaverage5') { 
    ROImap_fs5 <- get('ROImap_fs5'); n_vert=20484; 
    ROImap <- list(ROImap_fs5@data,ROImap_fs5@atlases)
  } else if (template=='fsaverage6') 
  { ROImap_fs6 <- get('ROImap_fs6'); n_vert=81924; 
  ROImap <- list(ROImap_fs6@data,ROImap_fs6@atlases)
  } else if (template=='fslr32k') 
  { ROImap_fslr32k <- get('ROImap_fslr32k'); n_vert=64984; 
  ROImap <- list(ROImap_fslr32k@data,ROImap_fslr32k@atlases)
  } else if (template=='CIT168') 
  { ROImap_hip <- get('ROImap_hip'); n_vert=14524; 
  ROImap <- list(ROImap_hip@data,ROImap_hip@atlases)
  } else { stop('The function currently only works with fsaverage5,  fsaverage6, fslr32k or CIT168 surfaces')}
  
  #checking template number based on dimensions (works for both matrices and vectors)
  if(template %in% c('fsaverage5','fsaverage6','fslr32k')) 
  {
    if (max(dim(t(parcel_data))) == 70) {atlas=1} 
    else if (max(dim(t(parcel_data))) == 148) {atlas=2} 
    else if (max(dim(t(parcel_data))) == 360) {atlas=3} 
    else if (max(dim(t(parcel_data))) == 100) {atlas=4} 
    else if (max(dim(t(parcel_data))) == 200) {atlas=5} 
    else if (max(dim(t(parcel_data))) == 400) {atlas=6}
    else { stop('The function could not identify what atlas your data was parcellated with, based on the number of columns or vectors (parcels). The cortical surfaces currently accept the aparc/Desikan-Killiany-70, Schaefer-100, Schaefer-200, Schaefer-400, Glasser-360, Destrieux-148 atlases.')}
  } else if (template=='CIT168')
  {
    if (max(dim(t(parcel_data))) == 10) {atlas=1}
    else { stop('Hippocampal CIT168 surfaces currently only accept the bigbrain 10-parcels atlas. The function could not identify the atlas from parcel_data, based on the number of columns or vectors (parcels).')}
  }
  
  
  if(length(dim(parcel_data))==2) #if parcel_data is a matrix
  {
    #init variables
    nregions=max(ROImap[[1]][,atlas])
    surf_dat=matrix(NA,nrow = NROW(parcel_data), ncol=n_vert)
    
    #mapping atlas label to fsaverage5 space
    for (sub in 1:NROW(parcel_data))
    {
      for (region in 1:nregions)  {surf_dat[sub,which(ROImap[[1]][,atlas]==region)]=parcel_data[sub,region]}      
    }
    
    
    ##########if parcel_data is a vector
  } else if(is.vector(parcel_data)==TRUE) 
  {
    #init variables
    nregions=max(ROImap[[1]][,atlas])
    surf_dat=rep(NA,n_vert)
    
    #mapping atlas label to the surface space
    for (region in 1:nregions)  {surf_dat[which(ROImap[[1]][,atlas]==region)]=parcel_data[region]}      
  }
  return(surf_dat)
}

############################################################################################################################
############################################################################################################################

#' @title fsaverage5 to fsaverage6
#'
#' @description Remaps vertex-wise surface data in fsaverage5 space to fsaverage6 space using the nearest neighbor approach 
#'
#' @param surf_data A N x V matrix object containing the surface data (N row for each subject, V for each vertex), in fsaverage5 (20484 vertices)  space. See also SURFvextract() output format. 
#'
#' @returns A matrix object containing vertex-wise surface data mapped in fsaverage6 space
#' @seealso \code{\link{fs6_to_fs5}}
#' @examples
#' CTv = runif(20484,min=0, max=100);
#' CTv_fs6 = fs5_to_fs6(CTv);
#' @export

#convert between fsaverage5 and fsaverage6 spacing
fs5_to_fs6=function(surf_data)
{
  #check length of vector
  if(length(surf_data)%%20484!=0) {stop("Length of surf_data is not a multiple of 20484")}
  
  #load atlas mapping surf_data
  fs6_to_fs5_map <- get('fs6_to_fs5_map')
  
  #mapping fsaverage5 to fsaverage6 space if surf_data is a vector length of 20484
  if(length(surf_data)==20484) {surf_data.fs6=surf_data[fs6_to_fs5_map]} 
  #mapping fsaverage5 to fsaverage6 space if surf_data is a Nx20484 matrix
  else {surf_data.fs6=surf_data[,fs6_to_fs5_map]}
  surf_data.fs6[is.na(surf_data.fs6)]=0
  return(surf_data.fs6)
}

#' @title fsaverage6 to fsaverage5
#'
#' @description Remaps vertex-wise surface data in fsaverage6 space to fsaverage5 space using the nearest neighbor approach
#'
#' @param surf_data A N x V matrix object containing the surface data (N row for each subject, V for each vertex), in fsaverage6 (81924 vertices) space. See also SURFvextract() output format.  
#'
#' @returns A matrix object containing vertex-wise surface data mapped in fsaverage5 space
#' @seealso \code{\link{fs5_to_fs6}}
#' @examples
#' surf_data = runif(81924,min=0, max=100);
#' fs5_data=fs6_to_fs5(surf_data)
#' @importFrom stats aggregate
#' @export


fs6_to_fs5=function(surf_data)
{
  #check length of vector
  if(max(dim(t(surf_data)))%%81924!=0) {stop("Length of surf_data is not a multiple of 81924")}
  
  #load atlas mapping surf_data
  fs6_to_fs5_map <- get('fs6_to_fs5_map')
  
  if(max(dim(t(surf_data)))==81924) #mapping fsaverage6 to fsaverage5 space if surf_data is a Nx81924 matrix
  {
    vert.idx=data.frame(fs6_to_fs5_map)
    
    if(inherits(surf_data, 'numeric'))
    {
      surf_data.fs5=aggregate(surf_data, list(vert.idx$fs6_to_fs5_map), FUN=mean)[,2] 
    }
    else if(inherits(surf_data, 'matrix'))
    {
      surf_data.fs5=matrix(NA,ncol=20484,nrow=nrow(surf_data))
      #if matrix, loops across rows
      for (i in 1:nrow(surf_data))
      {surf_data.fs5[i,]=aggregate(surf_data[i,], list(vert.idx$fs6_to_fs5_map), FUN=mean)[,2] 
      }
      
    }
  }
  return(surf_data.fs5)
}

############################################################################################################################
############################################################################################################################

#' @title Surface to volume
#'
#' @description Converts surface data to volumetric data (.nii file)
#'
#' @param surf_data A numeric vector or object containing the surface data, either in fsaverage5 (1 x 20484 vertices) or fsLR32k (1 x 64984 vertices) space. It can only be one row of vertices (not a cohort surface data matrix). 
#' @param filename A string object containing the desired name of the output .nii file (default is 'output.nii' in the R temporary directory (tempdir())).
#' @param VWR_check A boolean object specifying whether to check and validate system requirements. Default is TRUE.
#'
#' @returns A .nii volume file
#' @examples
#' CTv = runif(20484,min=0, max=100);
#' surf_to_vol(CTv, filename = paste0(tempdir(),'/volume.nii'), VWR_check=FALSE)
#' @importFrom reticulate import
#' @export

##converting surface to volumetric data and exporting it as a .nii file

surf_to_vol=function(surf_data, filename, VWR_check=TRUE)
{
  #Check required python dependencies. If files missing:
  #Will prompt the user to get them in interactive session 
  #Will stop if it's a non-interactive session 
  if (VWR_check == TRUE){
    message("Checking for VertexWiseR system requirements ... ")
    check = VWRfirstrun(requirement="conda/brainstat")
    if (!is.null(check)) {return(check)} 
  } else if(interactive()==FALSE) { return(message('Non-interactive sessions need requirement checks'))}
  
  if (missing("filename")) {
    message('No filename argument was given. The volume will be saved as "vol.nii" in R temporary directory (tempdir()).\n')
    filename=paste0(tempdir(),'/output.nii')
  }
  
  #check length of vector
  n_vert=length(surf_data)
  if(n_vert==20484) {template="fsaverage5"}
  else if (n_vert==64984) {template="fslr32k"} 
  else {stop("Only an surf_data vector with a length of 20484 (fsaverage5) or 64984 (fslr32k) is accepted")}
  
  #load python libraries
  interpolate=reticulate::import("brainstat.mesh.interpolate", delay_load = TRUE)
  nibabel=reticulate::import("nibabel", delay_load = TRUE)
  
  #convert and export .nii file
  stat_nii = interpolate$`_surf2vol`(template, surf_data)
  nibabel$save(stat_nii,filename)
  message(filename)
}
