## OTHER VERTEX-WISE FUNCTIONS
############################################################################################################################
############################################################################################################################
	
## permutation functions for random subject effects
  ## Paired/grouped data points are first shuffled within subjects, then these pairs/groups are shuffled between subjects
  perm_within_between=function(random)
  {
    ##for groups of 2 or more (subjects with 2 or more measurements)
    perm.idx=rep(NA, length(random))
    for(count in 2:max(table(random)))
    {
      if(length(which(table(random)==count))>0)
      {
        sub.id=as.numeric(which(table(random)==count))
          if(length(sub.id)>1)
          {
            ##between group shuffling
            recode.vec=sample(sub.id)
            vec.idx=1
            for(sub in sub.id)
            {
              perm.idx[which(random==sub)]=sample(which(random==recode.vec[vec.idx])) ##sample— within subject shuffling
              vec.idx=vec.idx+1
            }   
            remove(vec.idx,recode.vec)  
          } else 
          {
            ##if only one subject has a certain count, between subject shuffling will not be possible, only within-subject shuffling will be carried out
            perm.idx[which(random==sub.id)]=sample(which(random==sub.id)) ##sample— within subject shuffling
          }
      }
    }
    ##for subjects with a single measurement
    sub.idx=which(is.na(perm.idx))
    if(length(sub.idx)>1)
    {
      perm.idx[sub.idx]=sample(sub.idx)  
    } else 
    {
      perm.idx[sub.idx]=sub.idx
    }
    return(perm.idx)
  }

  ## Paired/grouped data points are shuffled within subjects, order of subjects in the dataset remains unchanged
  perm_within=function(random)
  {
    ##for groups of 2 or more (subjects with 2 or more measurements)
    perm.idx=rep(NA, length(random))
  
    for(count in 2:max(table(random)))
    {
      if(length(which(table(random)==count)>0))
      {
        sub.id=as.numeric(which(table(random)==count))
        for(sub in sub.id)
        {
          perm.idx[which(random==sub)]=sample(which(random==sub))
        }  
      }
    }
    return(perm.idx)
  }

  ## Paired/grouped data points are shuffled between subjects, order of data points within subjects remains unchanged.
  perm_between=function(random)
  {
    ##for groups of 2 or more (subjects with 2 or more measurements)
    perm.idx=rep(NA, length(random))
    for(count in 2:max(table(random)))
    {
      if(length(which(table(random)==count))>0)
      {
        sub.id=as.numeric(which(table(random)==count))
        if(length(sub.id)>1)
        {
          ##between group shuffling
          recode.vec=sample(sub.id)
          vec.idx=1
          for(sub in sub.id)
          {
            perm.idx[which(random==sub)]=which(random==recode.vec[vec.idx])
            vec.idx=vec.idx+1
          }   
          remove(vec.idx,recode.vec)  
        }
      }
    }
    ##for subjects with a single measurement
    sub.idx=which(is.na(perm.idx))
    if(length(sub.idx)>1)
    {
      perm.idx[sub.idx]=sample(sub.idx)  
    } else 
    {
      perm.idx[sub.idx]=sub.idx
    }
    return(perm.idx)
  }

############################################################################################################################
############################################################################################################################

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

  #Check required python dependencies. If files missing:
  #Will prompt the user to get them in interactive session 
  #Will stop if it's a non-interactive session 
  if (VWR_check == TRUE & length(sys.calls()) <= 1){
    message("Checking for VertexWiseR system requirements ... ")
    #brainstat must be installed for non hippocampal surf to run get_edgelist, check can be skipped for hippocampus (so 14524 vertices)
    if(max(dim(t(surf_data)))==14524) {
    check = VWRfirstrun('python/conda only')}
    else {check = VWRfirstrun(n_vert=max(dim(t(surf_data))))}
    if (!is.null(check)) {return(check)} 
  } else if(interactive()==FALSE) { return(message('Non-interactive sessions need requirement checks'))}
  
  #if surface_data is a path to an object, reads it
  if(inherits(surf_data,'character')==TRUE)
  {
    surf_data=readRDS(surf_data)
    #if also contained a subject list, only the surface data is kept
    if(inherits(surf_data,'list')==TRUE)
    {surf_data=surf_data[[2]]}
  }
  
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
  } else {stop("surf_data vector should only contain 20484 (fsaverage5), 81924 (fsaverage6) or 14524 (hippocampal vertices) columns")}

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
############################################################################################################################
############################################################################################################################

## Efficient way to extract t statistics from linear regression models to speed up the permutation process
## adapted from https://stackoverflow.com/questions/15820623/obtain-t-statistic-for-regression-coefficients-of-an-mlm-object-returned-by-l
extract.t=function(mod,row)
{
  p = mod$rank
  df.residual=NROW(mod$residuals)-NROW(mod$coefficients)
  rdf = df.residual
  Qr = mod$qr
  p1 = 1L:p
  r = mod$residuals
  rss = colSums(r^2)
  resvar = rss/rdf
  R = chol2inv(Qr[p1, p1, drop = FALSE])  
  se = (sqrt(diag(R) %*% t(resvar)))[row,]
  est = mod$coefficients[row,]
  tval = est/se 
  return(tval)
}
############################################################################################################################
############################################################################################################################
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph components

##find clusters using edgelist
getClusters=function(surf_data,edgelist)
{ 
  n_vert=length(surf_data)
  
  #listing out non-zero vertices
  vert=which(surf_data!=0)
  
  #matching non-zero vertices with adjacency matrices to obtain list of edges connecting between the non-zero vertices
  edgelist0=edgelist[!is.na(match(edgelist[,1],vert)),]
  if(length(edgelist0)>2)  {edgelist1=edgelist0[!is.na(match(edgelist0[,2],vert)),]} 
  else if (length(edgelist0)==2)  ##if only a single edge was identified, edgelist will no longer be a Nx2 matrix, hence need to reshape it into a matrix
  { 
    edgelist0=matrix(edgelist0,ncol=2,nrow=1)
    edgelist1=edgelist0[!is.na(match(edgelist0[,2],vert)),]
  } else {edgelist1=0}
  remove(surf_data,vert,edgelist0)
  
  if(length(edgelist1)>2) #if at least 2 edges are identified
  {
    #extracting cluster-related info from list of non-zero edges
    com=igraph::components(igraph::graph_from_data_frame(edgelist1, directed = FALSE))
    clust.size=com$csize
    
    #cluster mappings
    clust.map=rep(NA,n_vert)
    clust.map[as.numeric(names(com$membership))]=com$membership
  
  } else if(length(edgelist1)==2) #bypass cluster extraction procedure if only 1 edge is identified
  {
    clust.size=2
    clust.map=rep(NA,n_vert)
    clust.map[edgelist1]=1
  } else #bypass cluster extraction procedure if no edges are identified
  {
    clust.map="noclusters"
    clust.size="noclusters"
  }
  return(list(clust.map,clust.size,edgelist1))
}
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
#' @param atlas A numeric integer object corresponding to the atlas of interest. 1=aparc/Desikan-Killiany-70, 2=Destrieux-148, 3=Glasser-360, 4=Schaefer-100, 5=Schaefer-200, 6=Schaefer-400. For hippocampal surface, the 'bigbrain' 10-parcels hippocampal atlas is used by default and ignores the option.
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
#' @param surf_data A numeric vector or object containing the surface data, either in fsaverage5 (1 x 20484 vertices) or fsaverage6 (1 x 81924 vertices) space. It can only be one row of vertices (not a cohort surface data matrix). 
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
    else if (n_vert==81924) {template="fsaverage6"} 
    else {stop("Only an surf_data vector with a length of 20484 (fsaverage5) or 81924 (fsaverage6) is accepted")}
  
  #load python libraries
    interpolate=reticulate::import("brainstat.mesh.interpolate", delay_load = TRUE)
    nibabel=reticulate::import("nibabel", delay_load = TRUE)

  #convert and export .nii file
    stat_nii = interpolate$`_surf2vol`(template, surf_data)
    nibabel$save(stat_nii,filename)
    message(filename)
  }

############################################################################################################################
############################################################################################################################
#' @title Decode surface data
#'
#' @description Correlates the significant clusters of an earlier vertex-wise analysis with a database of task-based fMRI and voxel-based morphometric statistical maps and associate them with relevant key words. Decoding currently works with surfaces in fsaverage5 space only."
#'
#' @details The \href{https://nimare.readthedocs.io/en/stable/index.html}{'NiMARE'} python module is used for the imaging decoding and is imported via the reticulate package. The function also downloads the \href{https://github.com/neurosynth/neurosynth-data}{'Neurosynth' database} in the package's inst/extdata directory (~8 Mb) for the analysis.
#'
#' @param surf_data A numeric vector or object containing the surface data,  in fsaverage5 (1 x 20484 vertices). It can only be one row of vertices (not a cohort surface data matrix). 
#' @param contrast A string object indicating whether to decode the positive or negative mask ('positive' or 'negative')
#' @param VWR_check A boolean object specifying whether to check and validate system requirements. Default is TRUE.
#'
#' @returns A data.frame object listing the keywords and their Pearson's R values
#' @examples
#' CTv = rbinom(20484, 1, 0.001) 
#' decoding = decode_surf_data(CTv, 'positive', VWR_check=FALSE);
#' head(decoding)
#' @importFrom reticulate import r_to_py
#' @export

##CT image decoding
decode_surf_data=function(surf_data,contrast="positive", VWR_check=TRUE) 
{
  
  #Check required python dependencies. If files missing:
  #Will prompt the user to get them in interactive session 
  #Will stop if it's a non-interactive session
  if (VWR_check == TRUE){
    message("Checking for VertexWiseR system requirements ... ")
    check = VWRfirstrun(requirement="neurosynth")
    if (!is.null(check)) {return(check)} 
  } else if(interactive()==FALSE) { return(message('Non-interactive sessions need requirement checks'))}
  
  
  # Check if all values are positive
 if  (all(surf_data >= 0)==TRUE & contrast=="negative")
 {stop('No negative cluster was identified in the surf_data.')}
  # Check if all values are negative
 if  (all(surf_data <= 0)==TRUE & contrast=="positive")
 {stop('No positive cluster was identified in the surf_data.')}
  
  
  #if neurosynth database is installed
  if(file.exists(system.file('extdata','neurosynth_dataset.pkl.gz', package='VertexWiseR'))==TRUE)
  {
  ##checks length
    if(is.vector(surf_data)) {n_vert=length(surf_data)} else {n_vert=ncol(surf_data)}
    if(n_vert==20484) {template="fsaverage5"}
    else if (n_vert==81924) {stop("decoding of fsaverage6-space image is current not implemented, please resample the image to fsaverage5 space")} 
    else {stop("Only an surf_data vector with a length of 20484 (fsaverage5) is accepted")}

    #check contrast
    if(contrast != "positive" & contrast != "negative")  {stop("contrast has to be either positive or negative")} 
  
  message("Converting and interpolating the surface data ... ")
  
  ##import python libraries
  interpolate=reticulate::import("brainstat.mesh.interpolate", delay_load = TRUE)
  discrete=reticulate::import("nimare.decode", delay_load = TRUE)
  nimare.dataset=reticulate::import("nimare.dataset", delay_load = TRUE)
  
  ##selecting contrasts
  if(contrast=="positive")
  {
    surf_data[is.na(surf_data)]=0
    surf_data[surf_data<0]=0
    surf_data[surf_data>0]=1
  } else if (contrast=="negative")
  {
    surf_data[is.na(surf_data)]=0
    surf_data[surf_data>0]=0
    surf_data[surf_data<0]=1
  }

  ##convert surf_data vector to nii image
  stat_labels=reticulate::r_to_py(surf_data)
  stat_nii = interpolate$`_surf2vol`(template, stat_labels)
  

  ##running the decoding procedure
  neurosynth_dset = nimare.dataset$Dataset$load(system.file("extdata/neurosynth_dataset.pkl.gz", package='VertexWiseR'))
  message("\u2713 \n Correlating input image with images in the neurosynth database. This may take a while ... ")
  decoder = discrete$ROIAssociationDecoder(stat_nii)
  decoder$fit(neurosynth_dset)

  ##compiling the results
  decoder_df = data.matrix(decoder$transform())
  row.names(decoder_df)=gsub(pattern = "terms_abstract_tfidf__",x=row.names(decoder_df), replacement = "")
  result=data.frame(row.names(decoder_df),round(as.numeric(decoder_df),3))
  colnames(result)=c("keyword","r")
  result=result[order(-result$r),]
  message("\u2713 \n")
  return(result)
  } 
}  


############################################################################################################################
############################################################################################################################

#' @title Edge list fetcher
#' @description Functions that allows to download edge list for vertices of a given surface template from the BrainStat database
#' @param template A string object referring to a brain surface template. VertexWiseR currently works with: 'fsaverage5', 'fsaverage6' and 'fslr32k'.
#' @returns A Nx2 matrix object listing each vertex of the surface template and the vertices adjacent to it (making an edge together).
#' @examples
#' edgelist=get_edgelist("fslr32k")
#' @noRd
#' 
get_edgelist=function(template)
{
  #Load brainstat tools
  brainstat.datasets=reticulate::import("brainstat.datasets")  
  brainstat.mesh.utils=reticulate::import("brainstat.mesh.utils")
  
  #Read new python enviroment
  Renvironpath=paste0(tools::R_user_dir(package='VertexWiseR'),'/.Renviron')
  if (file.exists(Renvironpath)) {readRenviron(Renvironpath)}
  #Brainstat data, will either be stored in default $HOME path or 
  #custom if it's been set via VWRfirstrun()
  if (Sys.getenv('BRAINSTAT_DATA')=="")
  {brainstat_data_path=fs::path_home()} else if 
  (!Sys.getenv('BRAINSTAT_DATA')=="") 
  {brainstat_data_path=Sys.getenv('BRAINSTAT_DATA')}
  #convert path to pathlib object for brainstat
  data_dir=paste0(brainstat_data_path,'/brainstat_data/surface_data/')
  
  #Loads template surfaces
  surf_template=brainstat.datasets$fetch_template_surface(template, join=TRUE, data_dir=data_dir)
  
  #Returns edge list
  return(brainstat.mesh.utils$mesh_edges(surf_template)+1)
}

#' @title MNI coordinates fetcher
#' @description Functions that allows to download MNI coordinates of a given surface template
#' @param template A string object referring to a brain surface template. VertexWiseR currently works with: 'fsaverage5', 'fsaverage6' and 'fslr32k'.
#' @returns A matrix with X columns corresponding to the template's vertices and 3 rows corresponding to each vertex's X,Y,Z coordinates in MNI space
#' @examples
#' MNImap = get_MNIcoords("fslr32k")
#' @noRd

get_MNIcoords=function(template)
{
  #Load brainstat tools
  brainstat.datasets=reticulate::import("brainstat.datasets")  
  brainspace.mesh.mesh_elements=reticulate::import("brainspace.mesh.mesh_elements")
  
  #Read new python enviroment
  Renvironpath=paste0(tools::R_user_dir(package='VertexWiseR'),'/.Renviron')
  if (file.exists(Renvironpath)) {readRenviron(Renvironpath)}
  #Brainstat data, will either be stored in default $HOME path or 
  #custom if it's been set via VWRfirstrun()
  if (Sys.getenv('BRAINSTAT_DATA')=="")
  {brainstat_data_path=fs::path_home()} else if 
  (!Sys.getenv('BRAINSTAT_DATA')=="") 
  {brainstat_data_path=Sys.getenv('BRAINSTAT_DATA')}
  #convert path to pathlib object for brainstat
  data_dir=paste0(brainstat_data_path,'/brainstat_data/surface_data/')
  
  #Loads template surfaces
  surf.template=brainstat.datasets$fetch_template_surface(template=template, join=TRUE, data_dir=data_dir)
  
  #Returns MNI coordinates
  return(t(brainspace.mesh.mesh_elements$get_points(surf.template)))
}



########################################################################################################################################################################################################################################################

#' @title fs_stats()
#'
#' @description Extracts descriptive statistics, for the whole-brain and subcortical region-of-interests (ROI), within a FreeSurfer subjects directory. It reads them from the aseg.stats file, as generated by the default FreeSurfer preprocessing pipeline.
#'
#' @param sdirpath A string object indicating the path to the 'FreeSurfer' subjects directory. Default is the current working directory ("./").
#' @param sublist A string object indicating the path to the subject list generated by SURFvextract as 'sublist.txt' (optional). This allows users to retrieve stats only from a selected list of subjects. The subject list is a list with 1 subject ID per line.
#' @param ROImeasure A string object indicating what summary measure to extract for the subocrtical ROIs. Choices include: 'NVoxels', 'Volume_mm3', 'StructName', 'normMean', 'normStdDev', 'normMin', 'normMax', and 'normRange'. Default is 'Volume_mm3'.
#'
#' @returns A data.frame object with N columns per aseg.stats measures and N row per subjects.
#' @examples
#' fs_stats(sdirpath="freesurfer_subjdir")
#' @importFrom stringr str_split
#'@importFrom utils read.table
#' @export

fs_stats=function(sdirpath="./", sublist, ROImeasure='Volume_mm3') 
{
  
  
#check if sdirpath contains aseg.stats files 
  dircount = dir(path=sdirpath, recursive=TRUE, pattern="aseg.stats", include.dirs = TRUE, full.names = TRUE)
  if (length(dircount)==0) { return(message('aseg.stats files could not be found in the set sdirpath')) }
  
#decide subject_IDs to retrieve if sublist argument given
  if(!missing(sublist)) 
  {
    sublist=read.table(sublist)
    
    #get parent directory of stats/aseg.stats
    subdirs = stringr::str_split(dir(path=sdirpath, recursive=TRUE,
                         pattern="aseg.stats"), pattern='/stats/') 
    subdirs= unlist(lapply(subdirs, function(x) x[1]))

    #select only subject IDs match that directory name and use that
    dircount=dircount[which(subdirs %in% sublist[,1])==TRUE]

  }  

  
  
#creates main data.frame
FS_STAT=data.frame(matrix(ncol=1,nrow=length(dircount)))
colnames(FS_STAT) <- 'subject_ID'

#for each aseg stats, extracts measure values
  for (i in dircount) 
  {
    #read aseg.stats (read.delim to export the first Measures)
    stattable=read.delim(i, sep="\t")
    
      for (l in 1:nrow(stattable)) #for each line
      {  
        #if subjectname in the line store subject ID
        if (grepl("subjectname", stattable[l,])) 
        {
         subject_ID=stringr::str_split(stattable[l,],
                      pattern="subjectname ")[[1]][2]
         FS_STAT$subject_ID[which(dircount==i)]=subject_ID
        }
        
        #for each Measure, save name and value
        if (grepl("Measure", stattable[l,])) 
        {
          measure=stringr::str_split(stattable[l,],pattern=",")
          meas_name=stringr::str_split(measure[[1]][1], pattern="Measure ")[[1]][2]
          #if no column for measure name create it
          if (meas_name %in% colnames(FS_STAT)==FALSE)
          {  FS_STAT[[meas_name]] <- NA }
          #append measure value
          meas_val=as.numeric(gsub("\\D", "", measure))
          FS_STAT[[which(dircount==i),meas_name]] <- meas_val
        }
      }
    
    #read aseg.stats again (to easily get ROI mean volume)
    stattable2=read.table(i)
    
      for (l2 in 1:nrow(stattable2))
      {
        #get ROI name
        roi_name=stattable2[l2,5]
        #if no column for measure name create it
        if (roi_name %in% colnames(FS_STAT)==FALSE)
        {  FS_STAT[[roi_name]] <- NA }
        
        #select column based on user-specified ROImeasure:
        if(ROImeasure=='NVoxels'){r=3}
        else if (ROImeasure=='Volume_mm3') {r=4}
        else if (ROImeasure=='normMean') {r=6}
        else if (ROImeasure=='normStdDev') {r=7}
        else if (ROImeasure=='normMin') {r=8}
        else if (ROImeasure=='normMax') {r=9}
        else if (ROImeasure=='normRange') {r=10}
        
        #append measure value
        roi_val=stattable2[l2,r]
        FS_STAT[[which(dircount==i),roi_name]] <- roi_val
      }
    
  }
  

return(FS_STAT)  
}

############################################################################################################################
############################################################################################################################

#' @title Model structure check
#' @description Ensures the surface, contrast, model and/or random objects in the analyses have the appropriate structures to enable vertex analyses to be run properly on the given variables. Also provides a default smoothing procedure if required by the user.
#' @returns An error message if an issue or discrepancy is found.
#' @noRd

model_check=function(contrast, model, random, surf_data, smooth_FWHM)
{
  #If the contrast/model is a tibble (e.g., taken from a read_csv output)
  #converts the columns to regular data.frame column types
  if ('tbl_df' %in% class(contrast) == TRUE) {
    if (inherits(contrast[[1]],"character")==TRUE) {contrast = contrast[[1]]
    } else {contrast = as.numeric(contrast[[1]])}
  } 
  if ('tbl_df' %in% class(model) == TRUE) {
    model=as.data.frame(model)
    if (NCOL(model)==1) {model = model[[1]]
    } else { for (c in 1:NCOL(model)) { 
      if(inherits(model[,c],"double")==TRUE) {model[,c] = as.numeric(model[,c])}
    }  }
  }
  
  #numerise contrast
  if(inherits(contrast,"integer")==TRUE) {contrast=as.numeric(contrast)}
  
  #check if nrow is consistent for model and surf_data
  if(NROW(surf_data)!=NROW(model))  {stop(paste("The number of rows for surf_data (",NROW(surf_data),") and model (",NROW(model),") are not the same",sep=""))}
  
  #recode random variable to numeric
  if(!is.null(random)) { random=match(random,unique(random)) }
  
  ##checks
  #check contrast for consistency with the model data.frame
  if(NCOL(model)>1)
  {
    for(colno in 1:(NCOL(model)+1))
    {
      if(colno==(NCOL(model)+1))  {warning("contrast is not contained within model")}
      
      if(inherits(contrast,"character")==TRUE) 
      {
        if(identical(contrast,model[,colno]))  {break} 
      } else 
      {
        if(identical(suppressWarnings(as.numeric(contrast)),suppressWarnings(as.numeric(model[,colno]))))  {break}
      }
    }
  }  else
  {
    if(inherits(contrast,"character")==TRUE) 
    {
      if(identical(contrast,model))  {colno=1} 
      else  {warning("contrast is not contained within model")}
    } else
    {
      if(identical(as.numeric(contrast),as.numeric(model)))  {colno=1}
      else  {warning("contrast is not contained within model")}
    }
  }
  
  #incomplete data check
  idxF=which(complete.cases(model)==FALSE)
  if(length(idxF)>0)
  {
    message(paste("The model contains",length(idxF),"subjects with incomplete data. Subjects with incomplete data will be excluded from the current analysis\n"))
    model=model[-idxF,]
    contrast=contrast[-idxF]
    surf_data=surf_data[-idxF,]
    if(!is.null(random)) {random=random[-idxF]}
  }
  
  #check categorical and recode variable
  if(NCOL(model)>1)
  {
    for (column in 1:NCOL(model))
    {
      if(inherits(model[,column],"character")==TRUE)
      {
        if(length(unique(model[,column]))==2)
        {
          message(paste("The binary variable '",colnames(model)[column],"' will be recoded with ",unique(model[,column])[1],"=0 and ",unique(model[,column])[2],"=1 for the analysis\n",sep=""))
          
          recode=rep(0,NROW(model))
          recode[model[,column]==unique(model[,column])[2]]=1
          model[,column]=recode
          contrast=model[,colno]
        } else if(length(unique(model[,column]))>2)    {stop(paste("The categorical variable '",colnames(model)[column],"' contains more than 2 levels, please code it into binarized dummy variables",sep=""))}
      }      
    }
  } else
  {
    if(inherits(model,"character")==TRUE) 
    {
      if(length(unique(model))==2)
      {
        message(paste("The binary variable '",colnames(model),"' will be recoded such that ",unique(model)[1],"=0 and ",unique(model)[2],"=1 for the analysis\n",sep=""))
        
        recode=rep(0,NROW(model))
        recode[model==unique(model)[2]]=1
        model=recode
        contrast=model
      } else if(length(unique(model))>2)    {stop(paste("The categorical variable '",colnames(model),"' contains more than 2 levels, please code it into binarized dummy variables",sep=""))}
    }      
  }
  
  
  #check if surf_data is a multiple-rows matrix and NOT a vector
  if (is.null(nrow(surf_data)) | nrow(surf_data)==1)
  {stop("The surface data must be a matrix containing multiple participants (rows).")}
  
  ##smoothing
  n_vert=NCOL(surf_data)
  if(is.null(smooth_FWHM))
  {
    message("smooth_FWHM argument was not given. surf_data will not be smoothed here.\n")
  } else if(smooth_FWHM==0) 
  {
    message("smooth_FWHM set to 0: surf_data will not be smoothed here.\n")
  } else if(smooth_FWHM>0) 
  {
    message(paste("surf_data will be smoothed using a ",smooth_FWHM,"mm FWHM kernel", sep=""))
    surf_data=smooth_surf(surf_data, FWHM=smooth_FWHM)
  }
  surf_data[is.na(surf_data)]=0
  
  
  
  ##########################################
  #Output the right elements to be analysed
    model_summary=list(model=model, contrast=contrast,
                       random=random, surf_data=surf_data,
                       colno=colno)

  

 return(model_summary) 
}