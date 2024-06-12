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
#' @param surf_data A matrix object containing the surface data, see SURFvextract() or HIPvextract() output format
#' @param FWHM A numeric vector object containing the desired smoothing width in mm 
#'
#' @returns A matrix object with smoothed vertex-wise values
#' @examples
#' surf_data = readRDS(file = url(paste0("https://github.com",
#'"/CogBrainHealthLab/VertexWiseR/blob/main/inst/demo_data/",
#'"FINK_Tv_ses13.rds?raw=TRUE")))[1:3,]
#' surf_data_smoothed=smooth_surf(surf_data, 10);
#' @importFrom reticulate source_python
#' @export

## smooth surface data 
## FWHM input is measured in mm, which is subsequently converted into mesh units
smooth_surf=function(surf_data, FWHM)
{
  #Check required python dependencies. If files missing:
  #Will prompt the user to get them in interactive session 
  #Will stop if it's a non-interactive session 
  cat("Checking for VertexWiseR system requirements ... ")
  check = VWRfirstrun(requirement="miniconda only")
  if (!is.null(check)) {return(check)} else {cat("\u2713 \n")}
  
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
    edgelist<- get('edgelistfs5') 
    FWHM=FWHM/3.5 #converting mm to mesh units
  } else if(n_vert==81924) 
  {
    edgelist<- get('edgelistfs6') 
    FWHM=FWHM/1.4 #converting mm to mesh units
  } else if(n_vert==14524) 
  {
    edgelist<- get('edgelistHIP') 
    FWHM=FWHM/0.5 #converting m to mesh units
  } else {stop("surf_data vector should only contain 20484 (fsaverage5), 81924 (fsaverage6) or 14524 (hippocampal vertices) columns")}
  
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
getClusters=function(surf_data)
{ 
  n_vert=length(surf_data)
  
  #listing out non-zero vertices
  vert=which(surf_data!=0)
  
  #visible binding for edgelist object
  edgelist <- get("edgelist")
  
  ##narrow down to left or right hemisphere only to speed up matching process
  if(max(vert,na.rm = T)<=(n_vert/2))
  {
    nedgerow=min(which(edgelist>(n_vert/2)))-1
    edgelist=edgelist[1:nedgerow,]
  } else if(min(vert,na.rm = T)>(n_vert/2))
  {
    nedgerow=min(which(edgelist>(n_vert/2)))-1
    edgelist=edgelist[(nedgerow+1):NROW(edgelist),]
  }
  
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
    com=igraph::components(igraph::graph_from_data_frame(edgelist1, directed = F))
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
  return(list(clust.map,clust.size))
}
############################################################################################################################
############################################################################################################################

#' @title Surface to atlas
#'
#' @description Returns the mean or sum of vertex-wise surface data for each ROI of a selected atlas
#' @details The function currently works with the aparc/Desikan-Killiany-70, Destrieux-148, Glasser-360, Schaefer-100, Schaefer-200, Schaefer-400 atlases. ROI to vertex mapping data were obtained from the \href{https://github.com/MICA-MNI/ENIGMA/tree/master/enigmatoolbox/datasets/parcellations}{'ENIGMA toolbox'} ; data for Destrieux came from \href{https://github.com/nilearn/nilearn/blob/a366d22e426b07166e6f8ce1b7ac6eb732c88155/nilearn/datasets/atlas.py}{ 'Nilearn' 's nilearn.datasets.fetch_atlas_surf_destrieux}
#' 
#' For hippocampal data, the function currently works with the "bigbrain" atlas integrated in 'HippUnfold.' See also \href{https://doi.org/10.1016/j.neuroimage.2019.116328}{DeKraker et al., 2020}.
#'
#' @param surf_data A matrix object containing the surface data in fsaverage5 (20484 vertices), fsaverage6 (81924 vertices) or hippocampal (14524 vertices) space. See also Hipvextract() or SURFvextract() output format. 
#' @param atlas A numeric integer object corresponding to the atlas of interest. 1=aparc, 2=Destrieux-148, 3=Glasser-360, 4=Schaefer-100, 5=Schaefer-200, 6=Schaefer-400. For hippocampal surface, the 'bigbrain' hippocampal atlas is used by default and ignores the option.
#' @param mode A string indicating whether to extract the sum ('sum') or the average ('mean') of the ROI vertices values. Default is 'mean'.
#'
#' @returns A matrix object with ROI as column and corresponding average vertex-wise values as row
#' @seealso \code{\link{atlas_to_surf}}
#' @examples
#' CTv = runif(20484,min=0, max=100)
#' surf_to_atlas(CTv, 1)
#' @export

surf_to_atlas=function(surf_data,atlas,mode='mean') 
{  
  #check length of vector or ncol of matrix
  if(max(dim(t(surf_data)))!=20484 & max(dim(t(surf_data)))!=81924 
     & max(dim(t(surf_data)))!=14524) {stop("Length of surf_data is neither 20484, 81924, 14524: the object is not compatible with the function")}
  
  #atlas argument needed if not hippocampal data
  if(missing("atlas") & max(dim(t(surf_data)))!=14524) {stop("Please specify an atlas number among the following: 1=aparc, 2=Destrieux-148, 3=Glasser-360, 4=Schaefer-100, 5=Schaefer-200, 6=Schaefer-400")}
  
  #mapping fsaverage5 space vertice to atlas (Nx20484 vertices)
  if(max(dim(t(surf_data)))==20484) 
  {
    #load atlas mapping surf_data
    ROImap <- get('ROImap_fs5')
    #init variables
    nregions=max(ROImap[[1]][,atlas])
    #set NAs to 0
    surf_data[is.na(surf_data)]=0 
     
    if (is.vector(surf_data)==T) {surf_data=rbind(matrix(surf_data,ncol=20484,nrow=1),NA); isavector=T} #if vector, converts to matrix, and adds empty NA row to make object 2 dims
     
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
  
  
  #mapping fsaverage6 space vertice to atlas (Nx81924 vertices)
  if(max(dim(t(surf_data)))==81924) 
  {
    #load atlas mapping surf_data
    ROImap <- get('ROImap_fs6')
    #init variables
    nregions=max(ROImap[[1]][,atlas])
    #set NAs to 0
    surf_data[is.na(surf_data)]=0 
    
    if (is.vector(surf_data)==T) {surf_data=rbind(matrix(surf_data,ncol=81924,nrow=1),NA); isavector=T} #if vector, converts to matrix, and adds empty NA row to make object 2 dims
    
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
  
  #mapping hippocampal space vertice to atlas (Nx14524 vertices)
  if(max(dim(t(surf_data)))==14524) 
  {
    #load atlas mapping surf_data
    ROImap <- get('ROImap_HIP')
    #init variables
    nregions=max(ROImap[[1]][,1])
    #set NAs to 0
    surf_data[is.na(surf_data)]=0 
    
    if (is.vector(surf_data)==T) {surf_data=rbind(matrix(surf_data,ncol=14524,nrow=1),NA); isavector=T} #if vector, converts to matrix, and adds empty NA row to make object 2 dims
    
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
#' @description Maps average parcellation surface values (e.g. produced with the surf_to_atlas() function) to the fsaverage5 or fsaverage6 space
#' @details The function currently works with the Desikan-Killiany-70, Schaefer-100, Schaefer-200, Schaefer-400, Glasser-360, or Destrieux-148 atlases. ROI to vertex mapping data for 1 to 4 were obtained from the \href{https://github.com/MICA-MNI/ENIGMA/tree/master/enigmatoolbox/datasets/parcellations}{'ENIGMA toolbox'} ; and data for 5 from \href{https://github.com/nilearn/nilearn/blob/a366d22e426b07166e6f8ce1b7ac6eb732c88155/nilearn/datasets/atlas.py}{'Nilearn' 's nilearn.datasets.fetch_atlas_surf_destrieux} . atlas_to_surf() will automatically detect the atlas based on the number of columns.
#'
#' @param parcel_data A matrix or vector object containing average surface measures for each region of interest, see the surf_to_atlas() output format. 
#' @param template A string object stating the surface space on which to map the data ('fsaverage5' or 'fsaverage6').
#'
#' @returns A matrix or vector object containing vertex-wise surface data mapped in fsaverage5 or fsaverage6 space
#' @seealso \code{\link{surf_to_atlas}}
#' @examples
#' parcel_data = t(runif(100,min=0, max=100));
#' surf_data = atlas_to_surf(parcel_data, template='fsaverage5');
#' @export

atlas_to_surf=function(parcel_data, template) 
  {
    #load atlas mapping surface data
  if (template=='fsaverage5') { 
    ROImap <- get('ROImap_fs5'); n_vert=20484; 
  } else if (template=='fsaverage6') 
  { ROImap <- get('ROImap_fs6'); n_vert=81924; 
  } else { stop('The function currently only works with fsaverage5 and fsaverage6')}
  
    
 if(length(dim(parcel_data))==2) #if parcel_data is a matrix
  {
   if (ncol(parcel_data) == 70) {atlas=1} 
     else if (ncol(parcel_data) == 148) {atlas=2} 
     else if (ncol(parcel_data) == 360) {atlas=3} 
     else if (ncol(parcel_data) == 100) {atlas=4} 
     else if (ncol(parcel_data) == 200) {atlas=5} 
     else if (ncol(parcel_data) == 400) {atlas=6}
    else { stop('The function could not identify what atlas your data was parcellated with, based on the number of columns (parcels). The function currently works with the aparc/Desikan-Killiany-70, Schaefer-100, Schaefer-200, Schaefer-400, Glasser-360, or Destrieux-148 atlases.')}
    
    #init variables
    nregions=max(ROImap[[1]][,atlas])
    surf_dat=matrix(NA,nrow = NROW(parcel_data), ncol=n_vert)
    
    #mapping atlas label to fsaverage5 space
    for (sub in 1:NROW(parcel_data))
    {
      for (region in 1:nregions)  {surf_dat[sub,which(ROImap[[1]][,atlas]==region)]=parcel_data[sub,region]}      
    }
  } else if(is.vector(parcel_data)==T) #if parcel_data is a vector
  {
    if (length(parcel_data) == 70) {atlas=1} 
    else if (length(parcel_data) == 148) {atlas=2} 
    else if (length(parcel_data) == 360) {atlas=3} 
    else if (length(parcel_data) == 100) {atlas=4} 
    else if (length(parcel_data) == 200) {atlas=5} 
    else if (length(parcel_data) == 400) {atlas=6}
    else { stop('The function could not identify what atlas your data was parcellated with, based on the number of columns (parcels). The function currently works with the aparc/Desikan-Killiany-70, Schaefer-100, Schaefer-200, Schaefer-400, Glasser-360, or Destrieux-148 atlases.')}
    
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
#' @param surf_data A numeric vector or matrix object containing the surface data, see SURFvextract() output format. 
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
  return(surf_data.fs6)
}

#' @title fsaverage6 to fsaverage5
#'
#' @description Remaps vertex-wise surface data in fsaverage6 space to fsaverage5 space using the nearest neighbor approach
#'
#' @param surf_data A numeric vector or matrix object containing the surface data, see SURFvextract() output format. 
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

#' @title Surface plotter
#'
#' @description Plots surface data in a grid with one or multiple rows in a .png file
#'
#' @param surf_data  A numeric vector (length of V) or a matrix (N rows x V columns), where N is the number of subplots, and V is the number of vertices. It can be the output from SURFvextract() as well as masks or vertex-wise results outputted by analyses functions.
#' @param filename A string object containing the desired name of the output .png file.
#' @param title A string object for setting the title in the plot. Default is none. For titles that too long to be fully displayed within the plot, we recommend splitting them into multiple lines by inserting "\\n".
#' @param surface A string object containing the name of the type of cortical surface background rendered. Possible options include "white", "smoothwm","pial" and "inflated" (default). The surface parameter is ignored for hippocampal surface data.
#' @param cmap A string object specifying the name of an existing colormap or a vector of hexadecimal color codes to be used as a custom colormap. The names of existing colormaps are listed in the \href{https://matplotlib.org/stable/gallery/color/colormap_reference.html}{'Matplotlib' plotting library}. 
#' 
#' Default cmap is set to `"Reds"` for positive values, `"Blues_r"` for negative values and `"RdBu"` when both positive and negative values exist. 
#' @param limits A combined pair of numeric vector composed of the lower and upper color scale limits of the plot. If the limits are specified, the same limits will be applied to all subplots. When left unspecified, the same symmetrical limits c(-max(abs(surf_dat),max(abs(surf_dat))) will be used for all subplots. If set to NULL, each subplot will have its own limits corresponding to their min and max values
#' @param colorbar A logical object stating whether to include a color bar in the plot or not (default is TRUE).
#' @param size A combined pair of numeric vector indicating the image dimensions (width and height in pixels). Default is c(1920,400) for whole-brain surface and c(400,200) for hippocampal surface.
#' @param zoom A numeric value for adjusting the level of zoom on the figures. Default is 1.25 for whole-brain surface and 1.20 for hippocampal surface.
#'
#' @returns Outputs the plot as a .png image
#' @examples
#' results = runif(20484,min=0, max=1);
#' plot_surf(surf_data = results, filename='output.png',title = 
#' 'Cortical thickness', surface = 'inflated', cmap = 'Blues')
#' @importFrom reticulate tuple import np_array source_python
#' @importFrom grDevices col2rgb
#' @export

plot_surf=function(surf_data, filename, title="",surface="inflated",cmap,limits, colorbar=T, size, zoom)
{
  #Check required python dependencies. If files missing:
  #Will prompt the user to get them in interactive session 
  #Will stop if it's a non-interactive session 
  cat("Checking for VertexWiseR system requirements ...")
  check = VWRfirstrun(n_vert=max(dim(t(surf_data))))
  if (!is.null(check)) {return(check)} else {cat("\u2713 \n")}
  
  if (missing("filename")) {
    cat('No filename argument was given. The plot will be saved as "plot.png" in R temporary directory (tempdir()).\n')
    filename=paste0(tempdir(),'/plot.png')
  }
  
  #format title for single row
  if(is.null(nrow(surf_data)))
  {
    title=list('left'=list(title))
    rows=1
    surf_data=as.numeric(surf_data)
  } else {rows=nrow(surf_data)}

  #in a multi-row data scenario: insert a dummy title if title is missing  or repeat the title nrow times
  if(rows>1) 
    {
       if(missing("title")) {title=rep(NULL,rows)}
       else if (missing("title")) {title=rep(title,rows)}
    }
	  
  #check length of vector
  n_vert=length(surf_data)
  if(n_vert%%20484==0) {template="fsaverage5"}
  else if (n_vert%%81924==0) {template="fsaverage6"} 
  else if (n_vert%%14524!=0) {stop("surf_data vector should only contain 20484 (fsaverage5), 81924 (fsaverage6) or 14524 (hippocampal vertices) columns")}

  #if cmap is missing, select cmaps depending on whether the image contains positive only or negative only values
  if(missing("cmap"))
  {
    if(range(surf_data,na.rm = T)[1]>=0)  {cmap="Reds"}
    else if (range(surf_data,na.rm = T)[2]<=0)  {cmap="Blues_r"}
    else  {cmap="RdBu_r"}  
  }
	
  #custom cmap— if a vector of hex color codes is specified
  if(inherits(cmap,"colors")==T)
  {
    matplotlib=reticulate::import("matplotlib")
    
    custom_colors=t(col2rgb(cmap)/255) # convert hex color codes to RGB codes, then divide by 255 to convert to RGBA codes
    
    #save as python cmap object
    mymap = matplotlib$colors$LinearSegmentedColormap$from_list('my_colormap', custom_colors)
    matplotlib$colormaps$unregister(name = "custom_cmap")
    matplotlib$colormaps$register(cmap = mymap,name="custom_cmap")
    cmap="custom_cmap"  
  }
	
  #setting color scale limits
   maxlimit=max(abs(range(surf_data,na.rm = T)))
    if(missing("limits")) 
    {
      if(range(surf_data,na.rm = T)[1]>=0) {limits=reticulate::tuple(0,range(surf_data,na.rm = T)[2])} ##if image contains all positive values
      else if(range(surf_data,na.rm = T)[2]<=0) {limits=reticulate::tuple(range(surf_data,na.rm = T)[1],0)} ##if image contains all negative values
      else {limits=reticulate::tuple(-maxlimit,maxlimit)} ##symmetrical limits will be used if image contains both positive and negative values
    } else {
      ##user specified limits
      if(!is.null(limits))
      {
        limits=reticulate::tuple(limits[1],limits[2])  
      }   
    }
  
  if(n_vert%%14524!=0)
  {
  ##cortical surface fplots
    #import python libraries
    brainstat.datasets=reticulate::import("brainstat.datasets", delay_load = TRUE)  
    brainspace.plotting=reticulate::import("brainspace.plotting", delay_load = TRUE)  
    
    #loading fsaverage surface
    left=brainstat.datasets$fetch_template_surface(template, join=F, layer=surface)[1]
    right=brainstat.datasets$fetch_template_surface(template, join=F, layer=surface)[2]
    
    #default cortical size and zoom parametes
    if(missing("size")) { size=c(1920,rows*400)}
    if(missing("zoom")) { zoom=1.25 }
    
    surf_plot=brainspace.plotting$plot_hemispheres(left[[1]], right[[1]],  array_name=reticulate::np_array(surf_data),cmap=cmap, 
                                                size=reticulate::tuple(as.integer(size)),nan_color=reticulate::tuple(0.7, 0.7, 0.7, 1),
                                                return_plotter=T,background=reticulate::tuple(as.integer(c(1,1,1))),zoom=zoom,color_range=limits,
                                                label_text=title,interactive=F, color_bar=colorbar,  transparent_bg=FALSE)  ##disabling interactive mode because this causes RStudio to hang
  } else
  {
    #Solves the "no visible binding for global variable" issue
    . <- surfplot_canonical_foldunfold  <- NULL 
    internalenv <- new.env()
    assign("surfplot_canonical_foldunfold", surfplot_canonical_foldunfold, envir = internalenv)
    
  ##hippocampal plots
    #import python libraries
    reticulate::source_python(paste0(system.file(package='VertexWiseR'),'/python/hipp_plot.py'))
    
    #default hippocampal size and zoom parametes
    if(missing("size")) { size=c(400,200)}
    if(missing("zoom")) { zoom=1.2 }

    #reshaping surf_data into a 7262 x 2 x N array
    if(is.null(nrow(surf_data)))  {surf_data=cbind(surf_data[1:7262],surf_data[7263:14524])} #if N=1
    else  
      {
        surf_data.3d=array(NA,c(7262,2,nrow(surf_data))) #if N>1
        for (row in 1:nrow(surf_data))  {surf_data.3d[,,row]=cbind(surf_data[row,1:7262],surf_data[row,7263:14524])}
        surf_data=surf_data.3d
      }
    
    surf_plot=surfplot_canonical_foldunfold(surf_data,hipdat =get('hip_points_cells'),color_bar=colorbar,share="row",nan_color=reticulate::tuple(0.7, 0.7, 0.7, 1),size=as.integer(size), zoom=zoom,
                                         cmap=cmap,color_range=limits,label_text=title, return_plotter=T,interactive=F) ##disabling interactive mode because this causes RStudio to hang
  }
  #output plot as a .png image
  surf_plot$screenshot(filename=filename,transparent_bg = F)
}
############################################################################################################################
############################################################################################################################

#' @title Surface to volume
#'
#' @description Converts surface data to volumetric data (.nii file)
#'
#' @param surf_data A vector object containing the surface data, either in fsaverage5 or fsaverage6 space. It can only be one row of vertices (no cohort surface data matrix). 
#' @param filename A string object containing the desired name of the output .nii file (default is 'output.nii').
#'
#' @returns A .nii volume file
#' @examples
#' CTv = runif(20484,min=0, max=100);
#' surf_to_vol(CTv, filename = 'volume.nii')
#' @importFrom reticulate import
#' @export

##converting surface to volumetric data and exporting it as a .nii file

surf_to_vol=function(surf_data, filename)
  {
  #Check required python dependencies. If files missing:
  #Will prompt the user to get them in interactive session 
  #Will stop if it's a non-interactive session 
  cat("Checking for VertexWiseR system requirements ... ")
  check = VWRfirstrun(requirement="miniconda/brainstat")
  if (!is.null(check)) {return(check)} else {cat("\u2713 \n")}
  
  if (missing("filename")) {
    cat('No filename argument was given. The volume will be saved as "vol.nii" in R temporary directory (tempdir()).\n')
    filename=paste0(tempdir(),'/vol.nii')
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
    cat(filename)
  }

############################################################################################################################
############################################################################################################################
#' @title Decode surface data
#'
#' @description Correlates the significant clusters of an earlier vertex-wise analysis with a database of task-based fMRI and voxel-based morphometric statistical maps and associate them with relevant key words
#'
#' @details The \href{https://nimare.readthedocs.io/en/stable/index.html}{'NiMARE'} python module is used for the imaging decoding and is imported via the reticulate package. The function also downloads the \href{https://github.com/neurosynth/neurosynth-data}{'Neurosynth' database} in the package's inst/extdata direcotry (~8 Mb) for the analysis.
#'
#' @param surf_data a numeric vector with a length of 20484
#' @param contrast A string object indicating whether to decode the positive or negative mask ('positive' or 'negative')
#'
#' @returns A data.frame object listing the keywords and their Pearson's R values
#' @examples
#' CTv = rbinom(20484, 1, 0.001) 
#' decoding = decode_surf_data(CTv, 'positive');
#' head(decoding)
#' @importFrom reticulate import r_to_py
#' @export

##CT image decoding
decode_surf_data=function(surf_data,contrast="positive") 
{
  
  #Check required python dependencies. If files missing:
  #Will prompt the user to get them in interactive session 
  #Will stop if it's a non-interactive session
  cat("Checking for VertexWiseR system requirements ... ")
  check = VWRfirstrun(requirement="neurosynth")
  if (!is.null(check)) {return(check)} else {cat("\u2713 \n")}
  
  if(file.exists(system.file('extdata','neurosynth_dataset.pkl.gz', package='VertexWiseR'))==T)
  {
     
  ##checks length
    if(is.vector(surf_data)) {n_vert=length(surf_data)} else {n_vert=ncol(surf_data)}
    if(n_vert==20484) {template="fsaverage5"}
    else if (n_vert==81924) {stop("decoding of fsaverage6-space image is current not implemented, please resample the image to fsaverage5 space")} 
    else {stop("Only an surf_data vector with a length of 20484 (fsaverage5) is accepted")}

    #check contrast
    if(contrast != "positive" & contrast != "negative")  {stop("contrast has to be either positive or negative")} 
  
  cat("Converting and interpolating the surface data ... ")
  
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
  cat("\u2713 \n Correlating input image with images in the neurosynth database. This may take a while ... ")
  decoder = discrete$ROIAssociationDecoder(stat_nii)
  decoder$fit(neurosynth_dset)

  ##compiling the results
  decoder_df = data.matrix(decoder$transform())
  row.names(decoder_df)=gsub(pattern = "terms_abstract_tfidf__",x=row.names(decoder_df), replacement = "")
  result=data.frame(row.names(decoder_df),round(as.numeric(decoder_df),3))
  colnames(result)=c("keyword","r")
  result=result[order(-result$r),]
  cat("\u2713 \n")
  return(result)
  } 
}  


############################################################################################################################
############################################################################################################################
#' @title VertexWiseR system requirements installation
#'
#' @description Helps the user verify if VertexWisrR's system requirements are present and install them ('Miniconda', 'BrainStat' toolbox and libraries). If they are installed already, nothing will be overwritten. 
#'
#' @details VertexWiseR imports and makes use of the R package 'reticulate.' 'reticulate' is a package that allows R to borrow or translate Python functions into R. Using 'reticulate', the package calls functions from the 'BrainStat' Python module. For 'reticulate' to work properly with VertexWiseR, the latest version of 'Miniconda' needs to be installed with it — 'Miniconda' is a lightweight version of Python, specifically for use within 'RStudio'. Likewise, analyses of cortical surface require fsaverage templates as imported by 'BrainStat' The decode_surf_data() function also requires the 'Neurosynth' database to be downloaded.
#' @param requirement String that specifies a requirement to enquire about (for specific 'BrainStat' libraries: 'fsaverage5', 'fsaverage6', 'yeo_parcels'; for neurosynth database: "neurosynth"). Default is 'any' requirement and checks everything.
#' @param n_vert Numeric vector indicating the number of vertices of a given surface data so that only the required templates are asked for
#' @examples
#' VWRfirstrun()
#' @importFrom reticulate conda_binary py_module_available
#' @importFrom fs path_home
#' @importFrom methods is 
#' @importFrom utils menu read.csv write.csv
#' @export

VWRfirstrun=function(requirement="any", n_vert=0) 
{

#Loads the checklist which will store requirement checks results 
#and modify accordingly
  checklist = read.csv(paste0(system.file(package='VertexWiseR'),'/extdata/requirements_checklist.csv'), row.names = 1)
  
  
#First checks the n_vert argument. This ensures only the necessary fsaverage data is demanded:
  #are fsaverage5 templates in brainstat_data?
  if (n_vert==20484)
  {requirement='fsaverage5'}  
  #are fsaverage6 templates in brainstat_data?
  if  (n_vert==81924)
  {requirement='fsaverage6'}
  #is yeo parcellation data in brainstat_data?
  if (n_vert>0 & n_vert!=20484 & n_vert!=81924)
  {requirement='yeo_parcels'} 
  
  
if (interactive()==T) { #can only run interactively as it requires user's action
  
  #check if miniconda is installed
  if (is(tryCatch(reticulate::conda_binary(), error=function(e) e))[1] == 'simpleError')
  {
    checklist[1,]=FALSE;
    cat('Miniconda could not be found in the environment. \n')
    prompt = utils::menu(c("Yes", "No"), title=" Do you want miniconda to be installed now?")
    if (prompt==1) {reticulate::install_miniconda()}
    else { write.csv(checklist, file=paste0(system.file(package='VertexWiseR'),'/extdata/requirements_checklist.csv'), row.names = T)
      stop('VertexWiseR will not work properly without miniconda. reticulate::conda_list() should detect it on your system.\n\n')}
  } 
  else {checklist[1,]=TRUE} 
  
  #check if brainstat is installed
  if(!reticulate::py_module_available("brainstat") & requirement!="miniconda only") 
  {
    checklist[2,]=FALSE;
    cat('Brainstat could not be found in the environment. It is needed for vertex-wise linear models and the surface plotter to work.\n')
    prompt = utils::menu(c("Yes", "No"), title=" Do you want Brainstat to be installed now (~1.65 MB)? The NiMARE (~20.4 MB) and Brainspace (~84.2 MB) libraries are dependencies that will automatically be installed with it.")
    if (prompt==1){reticulate::py_install("brainstat",pip=TRUE)} 
    else {
      write.csv(checklist, file=paste0(system.file(package='VertexWiseR'),'/extdata/requirements_checklist.csv'), row.names = T)
      stop('VertexWiseR will not work properly without brainstat.\n\n')}
  } 
  else {checklist[2,]=TRUE} 
  
  
  #check if brainstat fsaverage/parcellation templates are installed (stops only if function needs it)
  if ((requirement=="any" | requirement=='fsaverage5')==T 
      & !file.exists(paste0(fs::path_home(),'/brainstat_data/surface_data/tpl-fsaverage/fsaverage5'))) 
  {     
    checklist[3,]=FALSE;
    cat('VertexWiseR could not find brainstat fsaverage5 templates in $home/brainstat_data/. They are needed if you want to analyse cortical surface in fsaverage5 space.')  
    prompt = utils::menu(c("Yes", "No"), title=" Do you want the fsaverage5 templates (~7.81 MB) to be downloaded now?")
    if (prompt==1){    
      brainstat.datasets.base=reticulate::import("brainstat.datasets.base", delay_load = TRUE)
      brainstat.datasets.base$fetch_template_surface("fsaverage5")
    } else if (requirement=='fsaverage5') {
      write.csv(checklist, file=paste0(system.file(package='VertexWiseR'),'/extdata/requirements_checklist.csv'), row.names = T)
      stop('VertexWiseR will not be able to analyse fsaverage5 data without the brainstat templates.\n\n')
    } else if (requirement=='any') {
      write.csv(checklist, file=paste0(system.file(package='VertexWiseR'),'/extdata/requirements_checklist.csv'), row.names = T)
      cat('VertexWiseR will not be able to analyse fsaverage5 data without the brainstat templates.\n\n')}
  } 
  else {checklist[3,]=TRUE} 
  
  if ((requirement=="any" | requirement=='fsaverage6')==T 
      & !file.exists(paste0(fs::path_home(),'/brainstat_data/surface_data/tpl-fsaverage/fsaverage6'))) 
  { 
    checklist[4,]=FALSE;
    cat('VertexWiseR could not find brainstat fsaverage6 templates in $home/brainstat_data/. They are needed if you want to analyse cortical surface in fsaverage6 space.')
   prompt = utils::menu(c("Yes", "No"), title=" Do you want the fsaverage6 templates (~31.2 MB) to be downloaded now?")
    
     if (prompt==1)
      { brainstat.datasets.base=reticulate::import("brainstat.datasets.base", delay_load = TRUE)
        brainstat.datasets.base$fetch_template_surface("fsaverage6")
      } else if (requirement=='fsaverage6') { 
        write.csv(checklist, file=paste0(system.file(package='VertexWiseR'),'/extdata/requirements_checklist.csv'), row.names = T)
        stop('VertexWiseR will not be able to analyse fsaverage6 data without the brainstat templates.\n\n')
      } else if (requirement=="any") {
        write.csv(checklist, file=paste0(system.file(package='VertexWiseR'),'/extdata/requirements_checklist.csv'), row.names = T)
        cat('VertexWiseR will not be able to analyse fsaverage6 data without the brainstat templates.\n\n')}
  } 
  else {checklist[4,]=TRUE} 
  
  if ((requirement=="any" | requirement=='fsaverage6' | requirement=='fsaverage5' | requirement=='yeo_parcels')==T 
      & !file.exists(paste0(fs::path_home(),'/brainstat_data/parcellation_data/')))
  {
     checklist[5,]=FALSE;
      cat('VertexWiseR could not find brainstat yeo parcellation data in $home/brainstat_data/. They are fetched by default by brainstat for vertex-wise linear models to run and cannot be ignored.')
      prompt = utils::menu(c("Yes", "No"), title=" Do you want the yeo parcellation data (~1.01 MB) to be downloaded now?")
      if (prompt==1){    
        brainstat.datasets.base=reticulate::import("brainstat.datasets.base", delay_load = TRUE)
        try(brainstat.datasets.base$fetch_parcellation(template="fsaverage",atlas="yeo", n_regions=7), silent=T)}  
        else if  (requirement=='fsaverage6' | requirement=='fsaverage5' | requirement=='yeo_parcels') 
        {
          write.csv(checklist, file=paste0(system.file(package='VertexWiseR'),'/extdata/requirements_checklist.csv'), row.names = T)
          stop('VertexWiseR will not be able to analyse cortical data without the parcellation data.\n\n')}
        else if (requirement=="any") 
        {
          write.csv(checklist, file=paste0(system.file(package='VertexWiseR'),'/extdata/requirements_checklist.csv'), row.names = T)
          cat('VertexWiseR will not be able to analyse cortical data without the parcellation data.\n\n')}
  }
  else {checklist[5,]=TRUE} 
  
  #Check if neurosynth database is present and download
  if ((requirement=="any" | requirement=='neurosynth')==T 
      & !file.exists(system.file('extdata','neurosynth_dataset.pkl.gz', package='VertexWiseR')))
  {
      checklist[6,]=FALSE;
      cat("\nneurosynth_dataset.pkl.gz is not detected in the package's external data directory (/inst/extdata). It is needed to be able to run decode_surf_data(). It can be downloaded from the github VertexWiseR directory.\n")
      prompt = utils::menu(c("Yes", "No"), title="Do you want the neurosynth database (7.5 MB) to be downloaded now?")
        if (prompt==1) {
          
          #function to check if url exists
          #courtesy of Schwarz, March 11, 2020, CC BY-SA 4.0:
          #https://stackoverflow.com/a/60627969
          valid_url <- function(url_in,t=2){
            con <- url(url_in)
            check <- suppressWarnings(try(open.connection(con,open="rt",timeout=t),silent=T)[1])
            suppressWarnings(try(close.connection(con),silent=T))
            ifelse(is.null(check),TRUE,FALSE)}
          
          #Check if URL works and avoid returning error but only print message as requested by CRAN:
          url="https://raw.githubusercontent.com/CogBrainHealthLab/VertexWiseR/main/inst/extdata/neurosynth_dataset.pkl.gz"
          if(valid_url(url)) {
            download.file(url="https://raw.githubusercontent.com/CogBrainHealthLab/VertexWiseR/main/inst/extdata/neurosynth_dataset.pkl.gz",destfile = paste0(system.file(package='VertexWiseR'),'/extdata/neurosynth_dataset.pkl.gz'))
          } else { 
            return("The neurosynth database (neurosynth_dataset.pkl.gz) failed to be downloaded from the github VertexWiseR directory. Please check your internet connection. Alternatively, you may visit https://github.com/CogBrainHealthLab/VertexWiseR/tree/main/inst/extdata and download the object manually.") #ends function
          } 
          
        #if user refuses, stops if required, just returns a message if optionnal at this stage
        } else if (requirement=="neurosynth") {
          write.csv(checklist, file=paste0(system.file(package='VertexWiseR'),'/extdata/requirements_checklist.csv'), row.names = T)
          stop("\ndecode_surf_data() can only work with the neurosynth database.\n") }       
         else if (requirement=="any") {
           write.csv(checklist, file=paste0(system.file(package='VertexWiseR'),'/extdata/requirements_checklist.csv'), row.names = T)
           cat("\ndecode_surf_data() can only work with the neurosynth database.\n")}
  }
  else {checklist[6,]=TRUE} 
  
  #update checklist csv
  write.csv(checklist, file=paste0(system.file(package='VertexWiseR'),'/extdata/requirements_checklist.csv'), row.names = T)
  
} 
else #if not interactive and any required file is missing, the script requires the user to run VWR interactively
{ 
  #creates the following object to warn upper functions that it's a non-interactive session when files are missing
  non_interactive="If system requirements are missing, VWRfirstrun() can only be run in an interactive R session to check for system requirements and to install them.";
  
  #miniconda missing?
  if (checklist[1,1] == FALSE) 
  {return(non_interactive)}
  
  #brainstat missing?
  if (checklist[2,1] == FALSE)
  {return(non_interactive)}
  
  #fsaverage5 missing?
  if ((requirement=="any" | requirement=='fsaverage5')==T 
      & checklist[3,1] == FALSE) 
  {return(non_interactive)}
  
  #fsaverage6 missing?
  if ((requirement=="any" | requirement=='fsaverage6')==T 
      & checklist[4,1] == FALSE)  
  {return(non_interactive)}
  
  #yeo parcels missing?
  if ((requirement=="any" | requirement=='fsaverage6' | requirement=='fsaverage5' | requirement=='yeo_parcels')==T 
      & checklist[5,1] == FALSE) 
  {return(non_interactive)}
  
   #neurosynth data missing?
   if ((requirement=="any" | requirement=='neurosynth')==T 
       & checklist[6,1] == FALSE) 
     {return(non_interactive)}
      
}

}

