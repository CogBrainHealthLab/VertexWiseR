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
#' @param surf_data A N x M matrix object containing the surface data (N row for each subject, M for each vertex), in fsaverage5 (20484 vertices), fsaverage6 (81924 vertices) or hippocampal (14524 vertices) space. See also Hipvextract() or SURFvextract() output format. 
#' @param FWHM A numeric vector object containing the desired smoothing width in mm 
#' @param VWR_check A boolean object specifying whether to check and validate system requirements. Default is TRUE.
#'
#' @returns A matrix object with smoothed vertex-wise values
#' @examples
#' surf_data = readRDS(file = url(paste0("https://github.com",
#'"/CogBrainHealthLab/VertexWiseR/blob/main/inst/demo_data/",
#'"FINK_Tv_ses13.rds?raw=TRUE")))[1:3,]
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
  if (VWR_check == TRUE){
    message("Checking for VertexWiseR system requirements ... ")
    check = VWRfirstrun(requirement="miniconda only")
    if (!is.null(check)) {return(check)} else {message("\u2713 \n")}
  } else if(interactive()==FALSE) { return(message('Non-interactive sessions need requirement checks'))}
  
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
#' For hippocampal data, the function currently works with the "bigbrain" atlas integrated in 'HippUnfold.' See also \doi{doi:10.1016/j.neuroimage.2019.116328}.
#'
#' @param surf_data A N x M matrix object containing the surface data (N row for each subject, M for each vertex), in fsaverage5 (20484 vertices), fsaverage6 (81924 vertices) or hippocampal (14524 vertices) space. See also Hipvextract() or SURFvextract() output format. 
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
  
  
  #mapping fsaverage6 space vertice to atlas (Nx81924 vertices)
  if(max(dim(t(surf_data)))==81924) 
  {
    #load atlas mapping surf_data
    ROImap <- get('ROImap_fs6')
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
  
  #mapping hippocampal space vertice to atlas (Nx14524 vertices)
  if(max(dim(t(surf_data)))==14524) 
  {
    #load atlas mapping surf_data
    ROImap <- get('ROImap_HIP')
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
  } else if(is.vector(parcel_data)==TRUE) #if parcel_data is a vector
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
#' @param surf_data A N x M matrix object containing the surface data (N row for each subject, M for each vertex), in fsaverage5 (20484 vertices)  space. See also SURFvextract() output format. 
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
#' @param A N x M matrix object containing the surface data (N row for each subject, M for each vertex), in fsaverage6 (81924 vertices) space. See also SURFvextract() output format.  
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
#' @param filename A string object containing the desired name of the output .png. Default is 'plot.png' in the R temporary directory (tempdir()).
#' @param title A string object for setting the title in the plot. Default is none. For titles that too long to be fully displayed within the plot, we recommend splitting them into multiple lines by inserting "\\n".
#' @param surface A string object containing the name of the type of cortical surface background rendered. Possible options include "white", "smoothwm","pial" and "inflated" (default). The surface parameter is ignored for hippocampal surface data.
#' @param cmap A string object specifying the name of an existing colormap or a vector of hexadecimal color codes to be used as a custom colormap. The names of existing colormaps are listed in the \href{https://matplotlib.org/stable/gallery/color/colormap_reference.html}{'Matplotlib' plotting library}. 
#' 
#' Default cmap is set to `"Reds"` for positive values, `"Blues_r"` for negative values and `"RdBu"` when both positive and negative values exist. 
#' @param limits A combined pair of numeric vector composed of the lower and upper color scale limits of the plot. If the limits are specified, the same limits will be applied to all subplots. When left unspecified, the same symmetrical limits c(-max(abs(surf_dat),max(abs(surf_dat))) will be used for all subplots. If set to NULL, each subplot will have its own limits corresponding to their min and max values
#' @param colorbar A logical object stating whether to include a color bar in the plot or not (default is TRUE).
#' @param size A combined pair of numeric vector indicating the image dimensions (width and height in pixels). Default is c(1920,400) for whole-brain surface and c(400,200) for hippocampal surface.
#' @param zoom A numeric value for adjusting the level of zoom on the figures. Default is 1.25 for whole-brain surface and 1.20 for hippocampal surface.
#' @param show.plot.window A logical object to determine if the generated plot is to be shown within RStudio's plot window
#' @param VWR_check A boolean object specifying whether to check and validate system requirements. Default is TRUE.
#'
#' @returns Outputs the plot as a .png image
#' @examples
#' results = runif(20484,min=0, max=1);
#' plot_surf(surf_data = results, 
#' filename=paste0(tempdir(),"/output.png"),
#' title = 'Cortical thickness', 
#' surface = 'inflated', cmap = 'Blues',
#' VWR_check=FALSE)
#' @importFrom reticulate tuple import np_array source_python
#' @importFrom grDevices col2rgb
#' @importFrom png readPNG
#' @importFrom grid grid.raster
#' @export

plot_surf=function(surf_data, filename, title="",surface="inflated",cmap,limits, colorbar=TRUE, size, zoom, show.plot.window=TRUE,VWR_check=TRUE)
{
  #Check required python dependencies. If files missing:
  #Will prompt the user to get them in interactive session 
  #Will stop if it's a non-interactive session 
  if (VWR_check == TRUE){
    message("Checking for VertexWiseR system requirements ...")
    check = VWRfirstrun(n_vert=max(dim(t(surf_data))))
    if (!is.null(check)) {return(check)} else {message("\u2713 \n")}
  } else if(interactive()==FALSE) { return(message('Non-interactive sessions need requirement checks'))}
  
  if (missing("filename")) {
    message('No filename argument was given. The plot will be saved as "plot.png" in R temporary directory (tempdir()).\n')
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
    if(range(surf_data,na.rm = TRUE)[1]>=0)  {cmap="Reds"}
    else if (range(surf_data,na.rm = TRUE)[2]<=0)  {cmap="Blues_r"}
    else  {cmap="RdBu_r"}  
  }
	
  #custom cmap— if a vector of hex color codes is specified
  if(inherits(cmap,"colors")==TRUE)
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
   maxlimit=max(abs(range(surf_data,na.rm = TRUE)))
    if(missing("limits")) 
    {
      if(range(surf_data,na.rm = TRUE)[1]>=0) {limits=reticulate::tuple(0,range(surf_data,na.rm = TRUE)[2])} ##if image contains all positive values
      else if(range(surf_data,na.rm = TRUE)[2]<=0) {limits=reticulate::tuple(range(surf_data,na.rm = TRUE)[1],0)} ##if image contains all negative values
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
    left=brainstat.datasets$fetch_template_surface(template, join=FALSE, layer=surface)[1]
    right=brainstat.datasets$fetch_template_surface(template, join=FALSE, layer=surface)[2]
    
    #default cortical size and zoom parametes
    if(missing("size")) { size=c(1920,rows*400)}
    if(missing("zoom")) { zoom=1.25 }
    
    surf_plot=brainspace.plotting$plot_hemispheres(left[[1]], right[[1]],  array_name=reticulate::np_array(surf_data),cmap=cmap, 
                                                size=reticulate::tuple(as.integer(size)),nan_color=reticulate::tuple(0.7, 0.7, 0.7, 1),
                                                return_plotter=TRUE,background=reticulate::tuple(as.integer(c(1,1,1))),zoom=zoom,color_range=limits,
                                                label_text=title,interactive=FALSE, color_bar=colorbar,  transparent_bg=FALSE)  ##disabling interactive mode because this causes RStudio to hang
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
                                         cmap=cmap,color_range=limits,label_text=title, return_plotter=TRUE,interactive=FALSE) ##disabling interactive mode because this causes RStudio to hang
  }
  #output plot as a .png image
  surf_plot$screenshot(filename=filename,transparent_bg = FALSE)
  if(show.plot.window==T)
	{
		img=png::readPNG(source =filename)
	  grid::grid.raster(img)
	}
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
    check = VWRfirstrun(requirement="miniconda/brainstat")
    if (!is.null(check)) {return(check)} else {message("\u2713 \n")}
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
    if (!is.null(check)) {return(check)} else {message("\u2713 \n")}
  } else if(interactive()==FALSE) { return(message('Non-interactive sessions need requirement checks'))}
  
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
#' @title VertexWiseR system requirements installation
#'
#' @description Helps the user verify if VertexWisrR's system requirements are present and install them ('Miniconda', 'BrainStat' toolbox and libraries). If they are installed already, nothing will be overwritten. 
#'
#' @details VertexWiseR imports and makes use of the R package 'reticulate.' 'reticulate' is a package that allows R to borrow or translate Python functions into R. Using 'reticulate', the package calls functions from the 'BrainStat' Python module. For 'reticulate' to work properly with VertexWiseR, 'Miniconda' needs to be installed with it — 'Miniconda' is a lightweight version of Python, specifically for use within 'RStudio'. If for a reason Miniconda cannot be installed, the function gives the choice to install a reticulate-suitable Python environment.
#' Vertex-wise statistical analyses of cortical surface require fsaverage and parcellation templates as imported by default in 'BrainStat'. 
#' The decode_surf_data() function also requires the 'Neurosynth' database to be downloaded.
#' @param requirement String that specifies a requirement to enquire about (for specific 'BrainStat' libraries: 'fsaverage5', 'fsaverage6', 'yeo_parcels'; for neurosynth database: "neurosynth"). Default is 'any' requirement and checks everything.
#' @param n_vert Numeric vector indicating the number of vertices of a given surface data so that only the required templates are asked for
#' @param promptless A boolean object specifying whether to prompt the user for action when system requirements are missing. If TRUE, VWRfirstrun() will simply inform of what is missing and will not prompt for action. Default is FALSE.
#' @return No returned value in interactive session. In non-interactive sessions, a string object informing that system requirements are missing.
#' @examples
#' VWRfirstrun()
#' @importFrom reticulate conda_binary py_module_available
#' @importFrom fs path_home
#' @importFrom methods is 
#' @importFrom utils menu
#' @importFrom rappdirs user_data_dir
#' @export

VWRfirstrun=function(requirement="any", n_vert=0, promptless=FALSE) 
{
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
  
# If custom installation paths have been defined by the user, source
# them from their define directory:
Renvironpath=paste0(tools::R_user_dir(package='VertexWiseR'),'/.Renviron')
if (file.exists(Renvironpath)) {readRenviron(Renvironpath)}

#default time limit to download is 60s which can be too short:
options(timeout=500); #set to 500s instead

if (interactive()==TRUE & promptless==FALSE) { 
  #can only run interactively as it requires user's action
  
  ##################################################################
  #check if miniconda (or suitable python environment) is installed
  if(is.null(reticulate::py_discover_config()) |
     is(tryCatch(reticulate::py_discover_config(), error=function(e) e))[1] == 'simpleError')
  {
    prompt = utils::menu(c("Yes, install Miniconda (Recommended)",
                           "Yes, install Python (if Miniconda could not be installed)", 
                           "No"), 
                         title="Miniconda or a suitable version of Python for reticulate could not be found in the environment. \n Do you want Miniconda or Python to be installed now?")
    
    if (prompt==1) #Install Miniconda
    {
      
      
      #define default path 
      #mimicking reticulate's non-exported miniconda_path_default()
      if (.Platform$OS.type == "mac") {
        # on macOS, different path for arm64 miniconda
        if (Sys.info()[["machine"]] == "arm64") {defaultpath="~/Library/r-miniconda-arm64"} else {defaultpath="~/Library/r-miniconda"}
        } else {
      # otherwise, use rappdirs default
      root <- normalizePath(rappdirs::user_data_dir(), winslash = "/",
                            mustWork = FALSE); 
      defaultpath=file.path(root, "r-miniconda")}
      
      
      
      #give the choices to specify path
      choice = utils::menu(c("Default", "Custom"), 
                           title=paste0("Miniconda's default installation path is ", defaultpath,". Type \"1\" or \"Default\" if you want to install Miniconda in the default Path. \nYou can, alternatively, type your own path (note that the Miniconda installer does not support paths containing spaces)."))
        
        if (choice==1) #Install Miniconda within default path
          { reticulate::install_miniconda()}
      
        else {        #Install Miniconda within custom path
        
        userpath <- readline("Enter the full path:")
        
        #will store path in .Renviron in tools::R_user_dir() 
        #location specified by CRAN, create it if not existing:
        envpath=tools::R_user_dir(package='VertexWiseR')
        if (!dir.exists(envpath)) {dir.create(envpath) }
        #make .Renviron file there and set conda/python paths:
        renviron_path <- file.path(envpath, ".Renviron")
        env_vars <- c(
          paste0('RETICULATE_PYTHON="',userpath,'python.exe"'),
          paste0('RETICULATE_MINICONDA_PATH="',userpath,'"'))
        # Write to the .Renviron file
        cat(paste(env_vars, collapse = "\n"), file = renviron_path, 
            sep = "\n", append = TRUE)
        #now everytime VWRfirstrun() is called, the .Renviron file 
        #is read and the custom path accessed:
        readRenviron(renviron_path)
        message(paste0("FYI, your custom path has been saved in ", 
                       renviron_path))
        
        #Install miniconda in the new path
        message('Installing miniconda ...')
        reticulate::install_miniconda() #uses miniconda_path() which relies on RETICULATE_MINICONDA_PATH defined above
        reticulate::py_install("numpy==1.26.4") 
        #posterior numpy versions break python functions
        }
      
  }
    else if (prompt==2) #Install Python instead of Miniconda
    { 
    
      #give the choices to specify path
      choice = utils::menu(c("Default", "Custom"), 
                           title=paste0("Reticulate's Python default installation is done within 'r-reticulate' in your OS libraries, via install_python(). Type \"1\" or \"Default\" if you want the default installation. \nYou can, alternatively, specify the path where Python will be installed."))
      
      if (choice==1) #Install Python within default path
        {reticulate::install_python()}
      
      else {        #Install Python within custom path
      
      userpath <- readline("Enter the full path:")
      
      message(paste('Installing Python in',userpath,'...\n'))
      
      python_custominstall(userpath)
      
      #Read new python enviroment
      Renvironpath=paste0(tools::R_user_dir(package='VertexWiseR'),'/.Renviron')
      if (file.exists(Renvironpath)) {readRenviron(Renvironpath)}
      
      message('A specific version of the Numpy package (<= 1.26.4) is more stable for analyses, it will now be installed in the Python libraries along with its dependencies.\n')
      system('pip install numpy==1.26.4') #pip instead of install_py as it will use a virtual environment
      #posterior numpy versions break python functions
      }
  }
  else { stop('VertexWiseR will not work properly without Miniconda or a suitable version of Python for reticulate.\n\n')}
  } 
  
  
  
  ##################################################################
  ###if Python installation can be found by reticulate, check numpy version
  if(!is.null(reticulate::py_discover_config()) |
     !is(tryCatch(reticulate::py_discover_config(), error=function(e) e))[1] == 'simpleError')
  { 
    
    # check if numpy is in the default python environment
    listpackages=reticulate::py_list_packages()
    #check numpy version
    numpyv=listpackages$version[which(listpackages$package=="numpy")] 
    
    if (length(numpyv) > 0) #if numpy exists at all 
    {
        #warn to install 1.26.2 if current version is superior  
        numpyv=gsub("\\.", "", numpyv);
        
        if (as.numeric(numpyv)/(10 ^ (nchar(numpyv) - 1)) > 1.264)
        { 
         warning("The current Python environment's Numpy package is version > 1.26.4. This may cause issues with this package.")
        }
    }
   }
  
  ##################################################################
  ###check if brainstat is installed
    
  #check everywhere brainstat could be installed
  if( !reticulate::py_module_available("brainstat") 
      & requirement!="miniconda only") 
  {
    prompt = utils::menu(c("Yes", "No"), title="The Brainstat package could not be found in your Python/Conda environment. It is needed for vertex-wise linear models and the surface plotter to work. \n Do you want Brainstat to be installed now (~1.65 MB)? The NiMARE (~20.4 MB) and Brainspace (~84.2 MB) libraries and other BrainStat dependencies will automatically be installed within your Python library.")
    if (prompt==1)
    {	
      
      #if miniconda exists, use reticulate to install
      if (tryCatch(file.exists(reticulate::conda_binary()), error = function(e) FALSE)) {
        reticulate::py_install("brainstat",pip=TRUE) 
      }
      else { #if only Python, install via pip
        system('pip install brainstat') 
        #install_py would use make virtual environment
      }
      
    } 
    else {
      stop('VertexWiseR will not work properly without brainstat.\n\n')}
  } 
  
  ##################################################################
  #check if brainstat fsaverage/parcellation templates are installed (stops only if function needs it)
  if ((requirement=="any" | requirement=='fsaverage5')==TRUE 
      & !file.exists(paste0(fs::path_home(),'/brainstat_data/surface_data/tpl-fsaverage/fsaverage5'))) 
  {     
    prompt = utils::menu(c("Yes", "No"), title="VertexWiseR could not find brainstat fsaverage5 templates in $HOME/brainstat_data/. They are needed if you want to analyse cortical surface in fsaverage5 space. \n  Do you want the fsaverage5 templates (~7.81 MB) to be downloaded now?")
    if (prompt==1){    
      brainstat.datasets.base=reticulate::import("brainstat.datasets.base", delay_load = TRUE)
      brainstat.datasets.base$fetch_template_surface("fsaverage5")
    } else if (requirement=='fsaverage5') {
      stop('VertexWiseR will not be able to analyse fsaverage5 data without the brainstat templates.\n\n')
    } else if (requirement=='any') {
      warning('VertexWiseR will not be able to analyse fsaverage5 data without the brainstat templates.\n\n')}
  } 
  
  if ((requirement=="any" | requirement=='fsaverage6')==TRUE 
      & !file.exists(paste0(fs::path_home(),'/brainstat_data/surface_data/tpl-fsaverage/fsaverage6'))) 
  { 
   prompt = utils::menu(c("Yes", "No"), title="VertexWiseR could not find brainstat fsaverage6 templates in $home/brainstat_data/. They are needed if you want to analyse cortical surface in fsaverage6 space. \n Do you want the fsaverage6 templates (~31.2 MB) to be downloaded now?")
    
     if (prompt==1)
      { brainstat.datasets.base=reticulate::import("brainstat.datasets.base", delay_load = TRUE)
        brainstat.datasets.base$fetch_template_surface("fsaverage6")
      } else if (requirement=='fsaverage6') { 
        stop('VertexWiseR will not be able to analyse fsaverage6 data without the brainstat templates.\n\n')
      } else if (requirement=="any") {
        warning('VertexWiseR will not be able to analyse fsaverage6 data without the brainstat templates.\n\n')}
  } 
  
  if ((requirement=="any" | requirement=='fsaverage6' | requirement=='fsaverage5' | requirement=='yeo_parcels')==TRUE 
      & !file.exists(paste0(fs::path_home(),'/brainstat_data/parcellation_data/')))
  {
      prompt = utils::menu(c("Yes", "No"), title="VertexWiseR could not find brainstat yeo parcellation data in $home/brainstat_data/. They are fetched by default by brainstat for vertex-wise linear models to run and cannot be ignored. \n Do you want the yeo parcellation data (~1.01 MB) to be downloaded now?")
      if (prompt==1){    
        brainstat.datasets.base=reticulate::import("brainstat.datasets.base", delay_load = TRUE)
        try(brainstat.datasets.base$fetch_parcellation(template="fsaverage",atlas="yeo", n_regions=7), silent=TRUE)}  
        else if  (requirement=='fsaverage6' | requirement=='fsaverage5' | requirement=='yeo_parcels') 
        {
          stop('VertexWiseR will not be able to analyse cortical data without the parcellation data.\n\n')}
        else if (requirement=="any") 
        {
          warning('VertexWiseR will not be able to analyse cortical data without the parcellation data.\n\n')
        }
  }
  
  #Check if neurosynth database is present and download
  if ((requirement=="any" | requirement=='neurosynth')==TRUE 
      & !file.exists(system.file('extdata','neurosynth_dataset.pkl.gz', package='VertexWiseR')))
  {
      prompt = utils::menu(c("Yes", "No"), title=paste0(
        "\nneurosynth_dataset.pkl is not detected inside VertexWiseR's installed package directory (", 
        system.file('extdata','neurosynth_dataset.pkl.gz', package='VertexWiseR'), 
        "). It is needed to be able to run decode_surf_data(). It can be downloaded from the github VertexWiseR directory.\n\nDo you want the neurosynth database (7.5 MB) to be downloaded now?"))
        if (prompt==1) {
          
          #function to check if url exists
          #courtesy of Schwarz, March 11, 2020, CC BY-SA 4.0:
          #https://stackoverflow.com/a/60627969
          valid_url <- function(url_in,t=2){
            con <- url(url_in)
            check <- suppressWarnings(try(open.connection(con,open="rt",timeout=t),silent=TRUE)[1])
            suppressWarnings(try(close.connection(con),silent=TRUE))
            ifelse(is.null(check),TRUE,FALSE)}
          
          #Check if URL works and avoid returning error but only print message as requested by CRAN:
          url="https://raw.githubusercontent.com/CogBrainHealthLab/VertexWiseR/main/inst/extdata/neurosynth_dataset.pkl.gz"
          if(valid_url(url)) {
            download.file(url="https://raw.githubusercontent.com/CogBrainHealthLab/VertexWiseR/main/inst/extdata/neurosynth_dataset.pkl.gz",destfile = paste0(system.file(package='VertexWiseR'),'/extdata/neurosynth_dataset.pkl.gz'))
          } else { 
            warning("The neurosynth database (neurosynth_dataset.pkl.gz) failed to be downloaded from the github VertexWiseR directory. Please check your internet connection. Alternatively, you may visit https://github.com/CogBrainHealthLab/VertexWiseR/tree/main/inst/extdata and download the object manually.") #ends function
          } 
          
        #if user refuses, stops if required, just returns a message if optionnal at this stage
        } else if (requirement=="neurosynth") {
          stop("\ndecode_surf_data() can only work with the neurosynth database.\n") }       
         else if (requirement=="any") {
           warning("\ndecode_surf_data() can only work with the neurosynth database.\n")}
  }
  
} 

#############################################################
#If the session is non-interactive and any required file is missing, the script will stop and require the user to run VWRfirstrun() interactively.
#non-interactive checks proceed when VWRfirstrun() is called by a function with the argument VWR_check defined.

#If the session is interactive but has the promptless option set as TRUE, the same checks will be run to inform of what is missing while not prompting for any action from the user.

else if (exists("VWR_check") | 
         (interactive()==TRUE & promptless==TRUE)) 
{ 
  #creates the 'non_interactive' message to inform upper functions that the non-interactive session cannot proceed when files are missing
  non_interactive="Run VWRfirstrun() in an interactive R session to check for the missing system requirements and to install them."
    
  #miniconda or python missing?
  if(is.null(reticulate::py_discover_config()) |
     is(tryCatch(reticulate::py_discover_config(), error=function(e) e))[1] == 'simpleError')
  {
    missingobj="Miniconda or a suitable version of Python for reticulate could not be found in the environment.\n"
    
    if (interactive()==FALSE)
    { non_interactive=paste0(missingobj,non_interactive)
      return(non_interactive)
    } else {message(missingobj)}
  } 
  
  #brainstat missing
  if (!reticulate::py_module_available("brainstat") & requirement!="miniconda only") 
  {
    missingobj="The Brainstat package could not be found in your Python/Conda environment. It is needed for vertex-wise linear models and the surface plotter to work.\n";
    
    if (interactive()==FALSE)
    { non_interactive=paste0(missingobj,non_interactive)
      return(non_interactive)
    } else {message(missingobj)}
  } 
  
  
  #fsaverage5 missing
  if ((requirement=="any" | requirement=='fsaverage5')==TRUE & !file.exists(paste0(fs::path_home(),'/brainstat_data/surface_data/tpl-fsaverage/fsaverage5'))) 
  {
    missingobj="VertexWiseR could not find brainstat fsaverage5 templates in the brainstat_data/ directory. They are needed if you want to analyse cortical surface in fsaverage5 space.\n";
    
    if (interactive()==FALSE)
    { non_interactive=paste0(missingobj,non_interactive)
      return(non_interactive)
    } else {message(missingobj)}
  } 
  
  #fsaverage6 missing
  if ((requirement=="any" | requirement=='fsaverage6')==TRUE & !file.exists(paste0(fs::path_home(),'/brainstat_data/surface_data/tpl-fsaverage/fsaverage6')))  
  {
    missingobj="VertexWiseR could not find brainstat fsaverage6 templates in the brainstat_data/ directory. They are needed if you want to analyse cortical surface in fsaverage6 space.\n";
    
    if (interactive()==FALSE)
    { non_interactive=paste0(missingobj,non_interactive)
    return(non_interactive)
    } else {message(missingobj)}
  } 
  
  
  #yeo parcels missing
  if ((requirement=="any" | requirement=='fsaverage6' | requirement=='fsaverage5' | requirement=='yeo_parcels')==TRUE 
      & !file.exists(paste0(fs::path_home(),
                          '/brainstat_data/parcellation_data/'))) 
  {
    missingobj="VertexWiseR could not find brainstat yeo parcellation data in the brainstat_data/ directory. They are fetched by default by brainstat for vertex-wise linear models to run and cannot be ignored.\n";
    
    if (interactive()==FALSE)
    { non_interactive=paste0(missingobj,non_interactive)
    return(non_interactive)
    } else {message(missingobj)}
  } 
  
  #neurosynth data missing
  if ((requirement=="any" | requirement=='neurosynth')==TRUE & !file.exists(system.file('extdata','neurosynth_dataset.pkl.gz', package='VertexWiseR'))) 
  {
    missingobj=paste0("neurosynth_dataset.pkl is not detected inside VertexWiseR's installed package directory (", system.file('extdata','neurosynth_dataset.pkl.gz', package='VertexWiseR'), "). It is needed to be able to run decode_surf_data().\n");
    
    if (interactive()==FALSE)
    { non_interactive=paste0(missingobj,non_interactive)
    return(non_interactive)
    } else {message(missingobj)}
  } 
  
  #If nothing is missing, missingobj will not have been created
  #thus allowing the function to inform that no (specified) requirements are missing
  if (promptless==TRUE) { if(!exists('missingobj'))
  { message('No system requirements are missing. \u2713') }}
  
}

}

