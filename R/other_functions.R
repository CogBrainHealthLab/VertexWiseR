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
  R = chol2inv(Qr[p1, p1, drop = FALSE])  
  if(is.matrix(mod$coefficients))
  {
    rss = colSums(r^2)
    resvar = rss/rdf 
    se = (sqrt(diag(R) %*% t(resvar)))[row,]
    est = mod$coefficients[row,]
  } else {
    rss = sum(r^2)
    resvar = rss/rdf
    se = (sqrt(diag(R) * resvar))[row]
    est = mod$coefficients[row]
  }
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
  edgelist1 = edgelist[rowSums(matrix(edgelist %in% vert,ncol=2)) == 2, , drop = FALSE]
  
  if(length(edgelist1)>2) #if at least 2 edges are identified
  {
    #extracting cluster-related info from list of non-zero edges
    #graph_from_data_frame() works a lot faster with character data
    com=igraph::components(igraph::graph_from_edgelist(edgelist1, directed = FALSE))
    clust.size=com$csize
    
    #cluster mappings
    clust.idx=which(com$csize>1) # need to remove clusters with only a single vertex
    clust.map=rep(NA,n_vert)
    clust.map[1:max(edgelist1)]=match(com$membership,clust.idx) #graph_from_edgelist does not return com$membership for each of the n_vert vertices
  
    clust.size=clust.size[clust.idx]
  
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
  brainstat.datasets=reticulate::import("brainstat.datasets", delay_load = TRUE)  
  brainstat.mesh.utils=reticulate::import("brainstat.mesh.utils", delay_load = TRUE)
  
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

############################################################################################################################
############################################################################################################################
  
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
  brainstat.datasets=reticulate::import("brainstat.datasets", delay_load = TRUE)  
  brainspace.mesh.mesh_elements=reticulate::import("brainspace.mesh.mesh_elements", delay_load = TRUE)
  
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
  
  #if contrast or model is a data.frame with 1 column, the variable needs to be flattened for the inherits() checks to work properly
  if(inherits(contrast,"data.frame")==TRUE & NCOL(contrast)==1){contrast=contrast[[1]]}
  if(inherits(model,"data.frame")==TRUE & NCOL(model)==1){model=model[[1]]}

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
      else  {stop("contrast is not contained within model")}
    } else
    {
      if(identical(as.numeric(contrast),as.numeric(model)))  {colno=1}
      else  {stop("contrast is not contained within model")}
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
      if(inherits(model[,column],"character")==TRUE | inherits(model[,column],"factor")==TRUE)
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
    
    if(inherits(model,"character")==TRUE | inherits(model,"factor")==TRUE)
    {
      if(length(unique(model))==2)
      {
        message(paste("The model variable is binary and will be recoded such that ",unique(model)[1],"=0 and ",unique(model)[2],"=1 for the analysis\n",sep=""))
        
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
  
############################################################################################################################
############################################################################################################################

#function to automatically read the surface matrix from a list object outputted by the extracter function (containing both the surface matrix and the list of subjects)
#if it is a string path to the rds, loads the file first

get_surf_obj=function(surf_data)
{
  #if surface_data is a path to an object, reads it
  if(inherits(surf_data,'character')==TRUE)
  { #if fails to read the path to the RDS, return error
    if (is(tryCatch(readRDS(surf_data), error=function(e) e))[1] == 'simpleError') 
  {stop('The surf_data given is a string and was therefore assumed to be a path to a \'.rds\' surface data file. The path failed to be accessed by readRDS().')} 
   else 
   {surf_data=readRDS(file=surf_data)}
  }
  
  #if surf_data contains subject list only read surface object
  if(inherits(surf_data,'list')==TRUE)
  { if ('surf_obj' %in% names(surf_data))
    {surf_data=as.matrix(surf_data$surf_obj)} 
    else {stop('The surf_data given is a list, but the package does not know which element in the list is meant to be the surface matrix. Please name the element "surf_obj" or enter the matrix as surf_data.'
    )}
  }
  return(surf_data)
}
 
############################################################################################################################
############################################################################################################################
  
#' @title Database check for SubCortexMesh subcortical surface analyses
#'
#' @description Function to check whether the external subcortical template database has been downloaded in order to run plotting, smoothing and modelling functions on surface data derived from SubCortexMesh processing. It assists downloading the data from the github online repository. It is only triggered when a user inputs surf_data with the same number of vertices as SubCortexMesh-derived objects.
#' @details The database is downloaded and stored in the path indicated by system.file('extdata', package='VertexWiseR'). 
#' @param template A string indicating the subcortical segmentation scheme. Default and currently only option is 'aseg'.
#' @examples
#' scm_database_check()
#' @importFrom utils menu
#' @noRd
  
scm_database_check=function(template='aseg'){
  
  #package's external data path
  extdatadir=system.file('extdata', package='VertexWiseR')
  
  #Check if cells data are present and download
  if (!dir.exists(paste0(extdatadir,'/scm_database/',template)))
  {
    missingobj=1
    
    prompt = utils::menu(c("Yes", "No"), title=paste0(
      "\nThe database directory was not detected inside VertexWiseR's installed package directory (",paste0(extdatadir,'scm_database/',template,'/'), 
      "). It is needed to smooth, plot and analyse subcortical surfaces from SubCortexMesh. It can be downloaded from the github VertexWiseR directory.\n\nDo you want the scm_database database (~9.6 MB, will be ~19.9 MB unzipped) to be downloaded now?"))
    if (prompt==1) {
      
      #function to check if url exists
      #courtesy of Schwarz, March 11, 2020, CC BY-SA 4.0:
      #https://stackoverflow.com/a/60627969
      valid_url <- function(url_in,t=2){
        con <- url(url_in)
        check <- suppressWarnings(try(open.connection(con,open="rt",timeout=t),silent=TRUE)[1])
        suppressWarnings(try(close.connection(con),silent=TRUE))
        ifelse(is.null(check),TRUE,FALSE)}
      
      #downloader 
      #Check if URL works and avoid returning error but only print message as requested by CRAN:
      url=paste0("https://raw.githubusercontent.com/CogBrainHealthLab/VertexWiseR/main/inst/extdata/scm_database.zip")
      if(valid_url(url)) {
        download.file(url=url,
                      destfile = dest_path)
      } else { 
        stop(paste0("The scm_database.zip failed to be downloaded from the github VertexWiseR directory. Please check your internet connection. Alternatively, you may visit https://github.com/CogBrainHealthLab/VertexWiseR/tree/main/inst/extdata/ and download the object manually.")) 
      }
      
      #unzipper
      dest_path=paste0(extdatadir,'/scm_database.zip')
      out_dir=paste0(extdatadir,'/scm_database')
      #unzip or untar as unzip struggles with long paths
      zip_result <- suppressWarnings(try(unzip(zipfile = dest_path, exdir = out_dir), silent = TRUE))
      if (is.null(zip_result) | inherits(zip_result, "try-error")==TRUE) {untar(tarfile = dest_path, exdir = out_dir)}
      unlink(dest_path) #clear zip file after extraction
      
      if (!dir.exists(out_dir))
      {
        stop(paste0("The scm_database.zip failed to be unzipped via unzip() or untar(). You might want to try manually."))   
      }
      
      #if user refuses, stops if required, just returns a message if optional at this stage
    } else {
      stop("\nSubCortexMesh subcortical data can only be smoothed, plotted and analysed with VertexWiseR's template database.\n") 
    } 
    
  } #silent if database exists
}
  
############################################################################################################################
############################################################################################################################
  
#' @title Aseg subcortical plotting parameters
#'
#' @description For each subcortical region-of-interest (ROI), loads the corresponding points and cells data, and sets up default zoom, size for the plotter if not specified in plot_surf().
#' @param surf_data  A numeric vector (length of V) or a matrix (N rows x V columns), where N is the number of subplots, and V is the number of vertices. It can be the output of ASEGvextract(), or masks outputted by vertex-wise analyses function that match the number of vertices of applicable ROIs.
#' @param size A combined pair of numeric vector indicating the image dimensions (width and height in pixels).
#' @param zoom A numeric value for adjusting the level of zoom on the figures. Default is 1.25 for whole-brain surface and 1.20 for hippocampal surface.
#' @param twocbars A boolean object stating whether the plot will include two bars as per plot_overlay_surf() or not. It matters for the size to encompass them properly for some regions. 
#' @returns Cell data for the left and right hemispheres (only one for the Brain-Stem), and the size and zoom parameters.
#' @noRd

aseg_plot_parameters=function(surf_data,size,zoom,twocbars=FALSE)
{
  #the "aseg_celldir" path's existence should be checked by scm_database_check() already
  aseg_celldir=paste0(system.file(package='VertexWiseR'),'/extdata/scm_database/aseg_points_cells/')
  n_vert=max(dim(t(surf_data)))
             
  #accumbens nuclei
  if (n_vert==2044)
  {
    lh_vert=1022; rh_vert=1022 
    if(is.null(size) & twocbars==FALSE) {size=c(600,400)}
    if(is.null(size) & twocbars==TRUE) {size=c(800,400)}
    if(is.null(zoom)) {zoom=1.4}
  } 
  #amygdalae
  else if (n_vert==3430) 
  {
    lh_vert=1638; rh_vert=1792 
    if(is.null(size) & twocbars==FALSE) {size=c(600,400)}
    if(is.null(size) & twocbars==TRUE) {size=c(800,400)}
    if(is.null(zoom)) {zoom=1.4}
  }
  #caudate
  else if (n_vert==6940)
  {
    lh_vert=3440; rh_vert=3500     
    if(is.null(size) & twocbars==FALSE) {size=c(500,400)}
    if(is.null(size) & twocbars==TRUE) {size=c(500,400)}
    if(is.null(zoom)) {zoom=2.5}
  }
  #cerebella
  else if (n_vert==39214) 
  {
    lh_vert=19559; rh_vert=19664
    if(is.null(size) & twocbars==FALSE) {size=c(400,400)}
    if(is.null(size) & twocbars==TRUE) {size=c(600,400)}
    if(is.null(zoom)) {zoom=1.7}
  }
  #hippocampi
  else if (n_vert==8132) 
  {
    lh_vert=4046; rh_vert=4086 
    if(is.null(size) & twocbars==FALSE) {size=c(400,400)}
    if(is.null(size) & twocbars==TRUE) {size=c(600,400)}
    if(is.null(zoom)) {zoom=1.7}
  }
  #pallidum
  else if (n_vert==3200) 
  {
    lh_vert=1600; rh_vert=1600 
    if(is.null(size) & twocbars==FALSE) {size=c(400,400)}
    if(is.null(size) & twocbars==TRUE) {size=c(700,400)}
    if(is.null(zoom)) {zoom=1.7}
  }
  #putamen
  else if (n_vert==8394) 
  {
    lh_vert=4268; rh_vert=4126 
    if(is.null(size) & twocbars==FALSE) {size=c(400,400)}
    if(is.null(size) & twocbars==TRUE) {size=c(700,400)}
    if(is.null(zoom)) {zoom=1.7}
  }
  #thalamus
  else if (n_vert==7768) 
  {
    lh_vert=3936; rh_vert=3832 
    if(is.null(size) & twocbars==FALSE) {size=c(600,400)}
    if(is.null(size) & twocbars==TRUE) {size=c(600,400)}
    if(is.null(zoom)) {zoom=1.9}
  }
  #ventral diencephala
  else if (n_vert==7144) 
  {
    lh_vert=3550; rh_vert=3594 
    if(is.null(size) & twocbars==FALSE) {size=c(400,400)}
    if(is.null(size) & twocbars==TRUE) {size=c(600,400)}
    if(is.null(zoom)) {zoom=1.9}
  }
  #brain stem
  else if (n_vert==9452) 
  {
    lh_vert=9452; rh_vert=9452 
    if(is.null(size) & twocbars==FALSE) {size=c(400,400)}
    if(is.null(size) & twocbars==TRUE) {size=c(500,400)}
    if(is.null(zoom)) {zoom=1.9}
  }
  #all aseg
  else if (n_vert==95718) 
  {
    lh_vert=95718; rh_vert=95718 
    if(is.null(size)) {size=c(900,400)}
    if(is.null(zoom)) {zoom=1.55}
  }
  else
  {
    stop('the surf_data object given does not match the number of vertices of any subcortical template surfaces.')
  }
  
  #get cell data
  celldat=scm_database_fetcher(n_vert,'points_cells')
  lh_celldat_path=readRDS(celldat$lh_celldat_path)
  rh_celldat_path=readRDS(celldat$rh_celldat_path)
  
  return(list(lh_celldata=lh_celldat_path, lh_vert=lh_vert,
              rh_celldata=rh_celldat_path, rh_vert=rh_vert,
              size=size, zoom=zoom))
}


############################################################################################################################
############################################################################################################################

#' @title SubCortexMesh database fetcher
#'
#' @description For each subcortical region-of-interest (ROI), loads the corresponding edge list
#' @param n_vert A single number or integer, representing the number of vertices contained in the template surface data.
#' @param data_object A string object stating the data object to be extracted from the scm_database ('edgelist', 'ROImap', 'template', or 'points_cells').
#' @param template A string indicating the subcortical segmentation scheme. Default and currently only option is 'aseg'.
#' @details Number of vertices in the (bilateral) matrix for each region-of-interest:
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
#' @returns If data_object=='edgelist': A N x 2 matrix object listing each vertex of the template and the vertices adjacent to it (making an edge together).
#' If data_object=='ROImap': A matrix object with N vertices from a template and each parcellation number the vertices correspond to in 6 atlases (6 columns).
#' If data_object=='template': A string object indicating the path to the ROI's surface template in the scm_database for reticulate to retrieve it.
#' If data_object=='points_cells': A string object indicating the path to the ROI's points and cells data in the scm_database.
#' @noRd

scm_database_fetcher=function(n_vert, data_object, template='aseg'){
  #preset template dimensions
  if (n_vert==2044) {roilabel='accumbens'}
  else if (n_vert==3430) {roilabel='amygdala'}
  else if (n_vert==6940) {roilabel='caudate'}
  else if (n_vert==39214) {roilabel='cerebellum'}
  else if (n_vert==8132) {roilabel='hippocampus'}
  else if (n_vert==3200) {roilabel='pallidum'}
  else if (n_vert==8394) {roilabel='putamen'}
  else if (n_vert==7768) {roilabel='thalamus'}
  else if (n_vert==7144) {roilabel='ventraldc'}
  else if (n_vert==9452) {roilabel='brainstem'}
  else if (n_vert==95718) {roilabel='allaseg'}
    
  #the "database" path's existence should be checked by scm_database_check() already
  if(data_object=='edgelist')
  {
    database=system.file(paste0('extdata/scm_database/',template,'/aseg_edgelists'), package='VertexWiseR')
    load(paste0(database,'/edgelist_',roilabel,'.rdata'))
    edgelist <- get(paste0('edgelist_',roilabel))@data
    return(edgelist)
  }
  
  if(data_object=='ROImap')
  {
    database=system.file(paste0('extdata/scm_database/',template,'/aseg_ROImaps'), package='VertexWiseR')
    load(paste0(database,'/ROImap_',roilabel,'.rdata'))
    ROImap <- get(paste0('ROImap_',roilabel))
    return(ROImap)
  }
  
  if(data_object=='template')
  { 
    database=system.file(paste0('extdata/scm_database/',template,'/aseg_templates'), package='VertexWiseR')
    templatepath=paste0(database, '/', roilabel,'.vtk')
    return(templatepath)
  }
  
  if(data_object=='points_cells')
  {
    database=system.file('extdata/scm_database/',template,'/aseg_points_cells', 
                         package='VertexWiseR')
    if (roilabel=='brainstem') 
    {lh_celldat_path=paste0(database,"/",roilabel,"_points_cells.rds")
     rh_celldat_path=paste0(database,"/",roilabel,"_points_cells.rds")
    } else if (roilabel=='allaseg')
    {lh_celldat_path=paste0(database,"/",roilabel,"_points_cells.rds")
    rh_celldat_path=paste0(database,"/",roilabel,"_points_cells.rds")
    } else {
    lh_celldat_path=paste0(database,"/left-",roilabel,"_points_cells.rds")
    rh_celldat_path=paste0(database,"/right-",roilabel,"_points_cells.rds")
    }
    return(list(lh_celldat_path=lh_celldat_path, 
                rh_celldat_path=rh_celldat_path))
  }
}