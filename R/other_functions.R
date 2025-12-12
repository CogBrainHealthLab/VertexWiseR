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
 
