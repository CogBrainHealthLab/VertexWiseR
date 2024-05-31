############################################################################################################################
############################################################################################################################
#' @title Vertex-wise analysis
#'
#' @description Fits a linear or linear mixed model with the cortical or hippocampal surface data as the predicted outcome, and returns cluster-thresholded (Random field theory) t-stat map selected contrast.
#'
#' @details The function imports and adapts the \href{https://brainstat.readthedocs.io/en/master/_modules/brainstat/stats/SLM.html#SLM)}{brainstat python library}. 
#' 
#' Output definitions:
#' - `nverts`: number of vertices in the cluster
#' - `P`: p-value of the cluster
#' - `X, Y and Z`: MNI coordinates of the vertex with the highest t-statistic in the cluster.
#' - `tstat`: t statistic of the vertex with the highest t-statistic in the cluster
#' - `region`: the region this highest -statistic vertex is located in, as determined/labelled by the selected atlas 
#'
#' @param model An N X V data.frame object containing N rows for each subject and V columns for each predictor included in the model. This data.frame should not include the random effects variable.
#' @param contrast A numeric vector or object containing the values of the predictor of interest. The cluster-thresholded t-stat maps will be estimated only for this predictor
#' @param random An object containing the values of the random variable (optional)
#' @param surf_data A matrix object containing the surface data, see SURFvextract() or HIPvextract() output format. 
#' @param p A numeric object specifying the p-value to threshold the results (Default is 0.05)
#' @param atlas A numeric integer object corresponding to the atlas of interest. 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360, 5=Destrieux-148.
#' @param smooth_FWHM A numeric vector object specifying the desired smoothing width in mm 
#'
#' @returns A list object containing the cluster level results, thresholded t-stat map, and positive, negative and bidirectional cluster maps.
#' 
#' @examples
#' demodata = readRDS(system.file('demo_data/SPRENG_behdata_site1.rds', 
#' package = 'VertexWiseR'))[1:100,]
#'CTv = readRDS(file = url(paste0("https://github.com",
#'"/CogBrainHealthLab/VertexWiseR/blob/main/inst/demo_data/",
#'"SPRENG_CTv_site1.rds?raw=TRUE")))[1:100,]
#'
#'vertexwise_model=vertex_analysis(model=demodata[,c(2,7)], 
#'contrast=demodata[,7], surf_data = CTv, atlas=1,p = 0.05)
#'print(vertexwise_model$cluster_level_results)
#' @importFrom reticulate import r_to_py
#' @export

##vertex wise analysis with mixed effects
vertex_analysis=function(model,contrast, random, surf_data, p=0.05, atlas=1, smooth_FWHM)  ## atlas: 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360, 5=Destrieux-148; ignored for hippocampal surfaces
{
  
  #Check required python dependencies. If files missing:
  #Will prompt the user to get them in interactive session 
  #Will stop if it's a non-interactive session
  . <- non_interactive <- NULL 
  VWRrequirements(ncol(surf_data))
  if (exists("non_interactive")) { 
    return(cat(non_interactive)) }
  
  #If the contrast/model is a tibble (e.g., taken from a read_csv output)
  #converts the columns to regular data.frame column types
  if ('tbl_df' %in% class(contrast) == TRUE) {
    if (inherits(contrast[[1]],"character")==T) {contrast = contrast[[1]]
    } else {contrast = as.numeric(contrast[[1]])}
  } 
  if ('tbl_df' %in% class(model) == TRUE) {
    model=as.data.frame(model)
    if (NCOL(model)==1) {model = model[[1]]
    } else { for (c in 1:NCOL(model)) { 
      if(inherits(model[,c],"double")==T) {model[,c] = as.numeric(model[,c])}
    }  }
  }
  
  if(inherits(contrast,"integer")==T) {contrast=as.numeric(contrast)}
  
  #recode random variale to numeric
    if(!missing("random"))
    { #recoding subject variable
      random=match(random,unique(random))
    }
  
  ##checks
    #check contrast
    if(NCOL(model)>1)
    {
      for(colno in 1:(NCOL(model)+1))
      {
        if(colno==(NCOL(model)+1))  {warning("contrast is not contained within model")}
        
        if(inherits(contrast,"character")==T) 
        {
          if(identical(contrast,model[,colno]))  {break} 
        } else 
        {
          if(identical(suppressWarnings(as.numeric(contrast)),suppressWarnings(as.numeric(model[,colno]))))  {break}
        }
      }
    }  else
    {
      if(inherits(contrast,"character")==T) 
      {
        if(identical(contrast,model))  {colno=1} 
        else  {warning("contrast is not contained within model")}
      } else
      {
        if(identical(as.numeric(contrast),as.numeric(model)))  {colno=1}
        else  {warning("contrast is not contained within model")}
      }
    }
    
    #check if nrow is consistent for model and surf_data
    if(NROW(surf_data)!=NROW(model))  {stop(paste("The number of rows for surf_data (",NROW(surf_data),") and model (",NROW(model),") are not the same",sep=""))}
    
    #incomplete data check
    idxF=which(complete.cases(model)==F)
    if(length(idxF)>0)
    {
      cat(paste("model contains",length(idxF),"subjects with incomplete data. Subjects with incomplete data will be excluded in the current analysis\n"))
      model=model[-idxF,]
      contrast=contrast[-idxF]
      surf_data=surf_data[-idxF,]
      if(!missing(random)) {random=random[-idxF]}
    }
    
    #check categorical and recode variable
    if(NCOL(model)>1)
    {
      for (column in 1:NCOL(model))
      {
        if(inherits(model[,column],"character")==T) 
        {
          if(length(unique(model[,column]))==2)
          {
            cat(paste("The binary variable '",colnames(model)[column],"' will be recoded with ",unique(model[,column])[1],"=0 and ",unique(model[,column])[2],"=1 for the analysis\n",sep=""))
            
            recode=rep(0,NROW(model))
            recode[model[,column]==unique(model[,column])[2]]=1
            model[,column]=recode
            contrast=model[,colno]
          } else if(length(unique(model[,column]))>2)    {stop(paste("The categorical variable '",colnames(model)[column],"' contains more than 2 levels, please code it into binarized dummy variables",sep=""))}
        }      
      }
    } else
    {
      if(inherits(model,"character")==T) 
      {
        if(length(unique(model))==2)
        {
          cat(paste("The binary variable '",colnames(model),"' will be recoded such that ",unique(model)[1],"=0 and ",unique(model)[2],"=1 for the analysis\n",sep=""))
          
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
    
    #check length of CT data and load the appropriate fsaverage files
    n_vert=ncol(surf_data)
    if(n_vert==20484)
    {
    template="fsaverage5"
     ROImap <- get('ROImap_fs5')
    } else if (n_vert==81924)
    {
      template="fsaverage6"
    ROImap <- get('ROImap_fs6')
    } else if (n_vert==14524)
    {
      if(file.exists(system.file('extdata','hip_template.fs', package='VertexWiseR'))==F)
      {
        cat("\nhip_template.fs is not detected in the current working directory. The hippocampus surface template will be downloaded\n\n")
        
        #function to check if url exists
        #courtesy of Schwarz, March 11, 2020, CC BY-SA 4.0:
        #https://stackoverflow.com/a/60627969
        valid_url <- function(url_in,t=2){
          con <- url(url_in)
          check <- suppressWarnings(try(open.connection(con,open="rt",timeout=t),silent=T)[1])
          suppressWarnings(try(close.connection(con),silent=T))
          ifelse(is.null(check),TRUE,FALSE)}
        
        #Check if URL works and avoid returning error but only print message as requested by CRAN:
        url="https://raw.githubusercontent.com/CogBrainHealthLab/VertexWiseR/main/inst/extdata/hip_template.fs"
        if(valid_url(url)) {
          download.file(url, destfile=paste0(system.file(package='VertexWiseR'),'/extdata/hip_template.fs'),mode = "wb")
        } else { 
          cat("\nhip_template.fs could not be downloaded from the github VertexWiseR directory and the url may be broken. Please check your internet connection or visit https://github.com/CogBrainHealthLab/VertexWiseR/tree/main/inst/extdata to download the object.")
          return() #ends function
        }
      }
        
        brainspace.mesh.mesh_io=reticulate::import("brainspace.mesh.mesh_io", delay_load = TRUE)
        template=brainspace.mesh.mesh_io$read_surface(paste0(system.file(package='VertexWiseR'),'/extdata/hip_template.fs'))
        ROImap <- get('ROImap_HIP')
    } else {stop("data vector should only contain 20484 (fsaverage5), 81924 (fsaverage6) or 14524 (hippocampal vertices) columns")}
  
  ##smoothing
    n_vert=NCOL(surf_data)
    if(missing("smooth_FWHM"))
    {
     cat("smooth_FWHM argument was not given. surf_data will not be smoothed here.\n")
    } else if(smooth_FWHM==0) 
    {
     cat("smooth_FWHM set to 0: surf_data will not be smoothed here.\n")
    } else if(smooth_FWHM>0) 
    {
      cat(paste("surf_data will be smoothed using a ",smooth_FWHM,"mm FWHM kernel", sep=""))
      surf_data=smooth_surf(surf_data, FWHM=smooth_FWHM)
    }
    surf_data[is.na(surf_data)]=0
  
  ##import python libaries
  brainstat.stats=reticulate::import("brainstat.stats", delay_load = TRUE)

  ##fitting model
  #preparing mask for model
  mask=array(rep(T,NCOL(surf_data)))
  maskNA=which(colSums(surf_data != 0) == 0)
  mask[which(colSums(surf_data != 0) == 0)]=F
  
  #fit model
  if(missing("random")) {model0=brainstat.stats$terms$FixedEffect(model, "_check_categorical" = F)}
  else {model0=brainstat.stats$terms$MixedEffect(ran = as.factor(random),fix = model,"_check_categorical" = F)}
  model=brainstat.stats$SLM$SLM(model = model0,
                                contrast=contrast,
                                surf = template,
                                mask=mask,
                                correction=c("fdr", "rft"),
                                cluster_threshold=p)
  model$fit(surf_data)
  
  #extracting tstats
  tstat=model$t
  
  ##extracting positive results
  cluster_pos=reticulate::py_to_r(model$P[["clus"]][[1]]) #pulling out results from brainstat's output
  cluster_pos=cluster_pos[cluster_pos$P<p,] #removing clusters that are not significant
  
  #extracting positive cluster map
  pos_clusterIDmap=model$P$clusid[[1]]
  
  if(NROW(cluster_pos)==0) #if no sig clusters emerged
  {
    cluster_pos="No significant clusters"
    pos_clusterIDmap=rep(0, NCOL(surf_data))
  } else
  {
    #creating new result variables in the cluster_pos objects
    cluster_pos$P=round(cluster_pos$P,3)
    cluster_pos$P[cluster_pos$P==0]="<0.001"
    cluster_pos=cluster_pos[ , !(names(cluster_pos) %in% "resels")] #removing the 'resels' column from the original brainstat output
    cluster_pos$X=NA
    cluster_pos$Y=NA
    cluster_pos$Z=NA
    cluster_pos$tstat=NA
    cluster_pos$region=NA
    
    #entering results for each cluster
    for (clusno in cluster_pos$clusid)
    {
      clus_tstat=tstat
      clus_tstat[pos_clusterIDmap!=clusno]=0
      cluster_pos$tstat[clusno]=round(clus_tstat[which.max(clus_tstat)],2)
      cluster_pos[clusno,4:6]=round(model$coord[,which.max(abs(clus_tstat))],1)
      
      #identifying region by matching the indices
      idx_pos=ROImap[[1]][,atlas][which.max(clus_tstat)]
      if(idx_pos>0){cluster_pos$region[clusno]=ROImap[[2]][,atlas][idx_pos] } ##to deal with desikan atlas missing vertex mappings
      else {cluster_pos$region[clusno]="unknown (use another atlas)"}
      
      remove(clus_tstat,idx_pos)
    }
    #thresholding positive cluster map
    pos_clusterIDmap[pos_clusterIDmap>max(cluster_pos$clusid)]=0
  }
  
  ##extracting negative results
  cluster_neg=reticulate::py_to_r(model$P[["clus"]][[2]]) #pulling out results from brainstat's output
  cluster_neg=cluster_neg[cluster_neg$P<p,] #removing clusters that are not significant
  
  #extracting negative cluster map
  neg_clusterIDmap=model$P$clusid[[2]]
  if(NROW(cluster_neg)==0) #if no sig clusters emerged
  {
    cluster_neg="No significant clusters"
    neg_clusterIDmap=rep(0, NCOL(surf_data))
  } else
  { #creating new result variables in the cluster_pos objects
    cluster_neg$P=round(cluster_neg$P,3)
    cluster_neg$P[cluster_neg$P==0]="<0.001"
    cluster_neg=cluster_neg[ , !(names(cluster_neg) %in% "resels")] #removing the 'resels' column from the original brainstat output
    cluster_neg$X=NA
    cluster_neg$Y=NA
    cluster_neg$Z=NA
    cluster_neg$tstat=NA
    cluster_neg$region=NA
    
    #entering results for each cluster
    for (clusno in cluster_neg$clusid)
    {
      clus_tstat=tstat
      clus_tstat[neg_clusterIDmap!=clusno]=0
      cluster_neg$tstat[clusno]=round(clus_tstat[which.min(clus_tstat)],2)
      cluster_neg[clusno,4:6]=round(model$coord[,which.max(abs(clus_tstat))],1)
      
      #identifying region by matching the indices
      idx_neg=ROImap[[1]][,atlas][which.min(clus_tstat)]
      if(idx_neg>0){cluster_neg$region[clusno]=ROImap[[2]][,atlas][idx_neg] } ##to deal with desikan atlas missing vertex mappings
      else {cluster_neg$region[clusno]="unknown (use another atlas)"}
      
      remove(clus_tstat,idx_neg)
    }
    #thresholding negative cluster map
    neg_clusterIDmap[neg_clusterIDmap>max(cluster_neg$clusid)]=0
  }
  ##combining results from both clusters into a list object
  cluster_results=list(cluster_pos,cluster_neg)
  names(cluster_results)=c("Positive contrast","Negative contrast")
  
  ##combining positive and negative cluster maps
  posc = as.matrix(as.numeric(pos_clusterIDmap))
  negc = as.matrix(as.numeric(neg_clusterIDmap))*-1
  posc[negc!=0,] <- negc[negc!=0,]
  posc[posc==0 & negc==0,] <- NA
  bi_clusterIDmap = posc
  
  ##generating thresholded t-stat vector (for plotting)
  tstat[intersect(which(neg_clusterIDmap==0),which(pos_clusterIDmap==0))]=NA
  tstat[is.na(tstat)]=0
  tstat[maskNA]=NA
  #setting 0s to NA to make vertices with t=0 empty in plots
  tstat[tstat==0]=NA
  
  ##generating positive and negative masks
  posmask=array(rep(0,NCOL(surf_data)))
  posmask[which(tstat>0)]=1
  posmask = t(posmask)
  
  negmask=array(rep(0,NCOL(surf_data)))
  negmask[which(tstat<0)]=1
  negmask = t(negmask)
  
  #listing objects to return
  returnobj=(list(cluster_results,as.numeric(tstat),as.numeric(posmask),as.numeric(negmask),as.numeric(pos_clusterIDmap),as.numeric(neg_clusterIDmap), as.numeric(bi_clusterIDmap)))
  names(returnobj)=c("cluster_level_results","thresholded_tstat_map","pos_mask","neg_mask","pos_clusterIDmap","neg_clusterIDmap", "bi_clusterIDmap")
  return(returnobj)
}
