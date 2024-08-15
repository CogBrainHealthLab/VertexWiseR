## FUNCTIONS FOR VERTEX-WISE TFCE ANALYSIS
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

############################################################################################################################
############################################################################################################################
#' @title Vertex-wise analysis with TFCE (fixed effect)
#'
#' @description Fits a linear model with the cortical or hippocampal surface data as the predicted outcome, and returns t-stat and TFCE statistical maps for the selected contrast.
#' 
#' @details This TFCE method is adapted from the \href{https://github.com/nilearn/nilearn/blob/main/nilearn/mass_univariate/_utils.py#L7C8-L7C8}{ 'Nilearn' Python library}. 
#' 
#' @param model An N X V data.frame object containing N rows for each subject and V columns for each predictor included in the model
#' @param contrast A N x 1 numeric vector or object containing the values of the predictor of interest. Its length should equal the number of subjects in model (and can be a single column from model). The t-stat and TFCE maps will be estimated only for this predictor.
#' @param surf_data A N x M matrix object containing the surface data (N row for each subject, M for each vertex), in fsaverage5 (20484 vertices), fsaverage6 (81924 vertices) or hippocampal (14524 vertices) space. See also Hipvextract() or SURFvextract() output format. 
#' @param nperm A numeric integer object specifying the number of permutations generated for the subsequent thresholding procedures (default = 100)
#' @param tail A numeric integer object specifying whether to test a one-sided positive (1), one-sided negative (-1) or two-sided (2) hypothesis
#' @param nthread A numeric integer object specifying the number of CPU threads to allocate 
#' @param smooth_FWHM A numeric vector object specifying the desired smoothing width in mm 
#' @param VWR_check A boolean object specifying whether to check and validate system requirements. Default is TRUE.
#'
#' @returns A list object containing the t-stat and the TFCE statistical maps which can then be subsequently thresholded using TFCE_threshold()
#' @seealso \code{\link{RFT_vertex_analysis}}, \code{\link{TFCE_vertex_analysis_mixed}}, \code{\link{TFCE_threshold}}
#'  
#' @examples
#' demodata = readRDS(system.file('demo_data/SPRENG_behdata_site1.rds',
#'package = 'VertexWiseR'))[1:5,]
#'surf_data = readRDS(file = url(paste0("https://github.com",
#'"/CogBrainHealthLab/VertexWiseR/blob/main/inst/demo_data/",
#'"SPRENG_CTv_site1.rds?raw=TRUE")))[1:5,]
#'model=demodata[,c(2,7)]
#'contrast=demodata[,7]
#'
#' TFCEpos=TFCE_vertex_analysis(model, contrast, surf_data, tail=1, 
#' nperm=5, nthread = 2, VWR_check=FALSE)
#' 
#' #To threshold the results, you may then run:
#' #results=TFCE_threshold(TFCEpos, p=0.05, atlas=1)
#' #results$cluster_level_results
#'
#' @importFrom reticulate import r_to_py
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom doSNOW registerDoSNOW
#' @importFrom stats .lm.fit complete.cases cor
#' @importFrom utils download.file install.packages installed.packages setTxtProgressBar txtProgressBar getFromNamespace
#' @export


##Main function

TFCE_vertex_analysis=function(model,contrast, surf_data, nperm=100, tail=2, nthread=10, smooth_FWHM, VWR_check=TRUE)
{
  
  #Check required python dependencies. If files missing:
  #Will prompt the user to get them in interactive session 
  #Will stop if it's a non-interactive session
  if (VWR_check == TRUE){
    message("Checking for VertexWiseR system requirements ... ")
    check = VWRfirstrun(requirement="miniconda only")
    if (!is.null(check)) {return(check)} else {message("\u2713 \n")}
  } else if(interactive()==FALSE) { return(message('Non-interactive sessions need requirement checks'))}
  
  
  #If the contrast/model is a tibble (e.g., taken from a read_csv output)
  #converts the columns to regular data.frame column types
  if ('tbl_df' %in% class(contrast) == TRUE) {
    if (inherits(contrast, "character") == TRUE) {contrast = contrast[[1]]
    } else {contrast = as.numeric(contrast[[1]])}
  } 
  if ('tbl_df' %in% class(model) == TRUE) {
    model=as.data.frame(model)
    if (NCOL(model)==1) {model = model[[1]]
    } else { for (c in 1:NCOL(model)) { 
      if(inherits(model[,c], "double")==TRUE) {model[,c] = as.numeric(model[,c])}
    }  }
  }
  
  if(inherits(contrast,"integer")==TRUE) {contrast=as.numeric(contrast)}
  ##load other vertex-wise functions (not needed when package is in library)
  #source("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/R/otherfunc.R?raw=TRUE")
  
  ##checks
  #check if nrow is consistent for model and surf_data
  if(NROW(surf_data)!=NROW(model))  {stop(paste("The number of rows for surf_data (",NROW(surf_data),") and model (",NROW(model),") are not the same",sep=""))}
  
  #incomplete data check
  idxF=which(complete.cases(model)==FALSE)
  if(length(idxF)>0)
  {
    message(paste("model contains",length(idxF),"subjects with incomplete data. Subjects with incomplete data will be excluded in the current analysis\n"))
    model=model[-idxF,]
    contrast=contrast[-idxF]
    surf_data=surf_data[-idxF,]
  }
  
  #check contrast
  if(NCOL(model)>1)
  {
    for(colno in 1:(NCOL(model)+1))
    {
      if(colno==(NCOL(model)+1))  {warning("contrast is not contained within model")}
      
      if (inherits(contrast, "character")==TRUE) 
      {
        if(identical(contrast,model[,colno]))  {break} 
      } else 
      {
        if(identical(suppressWarnings(as.numeric(contrast)),suppressWarnings(as.numeric(model[,colno]))))  {break}
      }
    }
  }  else
  {
    if (inherits(contrast,"character")==TRUE) 
    {
      if(identical(contrast,model))  {colno=1} 
      else  {warning("contrast is not contained within model")}
    } else
    {
      if(identical(as.numeric(contrast),as.numeric(model)))  {colno=1}
      else  {warning("contrast is not contained within model")}
    }
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
    if (inherits(model,"character")==TRUE) 
    {
      if(length(unique(model))==2)
      {
        message(paste("The binary variable '",colnames(model),"' will be recoded such that ",unique(model)[1],"=0 and ",unique(model)[2],"=1 for the analysis\n",sep=""))
        
        recode=rep(0,NROW(model))
        recode[model==unique(model)[2]]=1
        model=recode
        model=model
      } else if(length(unique(model))>2)    {stop(paste("The categorical variable '",colnames(model),"' contains more than 2 levels, please code it into binarized dummy variables",sep=""))}
    }      
  }
  
  #make internal invironment to save edgelist
  edgelistenv <- new.env()
  
  #check if surf_data is a multiple-rows matrix and NOT a vector
  if (is.null(nrow(surf_data)) | nrow(surf_data)==1)
  {stop("The surface data must be a matrix containing multiple participants (rows).")}
  
  #check length of surface data and load the appropriate edgelist files
  n_vert=ncol(surf_data)
  if(n_vert==20484)  {
    edgelist_fs5 <- get('edgelist_fs5')
    edgelist <- edgelist_fs5@data
    assign("edgelist", edgelist, envir = edgelistenv)
  }
  else if (n_vert==81924)  {
    edgelist_fs6 <- get('edgelist_fs6')
    edgelist <- edgelist_fs6@data
    assign("edgelist", edgelist, envir = edgelistenv)
  }
  else if (n_vert==14524)  {
    edgelist_hip <- get('edgelist_hip')
    edgelist <- edgelist_hip@data
    assign("edgelist", edgelist, envir = edgelistenv)
  }
  else {stop("The surf_data can only be a matrix with 20484 (fsaverage5), 81924 (fsaverage6) or 14524 (hippocampal vertices) columns.")}
  
  #check for collinearity
  if(NCOL(model)>1)
  {
    cormat=cor(model,use = "pairwise.complete.obs")
    cormat.0=cormat
    cormat.0[cormat.0==1]=NA
    if(max(abs(cormat.0),na.rm = TRUE) >0.5)
    {
      warning(paste("correlations among variables in model are observed to be as high as ",round(max(abs(cormat.0),na.rm = TRUE),2),", suggesting potential collinearity among predictors.\nAnalysis will continue...\n",sep=""))
    }
  }
  
  ##smoothing
  n_vert=NCOL(surf_data)
  if(missing("smooth_FWHM"))
  {
    message("smooth_FWHM argument was not given. surf_data will not be smoothed here.\n")
  } else if(smooth_FWHM==0)   
  {
    message("smooth_FWHM set to 0: surf_data will not be smoothed here.\n")
  } else if(smooth_FWHM>0) 
  {
    message(paste("surf_data will be smoothed using a ", smooth_FWHM,"mm FWHM kernel\n", sep=""))
    surf_data=smooth_surf(surf_data, FWHM=smooth_FWHM)
  }
  
  ##unpermuted model
  model=data.matrix(model)
  mod=.lm.fit(y=surf_data,x=data.matrix(cbind(1,model)))
  
  #extract tstat and calculate tfce image
  start=Sys.time()
  message("Estimating unpermuted TFCE image...")
  
  tmap.orig=extract.t(mod,colno+1)
  TFCE.orig=suppressWarnings(TFCE.multicore(data = tmap.orig,tail = tail,nthread=nthread, envir=edgelistenv, edgelist=edgelist))
  remove(mod)
  
  end=Sys.time()
  message(paste("Completed in",round(difftime(end,start, units="secs"),1),"secs\nEstimating permuted TFCE images...\n",sep=" "))
  
  ##permuted models
  #generating permutation sequences  
  permseq=matrix(NA, nrow=NROW(model), ncol=nperm)
  for (perm in 1:nperm)  {permseq[,perm]=sample.int(NROW(model))}
  
  #activate parallel processing
  unregister_dopar = function() {
    .foreachGlobals <- utils::getFromNamespace(".foreachGlobals", "foreach"); env =  .foreachGlobals;
    rm(list=ls(name=env), pos=env)
  }
  unregister_dopar()
  
  cl=parallel::makeCluster(nthread)
  
  #The functions from otherfunc.r below are not exported for users in library(VertexWiseR). For vertTFCE to use them, They are manually exported below beforehand:
  TFCE = utils::getFromNamespace("TFCE", "VertexWiseR")
  extract.t = utils::getFromNamespace("extract.t", "VertexWiseR")
  getClusters = utils::getFromNamespace("getClusters", "VertexWiseR")
  
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl, c("edgelist"), envir=edgelistenv)
  `%dopar%` = foreach::`%dopar%`
  
  #progress bar
  doSNOW::registerDoSNOW(cl)
  pb=txtProgressBar(max = nperm, style = 3)
  progress=function(n) setTxtProgressBar(pb, n)
  opts=list(progress = progress)
  
  ##fitting permuted regression model and extracting t-stats in parallel streams
  start=Sys.time()
  
  TFCE.max=foreach::foreach(perm=1:nperm, .combine="rbind",.export=c("getClusters","edgelist"), .options.snow = opts)  %dopar%
    {
      ##commented out alternative method of permutation— permuting only the contrast variable
      #model.permuted=model
      #model.permuted[,colno]=model.permuted[permseq[,perm],colno] ##permute only the contrast
      #mod.permuted=lm(surf_data~data.matrix(model.permuted))
      
      mod.permuted=.lm.fit(y=surf_data[permseq[,perm],],x=data.matrix(cbind(1,model)))
      tmap=extract.t(mod.permuted,colno+1)
      
      . <- model.permuted <- NULL #visible binding needed if commented out 
      remove(mod.permuted,model.permuted)
      return(max(abs(suppressWarnings(TFCE(data = tmap,tail = tail,
                                           edgelist=edgelist)))))
    }
  end=Sys.time()
  message(paste("\nCompleted in ",round(difftime(end, start, units='mins'),1)," minutes \n",sep=""))
  unregister_dopar()
  
  
  ##saving list objects
  returnobj=list(tmap.orig,TFCE.orig, TFCE.max,tail)
  names(returnobj)=c("t_stat","TFCE.orig","TFCE.max","tail")
  
  return(returnobj)
}
############################################################################################################################
############################################################################################################################

##TFCE single core— for estimating permuted TFCE statistics
##adapted from nilearn python library: https://github.com/nilearn/nilearn/blob/main/nilearn/mass_univariate/_utils.py#L7C8-L7C8
TFCE=function(data,tail=tail,edgelist)
{
  
  #selecting tail type
  if (tail==2) 
  {
    signs = c(-1, 1)
    max_score = max(abs(data),na.rm = TRUE)
  } else if(tail==1)
  {
    signs = 1
    max_score = max(data,na.rm = TRUE)
  } else if(tail==-1)
  {
    signs = -1
    max_score = max(-data,na.rm = TRUE)
  }
  
  #define TFCE parameters
  step=max_score / 100 #calculating number of steps for TFCE estimation
  score_threshs = seq(step, max_score, by = step) #Set based on determined step size
  
  #loop across different signs(i.e., for two tailed test)
  for (sign.idx in 1:length(signs)) 
  {
    temp_data = data * signs[sign.idx]
    tfce=rep(0,length(temp_data))
    
    #loop across different score_threshs values for TFCE estimation
    for(thresh.no in 1:length(score_threshs))
    {
      temp_data[temp_data < score_threshs[thresh.no]] = 0
      if(length(which(temp_data>0))>1) #if less than 2 vertices, skip the following steps
      {
        if (thresh.no>1)
        {
          if(length(clust.dat[[3]])<3)
          {
            clust.dat[[3]]=matrix(clust.dat[[3]],ncol = 2)
          }
          clust.dat=getClusters(temp_data,clust.dat[[3]])
        } else
        {
          clust.dat=getClusters(temp_data,edgelist)  
        }
        
        if (clust.dat[[2]][1]!="noclusters") #if no clusters, skip the following steps
        {
          non_zero_inds = which(clust.dat[[1]] >0)
          labeled_non_zero = clust.dat[[1]][non_zero_inds]
          cluster_tfces = signs[sign.idx] * clust.dat[[2]] * (score_threshs[thresh.no] ^ 2) #using the E=1 , H=2 paramters for 2D (vertex-wise data)
          tfce_step_values = rep(0, length(clust.dat[[1]]))
          tfce_step_values[non_zero_inds] = cluster_tfces[labeled_non_zero]
          tfce[non_zero_inds]=tfce_step_values[non_zero_inds]+tfce[non_zero_inds] #cumulatively add up all TFCE values at each vertex
          remove(non_zero_inds,cluster_tfces,tfce_step_values, labeled_non_zero)
        }
      }
    }
    remove(clust.dat)
    #combine results from positive and negative tails if necessary 
    if(sign.idx==1){tfce_step_values.all=tfce}
    else if (sign.idx==2){tfce_step_values.all=tfce_step_values.all+tfce}
    remove(tfce)
  }
  return(tfce_step_values.all)
}
############################################################################################################################
############################################################################################################################

##TFCE multicore— for estimating unpermuted TFCE statistics
##adapted from nilearn python library: https://github.com/nilearn/nilearn/blob/main/nilearn/mass_univariate/_utils.py#L7C8-L7C8
TFCE.multicore=function(data,tail=tail,nthread,envir,edgelist)
{
  
  #selecting tail type
  if (tail==2) 
  {
    signs = c(-1, 1)
    max_score = max(abs(data),na.rm = TRUE)
  } else if(tail==1)
  {
    signs = 1
    max_score = max(data,na.rm = TRUE)
  } else if(tail==-1)
  {
    signs = -1
    max_score = max(-data,na.rm = TRUE)
  }
  
  #define TFCE parameters
  step=max_score / 100 #calculating number of steps for TFCE estimation
  score_threshs = seq(step, max_score, by = step) #Set based on determined step size
  
  #loop across different signs(i.e., for two tailed test)
  for (sign.idx in 1:length(signs)) 
  {
    temp_data = data * signs[sign.idx]
    tfce=rep(0,length(temp_data))
    
    #activate parallel processing
    unregister_dopar = function() {
      .foreachGlobals <- utils::getFromNamespace(".foreachGlobals", "foreach"); 
      env =  .foreachGlobals;
      rm(list=ls(name=env), pos=env)
    }
    unregister_dopar()
    
    #The functions from otherfunc.r below are not exported for users in library(VertexWiseR). For vertTFCE to use them, They are manually exported below beforehand:
    getClusters = utils::getFromNamespace("getClusters", "VertexWiseR")
    
    cl=parallel::makeCluster(nthread)
    parallel::clusterExport(cl, c("edgelist"), envir=envir)
    doParallel::registerDoParallel(cl)
    `%dopar%` = foreach::`%dopar%`
    
    #Solves the "no visible binding for global variable" issue
    . <- thresh.no <- NULL 
    internalenv <- new.env()
    assign("thresh.no", thresh.no, envir = internalenv)
    
    #parallel loop across different score_threshs values for TFCE estimation
    tfce=foreach::foreach(thresh.no=1:length(score_threshs), .combine="rbind")  %dopar%
      {
        temp_data[temp_data < score_threshs[thresh.no]] = 0
        if(length(which(temp_data>0))>1) #if less than 2 vertices, skip the following steps
        {
          clust.dat=getClusters(temp_data,edgelist)
          if (clust.dat[[2]][1]!="noclusters") #if no clusters, skip the following steps
          {
            non_zero_inds = which(clust.dat[[1]] >0)
            labeled_non_zero = clust.dat[[1]][non_zero_inds]
            cluster_tfces = signs[sign.idx] * clust.dat[[2]] * (score_threshs[thresh.no] ^ 2)
            tfce_step_values = rep(0, length(clust.dat[[1]]))
            tfce_step_values[non_zero_inds] = cluster_tfces[labeled_non_zero]
            return(tfce_step_values)
          }
        }
      }
    #suppressWarnings(closeAllConnections())
    
    #combine results from positive and negative tails if necessary 
    if(length(tfce)>length(temp_data))
    {
      tfce=colSums(tfce)
    } else if(length(tfce)==0)
    {
      tfce=0
    }
    if(sign.idx==1)
    {
      tfce_step_values.all=tfce
    } else if (sign.idx==2)
    {
      tfce_step_values.all=tfce_step_values.all+tfce
    }
  }
  parallel::stopCluster(cl)
  unregister_dopar()
  return(tfce_step_values.all)
}
############################################################################################################################
############################################################################################################################
#' @title Thresholding TFCE output
#'
#' @description Threshold TFCE maps from the TFCE_vertex_analysis() output and identifies significant clusters at the desired threshold. 
#' 
#' @param TFCEoutput An object containing the output from TFCE_vertex_analysis()
#' @param p A numeric object specifying the p-value to threshold the results (Default is 0.05)
#' @param atlas A numeric integer object corresponding to the atlas of interest. 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360, 5=Destrieux-148 (Default is 1)
#' @param k Cluster-forming threshold (Default is 20)
#'
#' @returns A list object containing the cluster level results, thresholded t-stat map, and positive, negative and bidirectional cluster maps.
#' @examples
#' model1_TFCE=readRDS(system.file('demo_data/model1_TFCE.rds', 
#' package = 'VertexWiseR'))
#' 
#' TFCEanalysis_output=TFCE_threshold(model1_TFCE, p=0.05, atlas=1)
#' TFCEanalysis_output$cluster_level_results
#' @export

TFCE_threshold=function(TFCEoutput, p=0.05, atlas=1, k=20)
{
  nperm=length(TFCEoutput$TFCE.max)
  
  #check if number of permutations is adequate
  if(nperm<1/p)  {warning(paste("Not enough permutations were carried out to estimate the p<",p," threshold precisely\nConsider setting an nperm to at least ",ceiling(1/p),sep=""))}
  
  #creating local environment
  internalenv <- new.env()
  
  #function to make sure edgelist is passed to getcluster
  with_env <- function(f, e=internalenv) {
    stopifnot(is.function(f))
    environment(f) <- e
    f
  }
  
  #check which template is used and load appropriate template files
  n_vert=length(TFCEoutput$t_stat)
  if(n_vert==20484) 
  {    
    edgelist_fs5 <- get('edgelist_fs5')
    edgelist <- edgelist_fs5@data
    assign("edgelist", edgelist, envir = internalenv)
    
    ROImap_fs5 <- get('ROImap_fs5')
    ROImap <- list(ROImap_fs5@data,ROImap_fs5@atlases)
    assign("ROImap", ROImap, envir = internalenv)
    
    MNImap_fs5 <- get('MNImap_fs5')
    MNImap <- MNImap_fs5@data
    assign("MNImap", MNImap, envir = internalenv)
  }
  else if (n_vert==81924) 
  {
    edgelist_fs6 <- get('edgelist_fs6')
    edgelist <- edgelist_fs6@data
    assign("edgelist", edgelist, envir = internalenv)
    
    ROImap_fs6 <- get('ROImap_fs6')
    ROImap <- list(ROImap_fs6@data,ROImap_fs6@atlases)
    assign("ROImap", ROImap, envir = internalenv)
    
    MNImap_fs6 <- get('MNImap_fs6')
    MNImap <- MNImap_fs6@data
    assign("MNImap", MNImap, envir = internalenv)
  } 
  else if (n_vert==14524) 
  {
    edgelist_hip <- get('edgelist_hip')
    edgelist <- edgelist_hip@data
    assign("edgelist", edgelist, envir = internalenv)
    
    ROImap_hip <- get('ROImap_hip')
    ROImap <- list(ROImap_hip@data,ROImap_hip@atlases)
    ROImap=list(data.matrix(ROImap[[1]]),ROImap[[2]])
    assign("ROImap", ROImap, envir = internalenv)
    
    MNImap <- get('MNImap_hip')
    MNImap <- MNImap@data
    assign("MNImap", MNImap, envir = internalenv)
  } 
  ##generating p map
  tfce.p=rep(NA,n_vert)
  TFCEoutput$t_stat[is.na(TFCEoutput$t_stat)]=0
  for (vert in 1:n_vert)  {tfce.p[vert]=length(which(TFCEoutput$TFCE.max>abs(TFCEoutput$TFCE.orig[vert])))/nperm}
  
  ##generating thresholded t-stat map
  TFCEoutput$t_stat[is.na(TFCEoutput$t_stat)]=0
  
  t_stat.thresholdedP=TFCEoutput$t_stat
  t_stat.thresholdedP[tfce.p>p]=0
  
  ##Cluster level results
  ##positive cluster
  if(TFCEoutput$tail==1 |TFCEoutput$tail==2)
  {
    #applying p thresholding
    pos.t_stat.thresholdedP=t_stat.thresholdedP
    pos.t_stat.thresholdedP[pos.t_stat.thresholdedP<0]=0
    
    if(length(which(pos.t_stat.thresholdedP!=0))>1) #skip if no clusters detected
    {
      
      pos.clusters0=with_env(getClusters)(pos.t_stat.thresholdedP,edgelist) ## 1st getClusters() to identify all clusters with no. vertices > 1
      #applying k thresholding
      pos.clustID.remove=which(pos.clusters0[[2]]<k)
      pos.clusters0[[1]][which(!is.na(match(pos.clusters0[[1]],pos.clustID.remove)))]=NA
      
      #generating mask
      pos.clusters=with_env(getClusters)(pos.clusters0[[1]],edgelist) ## 2nd getClusters() to identify all clusters from the k-thresholded clustermap
      pos.clusters[[1]][is.na(pos.clusters[[1]])]=0
      pos.mask=rep(0,n_vert)
      ROImap[[1]][, atlas]
      ROImap=list(data.matrix(ROImap[[1]]),ROImap[[2]])
      
      #results table
      pos.clustermap=rep(NA,n_vert)
      
      if(pos.clusters[[2]][1]!="noclusters")
      {
        pos.mask[pos.clusters[[1]]>0]=1
        pos.clust.results=data.frame(matrix(NA,nrow=length(pos.clusters[[2]]), ncol=8))
        colnames(pos.clust.results)=c("clusid","nverts","P","X","Y","Z","tstat","region")
        clust.idx=1
        
        for(clust.no in order(pos.clusters[[2]],decreasing = TRUE))
        {
          clust.vert.idx=which(pos.clusters[[1]]==clust.no)
          pos.clustermap[clust.vert.idx]=clust.idx
          
          pos.clust.results[clust.idx,1]=clust.idx
          pos.clust.results[clust.idx,2]=length(clust.vert.idx)
          max.vert.idx=clust.vert.idx[which(abs(TFCEoutput$t_stat[clust.vert.idx])==max(abs(TFCEoutput$t_stat[clust.vert.idx]),na.rm = TRUE))[1]]
          pos.clust.results[clust.idx,3]=round(tfce.p[max.vert.idx],3)
          if(pos.clust.results[clust.idx,3]==0) {pos.clust.results[clust.idx,3]=paste("<",1/nperm,sep="")}
          pos.clust.results[clust.idx,c(4,5,6)]=round(MNImap[,max.vert.idx],1)
          pos.clust.results[clust.idx,7]=round(abs(TFCEoutput$t_stat[max.vert.idx]),2)
          
          atlas.idx=ROImap[[1]][,atlas][max.vert.idx]
          if(atlas.idx>0) {pos.clust.results[clust.idx,8]=ROImap[[2]][,atlas][atlas.idx] } ##to deal with desikan atlas missing vertex mappings
          else {pos.clust.results[clust.idx,8]="unknown (use another atlas)"}
          
          clust.idx=clust.idx+1
        }
      } else if(pos.clusters[[2]][1]=="noclusters")
      { 
        pos.clust.results="No significant clusters"
        pos.clustermap="No significant clusters"
      }
    } else ## 2nd getClusters()
    {
      pos.clust.results="No significant clusters"
      pos.clustermap="No significant clusters"
      pos.mask=rep(0,n_vert)
    }
  } else if(TFCEoutput$tail==-1)
  {
    pos.clust.results="Positive contrast not analyzed, only negative one-tailed TFCE statistics were estimated"
    pos.clustermap="No significant clusters"
    pos.mask=rep(0,n_vert)
  } 
  
  ##negative cluster
  if(TFCEoutput$tail==-1 |TFCEoutput$tail==2)
  {
    #applying p thresholding
    neg.t_stat.thresholdedP=t_stat.thresholdedP
    neg.t_stat.thresholdedP[neg.t_stat.thresholdedP>0]=0
    
    if(length(which(neg.t_stat.thresholdedP!=0))>1) #skip if no clusters detected
    {
      neg.clusters0=with_env(getClusters)(neg.t_stat.thresholdedP,edgelist) ## 1st getCluster() to identify all clusters with no. vertices > 1
      #applying k thresholding
      neg.clustID.remove=which(neg.clusters0[[2]]<k)
      neg.clusters0[[1]][which(!is.na(match(neg.clusters0[[1]],neg.clustID.remove)))]=NA
      
      #generating mask
      neg.clusters=with_env(getClusters)(neg.clusters0[[1]],edgelist) ## 2nd getCluster() to identify all clusters from the k-thresholded clustermap
      neg.clusters[[1]][is.na(neg.clusters[[1]])]=0
      neg.mask=rep(0,n_vert)
      
      #results table
      neg.clustermap=rep(NA,n_vert)
      
      if(neg.clusters[[2]][1]!="noclusters")
      {
        neg.mask[neg.clusters[[1]]>0]=1
        neg.clust.results=data.frame(matrix(NA,nrow=length(neg.clusters[[2]]), ncol=8))
        colnames(neg.clust.results)=c("clusid","nverts","P","X","Y","Z","tstat","region")
        clust.idx=1
        
        for(clust.no in order(neg.clusters[[2]],decreasing = TRUE))
        {
          clust.vert.idx=which(neg.clusters[[1]]==clust.no)
          neg.clustermap[clust.vert.idx]=clust.idx
          
          neg.clust.results[clust.idx,1]=clust.idx
          neg.clust.results[clust.idx,2]=length(clust.vert.idx)
          max.vert.idx=clust.vert.idx[which(abs(TFCEoutput$t_stat[clust.vert.idx])==max(abs(TFCEoutput$t_stat[clust.vert.idx]),na.rm = TRUE))[1]]
          neg.clust.results[clust.idx,3]=round(tfce.p[max.vert.idx],3)
          if(neg.clust.results[clust.idx,3]==0) {neg.clust.results[clust.idx,3]=paste("<",1/nperm,sep="")}
          neg.clust.results[clust.idx,c(4,5,6)]=round(MNImap[,max.vert.idx],1)
          neg.clust.results[clust.idx,7]=round(abs(TFCEoutput$t_stat[max.vert.idx]),2)
          
          atlas.idx=ROImap[[1]][,atlas][max.vert.idx]
          if(atlas.idx>0) {neg.clust.results[clust.idx,8]=ROImap[[2]][,atlas][atlas.idx] } ##to deal with desikan atlas missing vertex mappings
          else {neg.clust.results[clust.idx,8]="unknown (use another atlas)"}
          
          clust.idx=clust.idx+1
        }
      } else if(neg.clusters[[2]][1]=="noclusters")
      { 
        neg.clust.results="No significant clusters"
        neg.clustermap="No significant clusters"
      } 
    } else  ## 2nd getClusters()
    {
      neg.clust.results="No significant clusters"
      neg.clustermap="No significant clusters"
      neg.mask=rep(0,n_vert)
    }
  }
  else if(TFCEoutput$tail==1)
  {
    neg.clust.results="Negative contrast not analyzed, only negative one-tailed TFCE statistics were estimated"
    neg.clustermap="No significant clusters"
    neg.mask=rep(0,n_vert)
  } 
  
  ##combining positive and negative cluster maps
  #when significant clusters exist in both directions
  if (inherits(pos.clustermap,"character")!=TRUE & inherits(neg.clustermap,"character")!=TRUE) {
    posc = as.matrix(as.numeric(pos.clustermap))
    negc = as.matrix(as.numeric(neg.clustermap))*-1
    posc[!is.na(negc!=0),] <- negc[!is.na(negc!=0),]
    posc[posc==0 & negc==0,] <- NA
    bi.clusterIDmap = posc
  } else if (inherits(pos.clustermap,"character")!=TRUE) {
    bi.clusterIDmap=pos.clustermap
  } else if (inherits(neg.clustermap,"character")!=TRUE) {
    bi.clusterIDmap=neg.clustermap 
  } else if (inherits(pos.clustermap,"character")==TRUE & inherits(neg.clustermap,"character")==TRUE) {
    bi.clusterIDmap="No significant clusters"
  }
  
  ##saving list objects
  cluster_level_results=list(pos.clust.results,neg.clust.results)
  names(cluster_level_results)=c("Positive contrast", "Negative contrasts")
  
  t_stat.thresholdedPK=TFCEoutput$t_stat*(pos.mask+neg.mask)
  #setting 0s to NA to make vertex with t=0 empty in plots
  t_stat.thresholdedPK[t_stat.thresholdedPK==0]=NA
  
  returnobj=list(cluster_level_results, t_stat.thresholdedPK,pos.mask,neg.mask, pos.clustermap, neg.clustermap, bi.clusterIDmap)  
  names(returnobj)=c("cluster_level_results","thresholded_tstat_map","pos_mask","neg_mask","pos_clusterIDmap","neg_clusterIDmap", "bi_clusterIDmap")
  returnobj$cluster_level_results
  
  return(returnobj)
}  