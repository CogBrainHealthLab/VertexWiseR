## FUNCTION FOR MIXED-EFFECT VERTEX-WISE TFCE ANALYSIS
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

############################################################################################################################
############################################################################################################################
#' @title Vertex-wise analysis with TFCE (mixed effect)
#'
#' @description Fits a linear mixed effects model with the cortical or hippocampal surface data as the predicted outcome, and returns t-stat and TFCE statistical maps for the selected contrast.
#' 
#' @details This TFCE method is adapted from the \href{https://github.com/nilearn/nilearn/blob/main/nilearn/mass_univariate/_utils.py#L7C8-L7C8}{nilearn python library}. 
#' 
#' @param model An N X V data.frame object containing N rows for each subject and V columns for each predictor included in the model.This data.frame should not include the random effects variable.
#' @param contrast A numeric vector or object containing the values of the predictor of interest. The t-stat and TFCE maps will be estimated only for this predictor
#' @param surf_data A matrix object containing the surface data, see SURFvextract() or HIPvextract()  output format. 
#' @param random An object or vector containing the values of the random variable 
#' @param nperm A numeric integer object specifying the number of permutations generated for the subsequent thresholding procedures (default = 100)
#' @param tail A numeric integer object specifying whether to test a one-sided positive (1), one-sided negative (-1) or two-sided (2) hypothesis
#' @param nthread A numeric integer object specifying the number of CPU threads to allocate 
#' @param smooth_FWHM A numeric vector object specifying the desired smoothing width in mm 
#' @param perm_type A string object specifying whether to permute the rows ("row"), between subjects ("between"), within subjects ("within") or between and within subjects ("within_between") for random subject effects. Default is "row". 
#'
#'
#' @returns A list object containing the t-stat and the TFCE statistical maps which can then be subsequently thresholded using TFCE.threshold()
#' @examples
#' demodata = readRDS(system.file('demo_data/SPRENG_behdata_site1.rds', package = 'VertexWiseR'))[1:5,]
#'surf_data = readRDS(file = url(paste0("https://github.com",
#'"/CogBrainHealthLab/VertexWiseR/blob/main/inst/demo_data/",
#'"SPRENG_CTv_site1.rds?raw=TRUE")))[1:5,]
#'
#'TFCE.pos=TFCE.vertex_analysis.mixed(model=demodata[,c(2,7)],
#'contrast=demodata[,7], surf_data,random=demodata[,1], 
#'nperm =5,tail = 1, nthread = 2)
#'
#' #To get significant clusters, you may then run:
#' #results=TFCE.threshold(TFCE.output=TFCE.pos, p=0.05, atlas=1)
#' #results$cluster_level_results
#'
#' @importFrom reticulate import r_to_py
#' @importFrom foreach foreach 
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom doSNOW registerDoSNOW
#' @export


##Main function

TFCE.vertex_analysis.mixed=function(model,contrast, surf_data, random, nperm=100, tail=2, nthread=10, smooth_FWHM, perm_type="row")
{
  #Check required python dependencies. If files missing:
  #Will prompt the user to get them in interactive session 
  #Will stop if it's a non-interactive session 
  check=VWRfirstrun(n_vert=max(dim(t(surf_data))))
  if (!is.null(check)) {return(check)}
  
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
      if(inherits(model[,c],"double")==T | inherits(model[,c],"integer")==T) {model[,c] = as.numeric(model[,c])}
    }  }
  }
  
  if(inherits(contrast,"integer")==T) {contrast=as.numeric(contrast)}
  
  ##checks
  #check random variable and recode to numeric
  if(missing("random"))  {stop("random variable is missing")}
  else 
  { #recoding subject variable
    random=match(random,unique(random))
  }
  
  #check if nrow is consistent for model and surf_data
  if(NROW(surf_data)!=NROW(model))  {stop(paste("The number of rows for surf_data (",NROW(surf_data),") and model (",NROW(model),") are not the same",sep=""))}
  if(length(random)!=NROW(model))  {stop(paste("The number of rows for random (",length(random),") and model (",NROW(model),") are not the same",sep=""))}
  
  #incomplete data check
  idxF=which(complete.cases(model)==F)
  if(length(idxF)>0)
  {
    cat(paste("model contains",length(idxF),"subjects with incomplete data. Subjects with incomplete data will be excluded in the current analysis\n"))
    model=model[-idxF,]
    contrast=contrast[-idxF]
    surf_data=surf_data[-idxF,]
    random=random[-idxF]
  }
  
  #check contrast
  if(NCOL(model)>1)
  {
    for(colno in 1:(NCOL(model)+1))
    {
      if(colno==(NCOL(model)+1))  {warning("contrast is not contained within model")}
      
      if (inherits(contrast,"character")==T) 
      {
        if(identical(contrast,model[,colno]))  {break} 
      } else 
      {
        if(identical(suppressWarnings(as.numeric(contrast)),suppressWarnings(as.numeric(model[,colno]))))  {break}
      }
    }
  }  else
  {
    if (inherits(contrast,"character")==T) 
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
      if(inherits(model[,column],"character")==T)
      {
        if(length(unique(model[,column]))==2)
        {
          cat(paste("The binary variable '",colnames(model)[column],"' will be recoded with ",unique(model[,column])[1],"=0 and ",unique(model[,column])[2],"=1 for the analysis\n",sep=""))
          
          recode=rep(0,NROW(model))
          recode[model[,column]==unique(model[,column])[2]]=1
          model[,column]=recode
          contrast=model[,colno]
        } else if(length(unique(model[,column]))>2)    {stop(paste("The categorical variable '",colnames(model)[column],"' contains more than 2 levels, please code it into binarized dummy variables. 
If it is your random variable and it is non-binarizable, do not include it in the 'model' object.",sep=""))}
      }      
    }
  } else
  {
    if (inherits(model,"character")==T) 
    {
      if(length(unique(model))==2)
      {
        cat(paste("The binary variable '",colnames(model),"' will be recoded such that ",unique(model)[1],"=0 and ",unique(model)[2],"=1 for the analysis\n",sep=""))
        
        recode=rep(0,NROW(model))
        recode[model==unique(model)[2]]=1
        model=recode
        model=model
      } else if(length(unique(model))>2)    {stop(paste("The categorical variable '",colnames(model),"' contains more than 2 levels, please code it into binarized dummy variables",sep=""))}
    }      
  }
  
  #creating local environment
  edgelistenv <- new.env()
  
  #check if surf_data is a multiple-rows matrix and NOT a vector
  if (is.null(nrow(surf_data)) | nrow(surf_data)==1)
  {stop("The surface data must be a matrix containing multiple participants (rows).")}
  
  #check length of CT data and load the appropriate edgelist files
  n_vert=ncol(surf_data)
  if(n_vert==20484)
  {
    edgelist<- get('edgelistfs5') 
    assign("edgelist", edgelist, envir = edgelistenv)
  }
  else if (n_vert==81924)
  {
    edgelist <- get('edgelistfs6') 
    assign("edgelist", edgelist, envir = edgelistenv)
  }
  else if (n_vert==14524)
  {
    edgelist <- get('edgelistHIP')
    assign("edgelist", edgelist, envir = edgelistenv)
  }
  else {stop("The surf_data can only be a matrix with 20484 (fsaverage5), 81924 (fsaverage6) or 14524 (hippocampal vertices) columns.")}
  
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
    cat(paste("surf_data will be smoothed using a ", smooth_FWHM,"mm FWHM kernel\n", sep=""))
    surf_data=smooth_surf(surf_data, FWHM=smooth_FWHM)
  }
  surf_data[is.na(surf_data)]=0
  
  ##unpermuted model
  #preparing mask for model
  inc.vert.idx=which(colSums(surf_data) != 0)
  
  #construct model
  start=Sys.time()
  cat("Estimating unpermuted TFCE image...")
  brainstat.stats.terms=reticulate::import("brainstat.stats.terms", delay_load = TRUE)
  brainstat.stats.SLM=reticulate::import("brainstat.stats.SLM", delay_load = TRUE)
  terms=brainstat.stats.terms$MixedEffect(ran = as.factor(random),fix = model,"_check_categorical" = F)
  
  model.fit=brainstat.stats.SLM$SLM(model = terms,contrast=contrast,correction="None",cluster_threshold=1)
  
  #fit model
  model.fit$fit(surf_data[,inc.vert.idx])
  
  #compute unpermuted TFCE stats
  tmap.orig=rep(0,n_vert)
  tmap.orig[inc.vert.idx]=as.numeric(model.fit$t)  
  
  TFCE.multicore = utils::getFromNamespace("TFCE.multicore", "VertexWiseR")
  TFCE.orig=TFCE.multicore(tmap.orig,tail=tail,nthread=ceiling(nthread/2), envir=edgelistenv) #divide by 2 because more threads made things slower
  
  end=Sys.time()
  
  cat(paste("Completed in",round(difftime(end,start, units="secs"),1),"secs\nEstimating permuted TFCE images...\n",sep=" "))

  ##permuted model
  #generating permutation sequences  
  set.seed(123)
  permseq=matrix(NA, nrow=NROW(model), ncol=nperm)
  
  if(perm_type=="within_between") {for (perm in 1:nperm)  {permseq[,perm]=perm_within_between(random)}} 
  else if(perm_type=="within") {for (perm in 1:nperm)  {permseq[,perm]=perm_within(random)}} 
  else if(perm_type=="between") {for (perm in 1:nperm)  {permseq[,perm]=perm_between(random)}} 
  else if(perm_type=="row") {for (perm in 1:nperm)  {permseq[,perm]=sample.int(NROW(model))}}
  
  #activate parallel processing
  unregister_dopar = function() {
    .foreachGlobals <- utils::getFromNamespace(".foreachGlobals", "foreach"); 
    env =  .foreachGlobals;
    rm(list=ls(name=env), pos=env)
  }
  unregister_dopar()
  
  getClusters = utils::getFromNamespace("getClusters", "VertexWiseR")
  TFCE = utils::getFromNamespace("TFCE", "VertexWiseR")
  
  cl=parallel::makeCluster(nthread)
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl, c("edgelist"), envir=edgelistenv)
  `%dopar%` = foreach::`%dopar%`
  
  #progress bar
  doSNOW::registerDoSNOW(cl)
  pb=txtProgressBar(max = nperm, style = 3)
  progress=function(n) setTxtProgressBar(pb, n)
  opts=list(progress = progress)
  
  #fitting permuted model and extracting max-TFCE values in parallel streams

  start=Sys.time()
  TFCE.max=foreach::foreach(perm=1:nperm, .combine="rbind",.export=c("edgelist","getClusters"),.options.snow = opts)  %dopar%
    {
      #load python libraries within each foreach loop
      brainstat.stats.terms=reticulate::import("brainstat.stats.terms", delay_load = TRUE)
      brainstat.stats.SLM=reticulate::import("brainstat.stats.SLM", delay_load = TRUE)
      
      ##fit permuted models
      terms=brainstat.stats.terms$MixedEffect(ran = random,fix = model,"_check_categorical" = F)
      model.fit=brainstat.stats.SLM$SLM(model = terms,contrast=contrast,correction="None",cluster_threshold=1)
      model.fit$fit(surf_data[permseq[,perm],inc.vert.idx])
      
      tmap.perm=rep(0,n_vert)
      tmap.perm[inc.vert.idx]=as.numeric(model.fit$t)  
        
      return(max(abs(suppressWarnings(TFCE(data = tmap.perm,tail = tail)))))
    }
  end=Sys.time()
  cat(paste("\nCompleted in ",round(difftime(end, start, units='mins'),1)," minutes \n",sep=""))
  unregister_dopar()
  
  
  ##saving list objects
  returnobj=list(tmap.orig,TFCE.orig, TFCE.max,tail)
  names(returnobj)=c("t_stat","TFCE.orig","TFCE.max","tail")
  
  return(returnobj)
}
