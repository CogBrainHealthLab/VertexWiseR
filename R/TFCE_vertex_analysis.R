## FUNCTIONS FOR VERTEX-WISE TFCE ANALYSIS
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

############################################################################################################################
############################################################################################################################
#' @title Vertex-wise analysis with threshold-free cluster enhancement (fixed effect)
#'
#' @description Fits a linear model with the cortical or hippocampal surface data as the predicted outcome, and returns t-stat and threshold-free cluster enhancement (TFCE) statistical maps for the selected contrast.
#' 
#' @details This TFCE method is adapted from the \href{https://github.com/nilearn/nilearn/blob/main/nilearn/mass_univariate/_utils.py#L7C8-L7C8}{ 'Nilearn' Python library}. 
#' 
#' @param model An N X P data.frame object containing N rows for each subject and P columns for each predictor included in the model
#' @param contrast A N x 1 numeric vector or object containing the values of the predictor of interest. Its length should equal the number of subjects in model (and can be a single column from model). The t-stat and TFCE maps will be estimated only for this predictor.
#' @param formula An optional string or formula object describing the predictors to be fitted against the surface data, replacing the model, contrast, or random arguments. If this argument is used, the formula_dataset argument must also be provided.
#' - The dependent variable (DV) is not needed, and the formula will start with ~. The DV will be the surface data value by default, but it can be swapped with contrast as IV via the "inverse" argument.
#' - The first independent variable in the formula will always be interpreted as the contrast of interest for which to estimate cluster-thresholded t-stat maps. 
#' - Only one random regressor can be given and must be indicated as '(1|variable_name)'.
#' @param formula_dataset An optional data.frame object containing the independent variables to be used with the formula (the IV names in the formula must match their column names in the dataset).
#' @param inverse A boolean object stating whether to set the surface data as predictor of the contrast variable, instead of as dependent variable (default is FALSE). Other covariates in the model remain independent variables. This makes modelling slower.
#' @param surf_data A N x V matrix object containing the surface data (N row for each subject, V for each vertex), in fsaverage5 (20484 vertices), fsaverage6 (81924 vertices), fslr32k (64984 vertices) or hippocampal (14524 vertices) space. See also Hipvextract(), SURFvextract() or FSLRvextract output formats. Alternatively, a string object containing the path to the surface object (.rds file) outputted by extraction functions may be given.
#' @param nperm A numeric integer object specifying the number of permutations generated for the subsequent thresholding procedures (default = 100)
#' @param tail A numeric integer object specifying whether to test a one-sided positive (1), one-sided negative (-1) or two-sided (2) hypothesis
#' @param nthread A numeric integer object specifying the number of CPU threads to allocate 
#' @param smooth_FWHM A numeric vector object specifying the desired smoothing width in mm. It should not be specified if the surf_data has been smoothed previously with smooth_surf(), because this result in surf_data being smoothed twice.
#' @param VWR_check A boolean object specifying whether to check and validate system requirements. Default is TRUE.
#'
#' @returns A list object containing the t-stat and the TFCE statistical maps which can then be subsequently thresholded using TFCE_threshold()
#' @seealso \code{\link{RFT_vertex_analysis}}, \code{\link{TFCE_vertex_analysis_mixed}}, \code{\link{TFCE_threshold}}
#'  
#' @examples
#' demodata = readRDS(system.file('demo_data/SPRENG_behdata_site1.rds',
#'package = 'VertexWiseR'))[1:5,]
#'CTv = readRDS(file = url(paste0("https://github.com",
#'"/CogBrainHealthLab/VertexWiseR/blob/main/inst/demo_data/",
#'"SPRENG_CTv_site1.rds?raw=TRUE")))[1:5,]
#'
#' TFCEpos=TFCE_vertex_analysis(model=demodata[,c("sex","age")], 
#' contrast=demodata[,"age"], surf_data=CTv, tail=1, 
#' nperm=5, nthread = 2, VWR_check=FALSE)
#' 
#' #To threshold the results, you may then run:
#' #results=TFCE_threshold(TFCEpos, p=0.05, atlas=1)
#' #results$cluster_level_results
#'
#' #Formula alternative:
#' #formula= as.formula("~ age + sex")
#' #TFCEpos=TFCE_vertex_analysis(formula=formula, 
#' #formula_dataset=demodata, surf_data=CTv, tail=1, 
#' #nperm=5, nthread = 2, VWR_check=FALSE) 
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

TFCE_vertex_analysis=function(model,contrast, formula, formula_dataset, inverse=FALSE, surf_data, nperm=100, tail=2, nthread=10, smooth_FWHM, VWR_check=TRUE)
{
  #gets surface matrix if is surf_data is a list or path
  surf_data=get_surf_obj(surf_data)
  
  #Check required python dependencies. If files missing:
  #Will prompt the user to get them in interactive session 
  #Will stop if it's a non-interactive session
  if (VWR_check == TRUE){
    message("Checking for VertexWiseR system requirements ... ")
    check = VWRfirstrun(requirement="python/conda only")
    if (!is.null(check)) {return(check)} 
  } else if(interactive()==FALSE) { return(message('Non-interactive sessions need requirement checks'))}
  
  #if the user chooses to use a formula, run the formula reader
  #and output appropriate objects
  if (!missing(formula) & !missing(formula_dataset))
  {
    formula_model=model_formula_reader(formula, formula_dataset) 
    model=formula_model$model
    contrast=formula_model$contrast
  } else if ((missing(formula) & !missing(formula_dataset)) | (!missing(formula) & missing(formula_dataset)))
  {stop('The formula and the formula_dataset arguments must both be provided to work.')}
  
  #run all checks for correct structure, recode variables when needed with
  #model_check()
  if (missing(smooth_FWHM)) {smooth_FWHM=NULL};
  model_summary=model_check(model=model, contrast=contrast, 
                            random=NULL, surf_data=surf_data,
                            smooth_FWHM=smooth_FWHM)
  model=model_summary$model;
  contrast=model_summary$contrast;
  surf_data=model_summary$surf_data;
  colno=model_summary$colno
  
  #check length of surface data and load the appropriate edgelist files
  n_vert=ncol(surf_data)
  if(n_vert==20484)  {
    edgelist <- get_edgelist('fsaverage5') 
    assign("edgelist", edgelist)
  }
  else if (n_vert==81924)  {
    edgelist <- get_edgelist('fsaverage6') 
    assign("edgelist", edgelist)
  }
  else if (n_vert==64984)  {
    edgelist <- get_edgelist('fslr32k') 
    assign("edgelist", edgelist)
  }
  else if (n_vert==14524)  {
    edgelist_hip <- get('edgelist_hip')
    edgelist <- edgelist_hip@data
    assign("edgelist", edgelist)
  }
  else {stop("The surf_data can only be a matrix with 20484 (fsaverage5), 81924 (fsaverage6), 64984 (fslr32k) or 14524 (hippocampal vertices) columns.")}
  
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
  
  
  #if model is inverse, prepare covariates without contrast (vertex will be
  #added to it as the contrast further down)
  #identify and drop contrast column if multiple variables
  if (inverse==TRUE)
  {  
    if (!is.null(dim(model)) & dim(as.data.frame(model))[2]!=1)
    {
      for (cont_col in 1:ncol(model)) 
      {if (all(as.vector(model[,cont_col])==contrast)){break}}
      invmodel=model[ , -cont_col, drop = FALSE]
    } else #if only one variable, empty placeholder to add surf_data as IV
    {invmodel=data.frame(matrix(ncol=0, 
                              nrow=nrow(surf_data)))}
  }
  
  ##unpermuted model
  model=data.matrix(model)
  #if inverse, contrast will be a DV, and vertex-wise surface data an IV
  if (inverse==FALSE)
  {
    mod=.lm.fit(y=surf_data,x=data.matrix(cbind(1,model)))
    tmap.orig=extract.t(mod,colno+1)
  } else
  { 
    #inverse models, vertex by vertex
    tmap.orig=c()
    for (vert in 1:ncol(surf_data)) 
    { #one lm for every vertex
      vertmodel=cbind(invmodel,surf_data[,vert])
      invmod=.lm.fit(y=contrast,x=data.matrix(cbind(1,vertmodel)))
      tmap.orig=c(tmap.orig, extract.t(invmod, ncol(vertmodel)+1))
    }
  }
  
  #extract tstat and calculate tfce image
  start=Sys.time()
  message("Estimating unpermuted TFCE image...")
  
  TFCE.orig=suppressWarnings(TFCE.multicore(data = tmap.orig,tail = tail,nthread=nthread, envir=environment(), edgelist=edgelist))
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
  
  doParallel::registerDoParallel(cl)
  #preload variables for the cluster workers
  parallel::clusterExport(cl, c("edgelist","surf_data","permseq","contrast"), envir=environment())
  `%dopar%` = foreach::`%dopar%`
  
  #progress bar
  doSNOW::registerDoSNOW(cl)
  pb=txtProgressBar(max = nperm, style = 3)
  progress=function(n) setTxtProgressBar(pb, n)
  opts=list(progress = progress)
  
  ##fitting permuted regression model and extracting t-stats in parallel streams
  start=Sys.time()
  
  TFCE.max=foreach::foreach(perm=1:nperm, .combine="c", .options.snow = opts)  %dopar%
    {
      ##commented out alternative method of permutation— permuting only the contrast variable
      #model.permuted=model
      #model.permuted[,colno]=model.permuted[permseq[,perm],colno] ##permute only the contrast
      #mod.permuted=lm(surf_data~data.matrix(model.permuted))
      
      #permuted model
      #if inverse, contrast will be a DV, and vertex-wise surface data an IV
      if (inverse==FALSE)
      {
        mod.permuted=.lm.fit(y=surf_data[permseq[,perm],],
                             x=data.matrix(cbind(1,model)))
        tmap=extract.t(mod.permuted,colno+1)
      }
      else
      { 
        #inverse models, vertex by vertex
        tmap=c()
        for (vert in 1:ncol(surf_data)) 
        { #one lm for every vertex
          vertmodel.permuted=cbind(invmodel,surf_data[,vert])
          invmod.permuted=.lm.fit(y=contrast[permseq[,perm]],
                               x=data.matrix(cbind(1,vertmodel.permuted)))
          tmap=c(tmap,extract.t(invmod.permuted,colno+1))
        }
      }
      
      . <- model.permuted <- NULL #visible binding needed if commented out 
      remove(mod.permuted,model.permuted)
      return(max(abs(suppressWarnings(TFCE(data = tmap,
                                           tail = tail,
                                           edgelist=edgelist)))))
    }
  end=Sys.time()
  message(paste("\nCompleted in ",round(difftime(end, start, units='mins'),1)," minutes \n",sep=""))
  parallel::stopCluster(cl)
  unregister_dopar()
  
  ##saving list objects
  returnobj=list(tmap.orig,
                 TFCE.orig, 
                 TFCE.max,
                 tail)
  names(returnobj)=c("t_stat",
                     "TFCE.orig",
                     "TFCE.max",
                     "tail")
  
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
  n_threshs=length(score_threshs)

  #loop across different signs(i.e., for two tailed test)
  for (sign.idx in 1:length(signs)) 
  {
    temp_data = data * signs[sign.idx]
    tfce=matrix(0,nrow=n_threshs, ncol=length(temp_data))
    
    #loop across different score_threshs values for TFCE estimation
    for(thresh.no in 1:n_threshs)
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
          tfce[thresh.no,non_zero_inds] = cluster_tfces[labeled_non_zero]
          remove(non_zero_inds,cluster_tfces,tfce_step_values, labeled_non_zero)
        }
      }
    }
    remove(clust.dat)
    #combine results from positive and negative tails if necessary 
    if(sign.idx==1){tfce_step_values.all=colSums(tfce)}
    else if (sign.idx==2){tfce_step_values.all=tfce_step_values.all+colSums(tfce)}
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
#' @param atlas A numeric integer object corresponding to the atlas of interest.  1=Desikan, 2=Destrieux-148, 3=Glasser-360, 4=Schaefer-100, 5=Schaefer-200, 6=Schaefer-400. Set to `1` by default. This argument is ignored for hippocampal surfaces.
#' @param k Cluster-forming threshold (Default is 20)
#' @param VWR_check A boolean object specifying whether to check and validate system requirements. Default is TRUE.
#'
#' @returns A list object containing the cluster level results, unthresholded t-stat map, thresholded t-stat map, and positive, negative and bidirectional cluster maps.
#' @examples
#' model1_TFCE=readRDS(file = url(paste0("https://github.com/CogBrainHealthLab",
#' "/VertexWiseR/blob/main/inst/demo_data/model1_TFCE.rds?raw=TRUE")))
#' 
#' TFCEanalysis_output=TFCE_threshold(model1_TFCE, p=0.05, atlas=1,
#' VWR_check=FALSE)
#' TFCEanalysis_output$cluster_level_results
#' @export

TFCE_threshold=function(TFCEoutput, p=0.05, atlas=1, k=20, VWR_check = TRUE)
{
  
  #Check required python dependencies. If files missing:
  #Will prompt the user to get them in interactive session 
  #Will stop if it's a non-interactive session
  if (VWR_check == TRUE){
    message("Checking for VertexWiseR system requirements ... ")
    check = VWRfirstrun(requirement="conda/brainstat")
    if (!is.null(check)) {return(check)} 
  } else if(VWR_check == FALSE & interactive()==FALSE) { return(message('Non-interactive sessions need requirement checks'))}
  
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
    edgelist <- get_edgelist('fsaverage5') 
    assign("edgelist", edgelist, envir = internalenv)
    
    ROImap_fs5 <- get('ROImap_fs5')
    ROImap <- list(ROImap_fs5@data,ROImap_fs5@atlases)
    assign("ROImap", ROImap, envir = internalenv)
    
    MNImap <- get_MNIcoords('fsaverage5')
    assign("MNImap", MNImap, envir = internalenv)
  }
  else if (n_vert==64984) 
  {
    edgelist <- get_edgelist('fslr32k') 
    assign("edgelist", edgelist, envir = internalenv)
    
    ROImap_fslr32k <- get('ROImap_fslr32k')
    ROImap <- list(ROImap_fslr32k@data,ROImap_fslr32k@atlases)
    assign("ROImap", ROImap, envir = internalenv)
    
    MNImap <- get_MNIcoords('fslr32k')
    assign("MNImap", MNImap, envir = internalenv)
  } 
  else if (n_vert==81924) 
  {
    edgelist <- get_edgelist('fsaverage6') 
    assign("edgelist", edgelist, envir = internalenv)
    
    ROImap_fs6 <- get('ROImap_fs6')
    ROImap <- list(ROImap_fs6@data,ROImap_fs6@atlases)
    assign("ROImap", ROImap, envir = internalenv)
    
    MNImap <- get_MNIcoords('fsaverage6')
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
    
    brainspace.mesh.mesh_io=reticulate::import("brainspace.mesh.mesh_io",
                                               delay_load = TRUE)
    template=brainspace.mesh.mesh_io$read_surface(paste0(system.file(package='VertexWiseR'),'/extdata/hip_template.fs'))
    MNImap <- t(template$Points)
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
  TFCEoutput$t_stat[TFCEoutput$t_stat==0]=NA
  
  returnobj=list(cluster_level_results, 
                 TFCEoutput$t_stat,
                 t_stat.thresholdedPK,
                 pos.mask,neg.mask, 
                 pos.clustermap, 
                 neg.clustermap, 
                 bi.clusterIDmap)  
  names(returnobj)=c("cluster_level_results",
                     "tstat_map",
                     "thresholded_tstat_map",
                     "pos_mask","neg_mask",
                     "pos_clusterIDmap",
                     "neg_clusterIDmap", 
                     "bi_clusterIDmap")
  returnobj$cluster_level_results
  
  return(returnobj)
}  
