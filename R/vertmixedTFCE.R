## FUNCTION FOR MIXED-EFFECT VERTEX-WISE TFCE ANALYSIS
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

############################################################################################################################
############################################################################################################################
#' @title Vertex-wise analysis with threshold-free cluster enhancement (mixed effect)
#'
#' @description Fits a linear mixed effects model with the cortical or hippocampal surface data as the predicted outcome, and returns t-stat and threshold-free cluster enhancement (TFCE) statistical maps for the selected contrast.
#' 
#' @details This TFCE method is adapted from the \href{https://github.com/nilearn/nilearn/blob/main/nilearn/mass_univariate/_utils.py#L7C8-L7C8}{'Nilearn' Python library}. 
#' 
#' @param model An N X V data.frame object containing N rows for each subject and V columns for each predictor included in the model.This data.frame should not include the random effects variable.
#' @param contrast A N x 1 numeric vector or object containing the values of the predictor of interest. Its length should equal the number of subjects in model (and can be a single column from model). The t-stat and TFCE maps will be estimated only for this predictor.
#' @param random A N x 1 numeric vector or object containing the values of the random variable (optional). Its length should be equal to the number of subjects in model (it should NOT be inside the model data.frame).
#' @param formula An optional string or formula object describing the predictors to be fitted against the surface data, replacing the model, contrast, or random arguments. If this argument is used, the formula_dataset argument must also be provided.
#' - The dependent variable is not needed, as it will always be the surface data values. 
#' - The first independent variable in the formula will always be interpreted as the contrast of interest for which to estimate cluster-thresholded t-stat maps. 
#' - Only one random regressor can be given and must be indicated as '(1|variable_name)'.
#' @param formula_dataset An optional data.frame object containing the independent variables to be used with the formula (the IV names in the formula must match their column names in the dataset).
#' @param surf_data A N x M matrix object containing the surface data (N row for each subject, M for each vertex), in fsaverage5 (20484 vertices), fsaverage6 (81924 vertices), fslr32k (64984 vertices) or hippocampal (14524 vertices) space. See also Hipvextract(), SURFvextract() or FSLRvextract output formats.
#' @param nperm A numeric integer object specifying the number of permutations generated for the subsequent thresholding procedures (default = 100)
#' @param tail A numeric integer object specifying whether to test a one-sided positive (1), one-sided negative (-1) or two-sided (2) hypothesis
#' @param nthread A numeric integer object specifying the number of CPU threads to allocate 
#' @param smooth_FWHM A numeric vector object specifying the desired smoothing width in mm 
#' @param perm_type A string object specifying whether to permute the rows ("row"), between subjects ("between"), within subjects ("within") or between and within subjects ("within_between") for random subject effects. Default is "row". 
#' @param VWR_check A boolean object specifying whether to check and validate system requirements. Default is TRUE.
#'
#' @returns A list object containing the t-stat and the TFCE statistical maps which can then be subsequently thresholded using TFCE_threshold()
#' 
#' @seealso \code{\link{RFT_vertex_analysis}},  \code{\link{TFCE_vertex_analysis}}, \code{\link{TFCE_threshold}}
#' 
#' @examples
#' demodata = readRDS(system.file('demo_data/SPRENG_behdata_site1.rds', package = 'VertexWiseR'))[1:5,]
#'CTv = readRDS(file = url(paste0("https://github.com",
#'"/CogBrainHealthLab/VertexWiseR/blob/main/inst/demo_data/",
#'"SPRENG_CTv_site1.rds?raw=TRUE")))[1:5,]
#'
#'TFCEpos=TFCE_vertex_analysis_mixed(model=demodata[,c("sex",
#'"age")], contrast=demodata[,"age"], random=demodata[,"id"], 
#'surf_data=CTv, nperm =5,tail = 1, nthread = 2, VWR_check=FALSE)
#'
#' #To get significant clusters, you may then run:
#' #results=TFCE_threshold(TFCEpos, p=0.05, atlas=1)
#' #results$cluster_level_results
#'
#' #Formula alternative:
#' #formula= as.formula("~ age + sex + (1|id)")
#' #TFCEpos=TFCE_vertex_analysis_mixed(formula=formula, 
#' #formula_dataset=demodata, surf_data=CTv, tail=1, 
#' #nperm=5, nthread = 2, VWR_check=FALSE) 
#'
#' @importFrom reticulate import r_to_py
#' @importFrom foreach foreach 
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom doSNOW registerDoSNOW
#' @export


##Main function

TFCE_vertex_analysis_mixed=function(model,contrast, random, formula, formula_dataset, surf_data, nperm=100, tail=2, nthread=10, smooth_FWHM, perm_type="row", VWR_check=TRUE)
{
  #Check required python dependencies. If files missing:
  #Will prompt the user to get them in interactive session 
  #Will stop if it's a non-interactive session
  if (VWR_check == TRUE){
    message("Checking for VertexWiseR system requirements ... ")
    check = VWRfirstrun(n_vert=max(dim(t(surf_data))))
    if (!is.null(check)) {return(check)} 
  } else if(interactive()==FALSE) { return(message('Non-interactive sessions need requirement checks'))}
  
  #if the user chooses to use a formula, run the formula reader
  #and output appropriate objects
  if (!missing(formula) & !missing(formula_dataset))
  {
    formula_model=model_formula_reader(formula, formula_dataset) 
    model=formula_model$model
    contrast=formula_model$contrast
    if (!is.null(formula_model$random)) 
    {random=formula_model$random}
  } else if ((missing(formula) & !missing(formula_dataset)) | (!missing(formula) & missing(formula_dataset)))
  {stop('The formula and the formula_dataset arguments must both be provided to work.')}
  
  #run all checks for correct structure, recode variables when needed with
  #model_check()
  if (missing(random)) {random=NULL};
  if (missing(smooth_FWHM)) {smooth_FWHM=NULL};
  model_summary=model_check(model=model, contrast=contrast, 
                            random=random, surf_data=surf_data,
                            smooth_FWHM=smooth_FWHM)
  model=model_summary$model;
  contrast=model_summary$contrast;
  surf_data=model_summary$surf_data;
  if (!is.null(model_summary$random)) {random=model_summary$random}

  #creating local environment
  edgelistenv <- new.env()
  
  #check length of CT data and load the appropriate edgelist files
  n_vert=ncol(surf_data)
  if(n_vert==20484)
  {
    edgelist <- get_edgelist('fsaverage5') 
    assign("edgelist", edgelist, envir = edgelistenv)
  }
  else if (n_vert==81924)
  {
    edgelist <- get_edgelist('fsaverage6') 
    assign("edgelist", edgelist, envir = edgelistenv)
  }
  else if (n_vert==64984)
  {
    edgelist <- get_edgelist('fslr32k') 
    assign("edgelist", edgelist, envir = edgelistenv)
  }
  else if (n_vert==14524)
  {
    edgelist_hip <- get('edgelist_hip')
    edgelist <- edgelist_hip@data
    assign("edgelist", edgelist, envir = edgelistenv)
  }
  else {stop("The surf_data can only be a matrix with 20484 (fsaverage5), 81924 (fsaverage6), 64984 (fslr32k) or 14524 (hippocampal vertices) columns.")}
  
  ##unpermuted model
  #preparing mask for model
  inc.vert.idx=which(colSums(surf_data) != 0)

  #construct model
  start=Sys.time()
  message("Estimating unpermuted TFCE image...")
  brainstat.stats.terms=reticulate::import("brainstat.stats.terms", delay_load = TRUE)
  terms=brainstat.stats.terms$MixedEffect(ran = as.factor(random),fix = model,"_check_categorical" = FALSE)
  
  #Solves the "no visible binding for global variable" issue
  . <- SLM <- NULL 
  #read version of SLM that allows to specify the directory for the
  #fetch_template_surface option
  reticulate::source_python(paste0(system.file(package='VertexWiseR'),'/python/brainstat.stats.SLM_VWR.py'))
  
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
  
  model.fit=SLM(model = terms,
                contrast=contrast,
                correction="None",
                cluster_threshold=1, 
                data_dir=data_dir)
  
  #fit will fetch parcellation data in a different place
  model.fit$data_dir=paste0(brainstat_data_path,'/brainstat_data/parcellation_data/')
  
  #fit model
  SLM$fit(model.fit,surf_data[,inc.vert.idx])
  
  #compute unpermuted TFCE stats
  tmap.orig=rep(0,n_vert)
  tmap.orig[inc.vert.idx]=as.numeric(model.fit$t)  
  
  TFCE.multicore = utils::getFromNamespace("TFCE.multicore", "VertexWiseR")
  TFCE.orig=TFCE.multicore(tmap.orig,tail=tail,nthread=ceiling(nthread/2), envir=edgelistenv, edgelist=edgelist) #divide by 2 because more threads made things slower
  
  end=Sys.time()
  
  message(paste("Completed in",round(difftime(end,start, units="secs"),1),"secs\nEstimating permuted TFCE images...\n",sep=" "))

  ##permuted model
  #generating permutation sequences  
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
  
  #preset path to SLM
  SLMpath= paste0(system.file(package='VertexWiseR'),'/python/brainstat.stats.SLM_VWR.py')
  
  #fitting permuted model and extracting max-TFCE values in parallel streams
  start=Sys.time()
  TFCE.max=foreach::foreach(perm=1:nperm, .combine="rbind",.export="getClusters",.options.snow = opts)  %dopar%
    {
      #load python libraries within each foreach loop
      brainstat.stats.terms=reticulate::import("brainstat.stats.terms", delay_load = TRUE)
      #read version of SLM that allows to specify the directory for the
      #fetch_template_surface option
      reticulate::source_python(SLMpath)
      
      #Brainstat data, will either be stored in default $HOME path or 
      #custom if it's been set via VWRfirstrun()
      if (Sys.getenv('BRAINSTAT_DATA')=="")
      {brainstat_data_path=fs::path_home()} else if 
      (!Sys.getenv('BRAINSTAT_DATA')=="") 
      {brainstat_data_path=Sys.getenv('BRAINSTAT_DATA')}
      #convert path to pathlib object for brainstat
      data_dir=paste0(brainstat_data_path,'/brainstat_data/surface_data/')
      
      ##fit permuted models
      terms=brainstat.stats.terms$MixedEffect(ran = random,fix = model,"_check_categorical" = FALSE)
      model.fit=SLM(model = terms,
                    contrast=contrast,
                    correction="None",
                    cluster_threshold=1, 
                    data_dir=data_dir)
      
      #fit will fetch parcellation data in a different place
      model.fit$data_dir=paste0(brainstat_data_path,'/brainstat_data/parcellation_data/')
      
      #fit model
      SLM$fit(model.fit,surf_data[permseq[,perm],inc.vert.idx])
      
      tmap.perm=rep(0,n_vert)
      tmap.perm[inc.vert.idx]=as.numeric(model.fit$t)  
        
      return(max(abs(suppressWarnings(TFCE(data = tmap.perm,tail = tail,
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
