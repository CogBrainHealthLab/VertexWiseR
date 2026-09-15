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
#' @param model An N X P data.frame object containing N rows for each subject and P columns for each predictor included in the model.This data.frame should not include the random effects variable.
#' @param contrast A N x 1 numeric vector or object containing the values of the predictor of interest. Its length should equal the number of subjects in model (and can be a single column from model). The t-stat and TFCE maps will be estimated only for this predictor.
#' @param random A N x 1 numeric vector or object containing the values of the random variable (optional). Its length should be equal to the number of subjects in model (it should NOT be inside the model data.frame).
#' @param formula An optional string or formula object describing the predictors to be fitted against the surface data, replacing the model, contrast, or random arguments. If this argument is used, the formula_dataset argument must also be provided.
#' - The dependent variable (DV) is not needed, and the formula will start with ~. The DV will be the surface data value by default, but it can be swapped with contrast as IV via the "inverse" argument.
#' - The first independent variable in the formula will always be interpreted as the contrast of interest for which to estimate cluster-thresholded t-stat maps. 
#' - Only one random regressor can be given and must be indicated as '(1|variable_name)'.
#' @param formula_dataset An optional data.frame object containing the independent variables to be used with the formula (the IV names in the formula must match their column names in the dataset).
#' @param inverse A boolean object stating whether to set the surface data as predictor of the contrast variable, instead of as dependent variable (default is FALSE). Other covariates in the model remain independent variables. Formula cannot be modified (with surf_data as predictor) to trigger inverse. This makes modelling substantially slower.
#' @param surf_data A N x V matrix object containing the surface data (N row for each subject, V for each vertex), in fsaverage5 (20484 vertices), fsaverage6 (81924 vertices), fslr32k (64984 vertices) or hippocampal (14524 vertices) space. See also Hipvextract(), SURFvextract(), FSLRvextract(), or SCMvextract() output formats. Alternatively, a string object containing the path to the surface object (.rds file) outputted by extraction functions may be given.
#' @param nperm An integer specifying the number of permutations generated for the subsequent thresholding procedures (default = 100)
#' @param tail An integer specifying whether to test a one-sided positive (1), one-sided negative (-1) or two-sided (2) hypothesis
#' @param nthread An integer specifying the number of CPU threads to allocate 
#' @param smooth_FWHM A numeric vector object specifying the desired smoothing width in mm. It should not be specified if the surf_data has been smoothed previously with smooth_surf(), because this results in surf_data being smoothed twice.
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
#' @importFrom mori share
#' @importFrom Rcpp evalCpp
#' @importFrom foreach foreach 
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom doSNOW registerDoSNOW
#' @export


##Main function

TFCE_vertex_analysis_mixed=function(model,contrast, random, formula, formula_dataset, inverse=FALSE, surf_data, nperm=100, tail=2, nthread=10, smooth_FWHM, perm_type="row", VWR_check=TRUE)
{
  #gets surface matrix if is surf_data is a list or path
  surf_data=get_surf_obj(surf_data)
  
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
  model=as.matrix(model_summary$model);
  contrast=as.numeric(model_summary$contrast);
  surf_data=model_summary$surf_data;
  colno = model_summary$colno
  
  if (!is.null(model_summary$random)) {random=model_summary$random}

  #check length of CT data and load the appropriate edgelist files
  n_vert=ncol(surf_data)
  if(n_vert==20484)
  {
    edgelist <- get_edgelist('fsaverage5') 
    assign("edgelist", edgelist)
  }
  else if (n_vert==81924)
  {
    edgelist <- get_edgelist('fsaverage6') 
    assign("edgelist", edgelist)
  }
  else if (n_vert==64984)
  {
    edgelist <- get_edgelist('fslr32k') 
    assign("edgelist", edgelist)
  }
  else if (n_vert==14524)
  {
    edgelist_hip <- get('edgelist_hip')
    edgelist <- edgelist_hip@data
    assign("edgelist", edgelist)
  }
  #subcortexmesh
  else if (n_vert %in% c(2044,3430,6940,39214,8132,3200,8394,7768,7144,9452, 95718, 2026,3592,7570,31466,8244,3548,7908,8542,9516,82412) #fslfirst template
    ){
    #specify template 
    if (n_vert %in% c(2044,3430,6940,39214,8132,3200,8394,7768,7144,9452, 95718))         {template='fsaverage'} 
    else if (n_vert %in% c(2026,3592,7570,31466,8244,3548,7908,8542,9516,82412))         {template='fslfirst'}
    
    scm_database_check(template) #check if database directory is present
    edgelist=scm_database_fetcher(n_vert,'edgelist',template)
    assign("edgelist", edgelist)
  }
  else {stop("The surf_data can only be a matrix with 20484 (fsaverage5), 81924 (fsaverage6), 64984 (fslr32k) or 14524 (hippocampal vertices) columns. For SubCortexMesh subcortices, please refer to ?SCMvextract().")}
  
  ##unpermuted model
  #preparing mask for model
  inc.vert.idx=which(colSums(surf_data) != 0)

  #construct model
  start=Sys.time()
  message("Estimating unpermuted TFCE image...")

  # Numeric group codes also work with shared-memory serialization.
  random = match(random, unique(random))
  
  # Explicit intercept: disable automatic intercept addition below.
  X_full = data.matrix(cbind(1, model))
  X_null = data.matrix(cbind(1, model[, -colno, drop = FALSE]))
  
  ## Observed full model
  if (!inverse) {
    fit_full = lme_fast(
      Y = surf_data[, inc.vert.idx, drop = FALSE],
      X = X_full,
      id = random,
      add_intercept = FALSE
    )
    
    t_included = as.numeric(fit_full$t_stat[colno + 1L, ])
    rm(fit_full)
    
  } else {
    # Behavioral contrast is the response.
    # Each vertex supplies the final fixed-effect predictor.
    t_included = numeric(length(inc.vert.idx))
    
    for (k in seq_along(inc.vert.idx)) {
      vert = inc.vert.idx[k]
      X_inverse = cbind(X_null, surf_data[, vert])
      
      fit_inverse = lme_fast(
        Y = matrix(contrast, ncol = 1L),
        X = X_inverse,
        id = random,
        add_intercept = FALSE
      )
      
      t_included[k] = fit_inverse$t_stat[ncol(X_inverse), 1L]
    }
    
    rm(fit_inverse, X_inverse)
  }
  
  if (length(t_included) != length(inc.vert.idx) ||
      any(!is.finite(t_included))) {
    stop("Invalid t-statistics in the observed model.")
  }
  
  ## Reduced model for residual permutations
  Y_null = if (!inverse) {
    surf_data[, inc.vert.idx, drop = FALSE]
  } else {
    matrix(contrast, ncol = 1L)
  }
  
  fit_null = lme_fast(
    Y = Y_null,
    X = X_null,
    id = random,
    add_intercept = FALSE
  )
  
  # Fixed-effects predictions in the original observation order.
  # Preserve the existing use of marginal residuals.
  fitted_null = X_null %*% fit_null$coefficients
  residuals_null = Y_null - fitted_null
  
  if (any(!is.finite(fitted_null)) ||
      any(!is.finite(residuals_null))) {
    stop("Reduced-model predictions or residuals are nonfinite.")
  }
  
  rm(fit_null, X_full, X_null, Y_null)
  
  # Restore the full vertex layout before computing TFCE.
  tmap.orig = numeric(n_vert)
  tmap.orig[inc.vert.idx] = t_included
  TFCE.orig = TFCE(tmap.orig, tail = tail, edgelist = edgelist)
  
  end=Sys.time()
  
  message(paste("Completed in",round(difftime(end,start, units="secs"),1),"secs\nEstimating permuted TFCE images...\n",sep=" "))

  ##permuted model
  #generating permutation sequences  
  permseq=matrix(NA, nrow=NROW(model), ncol=nperm)
  
  if(perm_type=="within_between") {for (perm in 1:nperm)  {permseq[,perm]=perm_within_between(random)}} 
  else if(perm_type=="within") {for (perm in 1:nperm)  {permseq[,perm]=perm_within(random)}} 
  else if(perm_type=="between") {for (perm in 1:nperm)  {permseq[,perm]=perm_between(random)}} 
  else if(perm_type=="row") {for (perm in 1:nperm)  {permseq[,perm]=sample.int(NROW(model))}}
  
  ##activate parallel processing
    # Convert read-only inputs to shared-memory objects.
    if (!requireNamespace("mori", quietly = TRUE)) {
      stop("Install mori before running the shared-memory analysis.")
    }
    
    fitted_null = share(fitted_null)
    residuals_null = share(residuals_null)
    permseq = share(permseq)
    edgelist = share(edgelist)
    model = share(model)
    random = share(random)
    inc.vert.idx = share(inc.vert.idx)
  
  # Original surface values are only needed as predictors in inverse mode.
  surf_data = if (inverse) share(surf_data) else NULL
  
  unregister_dopar = function() {
    .foreachGlobals <- utils::getFromNamespace(".foreachGlobals", "foreach"); 
    env =  .foreachGlobals;
    #rm(list=ls(name=env), pos=env) #handled by foreach::registerDoSEQ()
  }
  unregister_dopar()
  
  cl=parallel::makeCluster(nthread)
  doParallel::registerDoParallel(cl)
  # Initialize workers BEFORE exporting shared-memory objects.
  parallel::clusterEvalQ(cl, {
    Sys.setenv(
      OPENBLAS_NUM_THREADS = "1",
      OMP_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1"
    )
    
    loadNamespace("mori")
    library(VertexWiseR)
    NULL
  })
  
  parallel::clusterExport(
    cl,
    c(
      "edgelist", "surf_data", "permseq",
      "fitted_null", "residuals_null",
      "inc.vert.idx", "n_vert", "inverse",
      "model", "colno", "random"
    ),
    envir = environment()
  )
  
  `%dopar%` = foreach::`%dopar%`
  
  #progress bar
  doSNOW::registerDoSNOW(cl)
  pb=txtProgressBar(max = nperm, style = 3)
  progress=function(n) setTxtProgressBar(pb, n)
  opts=list(progress = progress)
  
  #fitting permuted model and extracting max-TFCE values in parallel streams
  start=Sys.time()
  TFCE.max = foreach::foreach(
    perm = seq_len(nperm),
    .combine = "c",
    .options.snow = opts,
    .packages = c("VertexWiseR", "mori"),
    .noexport = c("edgelist", "surf_data", "permseq","fitted_null", "residuals_null","inc.vert.idx", "n_vert", "inverse","model", "colno", "random"),
    .errorhandling = "stop"
  ) %dopar% {
    
    idx = permseq[, perm]
    
    # Permute marginal residuals and restore nuisance predictions.
    Y_perm = fitted_null + residuals_null[idx, , drop = FALSE]
    
    if (!inverse) {
      X_full = data.matrix(cbind(1, model))
      
      fit_perm = lme_fast(
        Y = Y_perm,
        X = X_full,
        id = random,
        add_intercept = FALSE
      )
      
      t_included = as.numeric(fit_perm$t_stat[colno + 1L, ])
      
    } else {
      # Permuted behavioral response; fixed surface predictors.
      X_nuisance = data.matrix(
        cbind(1, model[, -colno, drop = FALSE])
      )
      
      t_included = numeric(length(inc.vert.idx))
      
      for (k in seq_along(inc.vert.idx)) {
        vert = inc.vert.idx[k]
        X_inverse = cbind(X_nuisance, surf_data[, vert])
        
        fit_perm = tryCatch(
          lme_fast(
            Y = Y_perm,
            X = X_inverse,
            id = random,
            add_intercept = FALSE
          ),
          error = function(e) {
            stop(
              "Inverse model failed in permutation ", perm,
              " at vertex ", vert, ": ",
              conditionMessage(e),
              call. = FALSE
            )
          }
        )
        
        t_included[k] = fit_perm$t_stat[ncol(X_inverse), 1L]
      }
    }
    
    # Full-mesh edgelist requires the full vertex layout.
    tmap_perm = numeric(n_vert)
    tmap_perm[inc.vert.idx] = t_included
    
    tfce_perm = TFCE(
      data = tmap_perm,
      tail = tail,
      edgelist = edgelist
    )
    
    max(abs(tfce_perm))
  }
  end=Sys.time()
  message(paste("\nCompleted in ",round(difftime(end, start, units='mins'),1)," minutes \n",sep=""))
  parallel::stopCluster(cl)
  unregister_dopar()
  foreach::registerDoSEQ()
  
  ##saving list objects
  returnobj=list(tmap.orig,TFCE.orig, TFCE.max,tail)
  names(returnobj)=c("t_stat","TFCE.orig","TFCE.max","tail")
  
  return(returnobj)
}
