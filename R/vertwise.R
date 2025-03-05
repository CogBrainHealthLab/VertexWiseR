############################################################################################################################
############################################################################################################################
#' @title Vertex-wise analysis with random field theory cluster correction
#'
#' @description Fits a linear or linear mixed model with the cortical or hippocampal surface data as the predicted outcome, and returns cluster-thresholded (Random field theory) t-stat map selected contrast.
#'
#' @details The function imports and adapts the \href{https://brainstat.readthedocs.io/en/master/_modules/brainstat/stats/SLM.html#SLM)}{ 'BrainStat' Python library}. 
#' 
#' By default, false discovery rate correction is used together with the Random field theory (RFT) cluster correction. To look at data without any form of cluster correction, users can simply refer to the outputted 'tstat_map'.
#' 
#' Output definitions:
#' - `nverts`: number of vertices in the cluster
#' - `P`: p-value of the cluster
#' - `X, Y and Z`: MNI coordinates of the vertex with the highest t-statistic in the cluster.
#' - `tstat`: t statistic of the vertex with the highest t-statistic in the cluster
#' - `region`: the region this highest -statistic vertex is located in, as determined/labelled by the selected atlas 
#'
#' @param model An N X P data.frame object containing N rows for each subject and P columns for each predictor included in the model. This data.frame should not include the random effects variable.
#' @param contrast A N x 1 numeric vector or object containing the values of the predictor of interest. Its length should equal the number of subjects in model (and can be a single column from model). The cluster-thresholded t-stat maps will be estimated only for this predictor. 
#' @param random A N x 1 numeric vector or object containing the values of the random variable (optional). Its length should be equal to the number of subjects in model (it should NOT be inside the model data.frame).
#' @param formula An optional string or formula object describing the predictors to be fitted against the surface data, replacing the model, contrast, or random arguments. If this argument is used, the formula_dataset argument must also be provided.
#' - The dependent variable is not needed, as it will always be the surface data values. 
#' - The first independent variable in the formula will always be interpreted as the contrast of interest for which to estimate cluster-thresholded t-stat maps. 
#' - Only one random regressor can be given and must be indicated as '(1|variable_name)'.
#' @param formula_dataset An optional data.frame object containing the independent variables to be used with the formula (the IV names in the formula must match their column names in the dataset).
#' @param surf_data A N x V matrix object containing the surface data (N row for each subject, V for each vertex), in fsaverage5 (20484 vertices), fsaverage6 (81924 vertices), fslr32k (64984 vertices) or hippocampal (14524 vertices) space. See also Hipvextract(), SURFvextract() or FSLRvextract output formats. Alternatively, a string object containing the path to the surface object (.rds file) outputted by extraction functions may be given.
#' @param p A numeric object specifying the p-value to threshold the results (Default is 0.05)
#' @param atlas A numeric integer object corresponding to the atlas of interest.  1=Desikan, 2=Destrieux-148, 3=Glasser-360, 4=Schaefer-100, 5=Schaefer-200, 6=Schaefer-400. Set to `1` by default. This argument is ignored for hippocampal surfaces.
#' @param smooth_FWHM A numeric vector object specifying the desired smoothing width in mm. It should not be specified if the surf_data has been smoothed previously with smooth_surf(), because this result in surf_data being smoothed twice.
#' @param VWR_check A boolean object specifying whether to check and validate system requirements. Default is TRUE.
#'
#' @returns A list object containing the cluster level results, unthresholded t-stat map, thresholded t-stat map, positive, negative and bidirectional cluster maps, and a FDR-corrected p-value map.
#' 
#' @seealso \code{\link{TFCE_vertex_analysis}}, \code{\link{TFCE_vertex_analysis_mixed}}
#' 
#' @examples
#' demodata = readRDS(system.file('demo_data/SPRENG_behdata_site1.rds', 
#' package = 'VertexWiseR'))[1:100,]
#' CTv = readRDS(file = url(paste0("https://github.com",
#'"/CogBrainHealthLab/VertexWiseR/blob/main/inst/demo_data/",
#'"SPRENG_CTv_site1.rds?raw=TRUE")))[1:100,] 
#'
#' vertexwise_model=RFT_vertex_analysis(model=demodata[,c("sex","age")], 
#' contrast=demodata[,"age"], surf_data = CTv, atlas=1,p = 0.05, 
#' VWR_check=FALSE)
#' 
#' #Description of the output:
#' #vertexwise_model$cluster_level_results
#' 
#' #Formula alternative:
#' #formula= as.formula("~ age + sex")
#' #vertexwise_model=RFT_vertex_analysis(formula=formula, 
#' #formula_dataset=demodata, surf_data = CTv, atlas=1, p = 0.05, 
#' #VWR_check=FALSE)
#' @importFrom reticulate import r_to_py
#' @importFrom stats as.formula na.pass 
#' @export

##vertex wise analysis with mixed effects
RFT_vertex_analysis=function(model,contrast, random, formula, formula_dataset, surf_data, p=0.05, atlas=1, smooth_FWHM, VWR_check=TRUE)  ## atlas: 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360, 5=Destrieux-148; ignored for hippocampal surfaces
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
  
  #run all checks for correct structure, recode variables when needed with model_check()
  if (missing(random)) {random=NULL};
  if (missing(smooth_FWHM)) {smooth_FWHM=NULL};
  model_summary=model_check(model=model, contrast=contrast, 
                            random=random, surf_data=surf_data,
                            smooth_FWHM=smooth_FWHM)
  model=model_summary$model
  contrast=model_summary$contrast
  surf_data=model_summary$surf_data
  if (!is.null(model_summary$random)) {random=model_summary$random}
 
    #check length of CT data and load the appropriate fsaverage files
    n_vert=ncol(surf_data)
    if(n_vert==20484)
    {
      template="fsaverage5"
      ROImap_fs5 <- get('ROImap_fs5')
      ROImap <- list(ROImap_fs5@data,ROImap_fs5@atlases)
    } else if (n_vert==64984)
    {
      template="fslr32k"
      ROImap_fslr32k <- get('ROImap_fslr32k')
      ROImap <- list(ROImap_fslr32k@data,ROImap_fslr32k@atlases)
    } else if (n_vert==81924)
    {
      template="fsaverage6"
      ROImap_fs6 <- get('ROImap_fs6')
      ROImap <- list(ROImap_fs6@data,ROImap_fs6@atlases)
    } else if (n_vert==14524)
    {
      brainspace.mesh.mesh_io=reticulate::import("brainspace.mesh.mesh_io", delay_load = TRUE)
      template=brainspace.mesh.mesh_io$read_surface(paste0(system.file(package='VertexWiseR'),'/extdata/hip_template.fs'))
      ROImap_hip <- get('ROImap_hip')
      ROImap <- list(ROImap_hip@data,ROImap_hip@atlases)
    } else 
    {stop("data vector should only contain 20484 (fsaverage5), 81924 (fsaverage6), 64984 (fslr32k) or 14524 (hippocampal vertices) columns")}
  
  #Solves the "no visible binding for global variable" issue
  . <- SLM <- NULL 
    
  ##import python libaries
  reticulate::source_python(paste0(system.file(package='VertexWiseR'),'/python/brainstat.stats.SLM_VWR.py'))
  

  ##fitting model
  #preparing mask for model
  mask=array(rep(TRUE,NCOL(surf_data)))
  maskNA=which(colSums(surf_data != 0) == 0)
  mask[which(colSums(surf_data != 0) == 0)]=FALSE
  
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
  
  #define model to fit
  if(missing(random)) {model0=FixedEffect(model, "_check_categorical" = FALSE)}
  else {model0=MixedEffect(ran = as.factor(random),fix = model,"_check_categorical" = FALSE)}
  

  #read version of SLM that allows to specify the directory for the
  #fetch_template_surface option
  
  
  model=SLM(model = model0,
            contrast=contrast,
            surf = template,
            mask=mask,
            correction=c("fdr", "rft"),
            cluster_threshold=p,
            data_dir=data_dir)
  
  #fit will fetch parcellation data in a different place
  model$data_dir=paste0(brainstat_data_path,'/brainstat_data/parcellation_data/')
  
  #fit model
  SLM$fit(model,surf_data)
  
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
  
  #False discovery rate-corrected p value map
  fdrpmap=model$Q
  
  #listing objects to return
  returnobj=(list(cluster_results,as.numeric(model$t), as.numeric(tstat),as.numeric(posmask),as.numeric(negmask),as.numeric(pos_clusterIDmap),as.numeric(neg_clusterIDmap), as.numeric(bi_clusterIDmap), as.numeric(fdrpmap)))
  names(returnobj)=c("cluster_level_results","tstat_map","thresholded_tstat_map","pos_mask","neg_mask","pos_clusterIDmap","neg_clusterIDmap", "bi_clusterIDmap", "fdr_pmap")
  return(returnobj)
}
