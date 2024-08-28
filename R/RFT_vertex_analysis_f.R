
#' @title Vertex-wise analysis with random field theory cluster correction (formula)
#'
#' @description Alternative function to RFT_vertex_analysis() where a dataset and formula arguments can be provided instead of contrast or model.
#' @details  
#' 
#' Output definitions:
#' - `nverts`: number of vertices in the cluster
#' - `P`: p-value of the cluster
#' - `X, Y and Z`: MNI coordinates of the vertex with the highest t-statistic in the cluster.
#' - `tstat`: t statistic of the vertex with the highest t-statistic in the cluster
#' - `region`: the region this highest -statistic vertex is located in, as determined/labelled by the selected atlas 
#'
#' @param formula A string or formula object describing the predictors to be fitted against the surface data:
#' - The dependent variable can have any name, it will always be the surface data values. 
#' - The first independent variable in the formula will always be interpreted as the contrast of interest for which to estimate cluster-thresholded t-stat maps. 
#' - Only one random regressor can be given and must be indicated as '(1|variable_name)'.
#' @param dataset A data.frame object where the independent variables or predictors are stored (each IV's column names has to match the formula names).
#' @param surf_data A N x M matrix object containing the surface data (N row for each subject, M for each vertex), in fsaverage5 (20484 vertices), fsaverage6 (81924 vertices), fslr32k (64984 vertices) or hippocampal (14524 vertices) space. See also Hipvextract() or SURFvextract() output format. 
#' @param p A numeric object specifying the p-value to threshold the results (Default is 0.05)
#' @param atlas A numeric integer object corresponding to the atlas of interest. 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360, 5=Destrieux-148.
#' @param smooth_FWHM A numeric vector object specifying the desired smoothing width in mm 
#' @param VWR_check A boolean object specifying whether to check and validate system requirements. Default is TRUE.
#'
#' @returns A list object containing the cluster level results, thresholded t-stat map, and positive, negative and bidirectional cluster maps. Additionally, the formula and 'model' object based on the formula is included.
#' 
#' @seealso \code{\link{RFT_vertex_analysis}}
#' 
#' @examples
#' formula= as.formula("thickness ~ age + sex")
#' demodata = readRDS(system.file('demo_data/SPRENG_behdata_site1.rds', 
#'                               package = 'VertexWiseR'))
#' CTv = readRDS(file = url(paste0("https://github.com",
#' "/CogBrainHealthLab/VertexWiseR/blob/main/inst/demo_data/",
#' "SPRENG_CTv_site1.rds?raw=TRUE")))
#'
#' vertexwise_model=RFT_vertex_analysis_f(
#'   formula=formula, 
#'   dataset=demodata, 
#'   smooth_FWHM=10,
#'   surf_data = CTv, 
#'   p = 0.05,
#'   VWR_check=FALSE)
#' 
#' #Description of the output:
#' #vertexwise_model$cluster_level_results
#' 
#' @importFrom stringr str_split str_detect str_replace str_match
#' @importFrom stats lm model.matrix runif
#' @export

RFT_vertex_analysis_f=function(formula, dataset, surf_data, p=0.05, atlas=1, smooth_FWHM, VWR_check=TRUE) 
{
 
  #Check required python dependencies. If files missing:
  #Will prompt the user to get them in interactive session 
  #Will stop if it's a non-interactive session
  if (VWR_check == TRUE){
    message("Checking for VertexWiseR system requirements ... ")
    check = VWRfirstrun(n_vert=max(dim(t(surf_data))))
    if (!is.null(check)) {return(check)}
  } else if(interactive()==FALSE) { return(message('Non-interactive sessions need requirement checks'))}
  
  #If the dataset is a tibble (e.g., taken from a read_csv output)
  #converts the columns to regular data.frame column types
  if ('tbl_df' %in% class(dataset) == TRUE) {
    dataset=as.data.frame(dataset)
    if (NCOL(dataset)==1) {dataset = dataset[[1]]
    } else { for (c in 1:NCOL(dataset)) { 
      if(inherits(dataset[,c],"double")==TRUE) {dataset[,c] = as.numeric(dataset[,c])}
    }  }
  }
  
  #turns formula to string if formula object
  if (inherits(formula, 'formula'))
  {formula_str <- paste(deparse(formula), collapse = " ")
  } else if (inherits(formula, 'character'))
  {formula_str=formula
  } else
  {stop('The formula should be a formula object or a character string object.\n')}
  
  ###Run a dummy regression to define all variables from the
  ###formula (interactions, randoms etc.) and store them in a data.frame
  
  #get the DV name to create a dummy DV
  if (length(grep(x=formula_str, pattern="~"))==0)
  {stop('The formula needs a dependent variable defined.\n')}
  DV=stringr::str_split(formula_str, pattern = "~")[[1]][1]
  DV=gsub(" ", "", DV) #remove spaces
  data=cbind(runif(length(dataset[,1]), min=1, max=100),dataset); 
  colnames(data)[1]=DV;
  
  ##Run dummy model to create the data.frame
  #if has random variable, will be saved separately
  if (stringr::str_detect(formula_str, pattern="\\(1 \\|")==FALSE
      & stringr::str_detect(formula_str, pattern="\\(1\\|")==FALSE)
  { #dummy modelling of the formula 
    mod=lm(formula_str,data=data);
    #extracting variables as computed within formula
    model=model.matrix(mod)[,2:ncol(model.matrix(mod))] #without intercept
  } else
  {
  #extract the name of the random variable in (1|var) or equivalent
  pattern <- "\\s*\\(1\\s*\\|\\s*(.*?)\\s*\\)"
  random_var <- stringr::str_match(formula_str, pattern)[,2]

  #remove the random_var, run normal lm, and save it as the random object
  formula_str <- stringr::str_replace(formula_str, pattern, "")
  random=data[,random_var]
  colnames(random) <- random_var
  
  #run the model 
  mod=lm(formula_str,data=data);
  #extracting variables as computed within formula
  model=model.matrix(mod)[,2:ncol(model.matrix(mod))] #without intercept
  model=cbind.data.frame(model,)   #append random var

  }
  
  #The models do not accept variables with more than 2 levels
  #check from the formula and stop if one var has that issue
  if(length(mod$xlevels)!=0) ## if no categorical variables are entered in the formula, mod$xlevels will be empty
  {
   for (var in 1:length(mod$xlevels)) 
   {
     if (length(mod$xlevels[[var]]) > 2)
     {stop(paste("The categorical variable '", names(mod$xlevels[var]),"' contains more than 2 levels, please code it into binarized dummy variables.\n",sep=""))}
     else if (length(mod$xlevels[[var]]) == 2)
     {message(paste("The categorical variable '", names(mod$xlevels[var]),"' has been coded 0 for", mod$xlevels[[var]][1], "and 1 for", mod$xlevels[[var]][2], "respectively.\n"))}
   }
  }
  
  
  if (!exists('random_var')) { random = NULL }
  
  #Run the regular vertex analysis with the right objects
  RFTmodel=RFT_vertex_analysis(model=model,
                                     contrast=model[,1], #contrast is 1st IV
                                     random=random,
                                     surf_data=surf_data, 
                                     p=p, 
                                     atlas=atlas, 
                                     smooth_FWHM=smooth_FWHM, 
                                     VWR_check=FALSE)
  RFTmodel$formula=formula_str
  RFTmodel$model=model


return(RFTmodel)
}
