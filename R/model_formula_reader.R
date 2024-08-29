
#' @title Formula reader for RFT_vertex_analysis()
#'
#' @description Alternative function to RFT_vertex_analysis() where a dataset and formula arguments can be provided instead of contrast or model.
#' @param formula A string or formula object describing the predictors to be fitted against the surface data:
#' - The dependent variable is not needed, as it will always be the surface data values. 
#' - The first independent variable in the formula will always be interpreted as the contrast of interest for which to estimate cluster-thresholded t-stat maps. 
#' - Only one random regressor can be given and must be indicated as '(1|variable_name)'.
#' @param formula_dataset A data.frame object where the independent variables or predictors are stored (each IV's column names has to match the formula names).
#'
#' @returns A list object containing the model, contrast and random data.frame to be used in the RFT_vertex_analysis() modelling.
#' 
#' @seealso \code{\link{RFT_vertex_analysis}}
#' 
#' @examples
#' formula= as.formula("~ age + sex")
#' formula_dataset = readRDS(system.file('demo_data/SPRENG_behdata_site1.rds',
#'                               package = 'VertexWiseR'))
#' vertexwise_model=model_formula_reader(
#'   formula, 
#'   formula_dataset)
#' 
#' @importFrom stringr str_split str_detect str_replace str_match
#' @importFrom stats lm model.matrix runif as.formula na.pass 
#' @noRd

model_formula_reader=function(formula, formula_dataset) 
{
  #If the dataset is a tibble (e.g., taken from a read_csv output)
  #converts the columns to regular data.frame column types
  if ('tbl_df' %in% class(formula_dataset) == TRUE) {
    formula_dataset=as.data.frame(formula_dataset)
    if (NCOL(formula_dataset)==1) {formula_dataset = formula_dataset[[1]]
    } else { for (c in 1:NCOL(formula_dataset)) { 
      if(inherits(formula_dataset[,c],"double")==TRUE) {formula_dataset[,c] = as.numeric(formula_dataset[,c])}
    }  }
  }
  
  # Replace empty cells with NA
  formula_dataset[formula_dataset == ""] <- NA
  
  #turns formula to string if formula object
  if (inherits(formula, 'formula'))
  {formula_str <- paste(deparse(formula), collapse = " ")
  } else if (inherits(formula, 'character'))
  {formula_str=formula
  } else
  {stop('The formula should be a formula object or a character string object.\n')}
  
  ###Run a dummy regression to define all variables from the
  ###formula (interactions, randoms etc.) and store them in a data.frame
  
  #if there is no wave, adds it   
  if (length(grep(x=formula_str, pattern="~"))==0)
  {formula_str=paste('~',formula_str)}
  #if there is a DV defined in the formula, removes it
  DV=stringr::str_split(formula_str, pattern = "~")[[1]][1]
  formula_str=gsub(DV, "", formula_str) 
  #create dummy DV to formula as default
  formula_str=paste('surfmeasure',formula_str)
  data=cbind(runif(length(formula_dataset[,1]), min=1, max=100),formula_dataset); 
  colnames(data)[1]='surfmeasure';
  
  ##Run dummy model to create the data.frame
  #if has random variable, will be saved separately
  if (stringr::str_detect(formula_str, pattern="\\(1 \\|")==FALSE
      & stringr::str_detect(formula_str, pattern="\\(1\\|")==FALSE)
  { #dummy modelling of the formula 
    mod=lm(as.formula(formula_str),data=data);
    #extracting variables as computed within formula
    model=model.matrix(mod,na.action=na.pass)[,2:ncol(model.matrix(mod))] #without intercept
  } else
  {
  #extract the name of the random variable in (1|var) or equivalent
  pattern <- "\\s*\\(1\\s*\\|\\s*(.*?)\\s*\\)"
  random_var <- stringr::str_match(formula_str, pattern)[,2]

  #remove the random_var, run normal lm, and save it as the random object
  formula_str <- stringr::str_replace(formula_str, pattern, "")
  random=data[,random_var]
  
  #run the model 
  mod=lm(as.formula(formula_str),data=data);
  #extracting variables as computed within formula
  model=model.matrix(mod,na.action=na.pass)[,2:ncol(model.matrix(mod))] #without intercept
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
  
  #Run the regular vertex analysis with the right objects
    if (exists('random')) {
      RFTmodel=list(model=model,
                    contrast=data.matrix(model[,1]), 
                    random=random)
    } else {
      RFTmodel=list(model=model,
                    contrast=data.matrix(model[,1]))
    }


return(RFTmodel)
}
