#' @title Decode surface data
#'
#' @description Correlates the significant clusters of an earlier vertex-wise analysis with a database of task-based fMRI and voxel-based morphometric statistical maps and associate them with relevant key words. Decoding currently works with surfaces in fsaverage5 space only."
#'
#' @details The \href{https://nimare.readthedocs.io/en/stable/index.html}{'NiMARE'} python module is used for the imaging decoding and is imported via the reticulate package. The function also downloads the \href{https://github.com/neurosynth/neurosynth-data}{'Neurosynth' database} in the package's inst/extdata directory (~8 Mb) for the analysis.
#'
#' @param surf_data A numeric vector or object containing the surface data, in fsaverage5 (1 x 20484 vertices) or fsLR32k (1 x 64984 vertices) space. It can only be one row of vertices (not a cohort surface data matrix). 
#' @param contrast A string object indicating whether to decode the positive or negative mask ('positive' or 'negative')
#' @param VWR_check A boolean object specifying whether to check and validate system requirements. Default is TRUE.
#'
#' @returns A data.frame object listing the keywords and their Pearson's R values
#' @examples
#' CTv = rbinom(20484, 1, 0.001) 
#' decoding = decode_surf_data(CTv, 'positive', VWR_check=FALSE);
#' head(decoding)
#' @importFrom reticulate import r_to_py
#' @export

##CT image decoding
decode_surf_data=function(surf_data,contrast="positive", VWR_check=TRUE) 
{
  
  #Check required python dependencies. If files missing:
  #Will prompt the user to get them in interactive session 
  #Will stop if it's a non-interactive session
  if (VWR_check == TRUE){
    message("Checking for VertexWiseR system requirements ... ")
    check = VWRfirstrun(requirement="neurosynth")
    if (!is.null(check)) {return(check)} 
  } else if(interactive()==FALSE) { return(message('Non-interactive sessions need requirement checks'))}
  
  
  # Check if all values are positive
  if  (all(surf_data >= 0)==TRUE & contrast=="negative")
  {stop('No negative cluster was identified in the surf_data.')}
  # Check if all values are negative
  if  (all(surf_data <= 0)==TRUE & contrast=="positive")
  {stop('No positive cluster was identified in the surf_data.')}
  
  
  #if neurosynth database is installed
  if(file.exists(system.file('extdata','neurosynth_dataset.pkl.gz', package='VertexWiseR'))==TRUE)
  {
    ##checks length
    if(is.vector(surf_data)) {n_vert=length(surf_data)} else {n_vert=ncol(surf_data)}
    if(n_vert==20484) {template="fsaverage5"}
    else if (n_vert==64984) {template="fslr32k"}
    else {stop("Only an surf_data vector with a length of 20484 (fsaverage5) or 64984 (fslr32k) is accepted")}
    
    #check contrast
    if(contrast != "positive" & contrast != "negative")  {stop("contrast has to be either positive or negative")} 
    
    message("Converting and interpolating the surface data ... ")
    
    ##import python libraries
    interpolate=reticulate::import("brainstat.mesh.interpolate", delay_load = TRUE)
    discrete=reticulate::import("nimare.decode", delay_load = TRUE)
    nimare.dataset=reticulate::import("nimare.dataset", delay_load = TRUE)
    
    ##selecting contrasts
    if(contrast=="positive")
    {
      surf_data[is.na(surf_data)]=0
      surf_data[surf_data<0]=0
      surf_data[surf_data>0]=1
    } else if (contrast=="negative")
    {
      surf_data[is.na(surf_data)]=0
      surf_data[surf_data>0]=0
      surf_data[surf_data<0]=1
    }
    
    ##convert surf_data vector to nii image
    stat_labels=reticulate::r_to_py(surf_data)
    stat_nii = interpolate$`_surf2vol`(template, stat_labels)
    
    
    ##running the decoding procedure
    neurosynth_dset = nimare.dataset$Dataset$load(system.file("extdata/neurosynth_dataset.pkl.gz", package='VertexWiseR'))
    message("\u2713 \n Correlating input image with images in the neurosynth database. This may take a while ... ")
    decoder = discrete$ROIAssociationDecoder(stat_nii)
    decoder$fit(neurosynth_dset)
    
    ##compiling the results
    decoder_df = data.matrix(decoder$transform())
    row.names(decoder_df)=gsub(pattern = "terms_abstract_tfidf__",x=row.names(decoder_df), replacement = "")
    result=data.frame(row.names(decoder_df),round(as.numeric(decoder_df),3))
    colnames(result)=c("keyword","r")
    result=result[order(-result$r),]
    message("\u2713 \n")
    return(result)
  } 
}  