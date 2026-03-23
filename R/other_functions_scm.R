#collection of functions for the analysis of subcortical analyses from SubCortexMesh

############################################################################################################################
############################################################################################################################

#' @title Database check for SubCortexMesh subcortical surface analyses
#'
#' @description Function to check whether the external subcortical template database has been downloaded in order to run plotting, smoothing and modelling functions on surface data derived from SubCortexMesh processing. It assists downloading the data from the github online repository. It is only triggered when a user inputs surf_data with the same number of vertices as SubCortexMesh-derived objects.
#' @details The database is downloaded and stored in the path indicated by system.file('extdata', package='VertexWiseR'). 
#' @param template A string indicating the subcortical segmentation scheme. Current options are 'fsaverage' and 'fslfirst'.
#' @examples
#' scm_database_check()
#' @importFrom utils menu
#' @noRd

scm_database_check=function(template){
  
  if (missing(template)) {
    stop('A template argument must be provided ("fsaverage" or "fslfirst").')
  }
  if (template=='fsaverage'){datsize='19.9'}
  else if (template=='fslfirst'){datsize='17.3'}
  else {stop(paste0(template,' is not a recognized template. Options are "fsaverage" or "fslfirst".'))}
  
  #package's external data path
  extdatadir=system.file('extdata', package='VertexWiseR')
  if (!dir.exists(paste0(extdatadir,'/scm_database/'))) {dir.create(paste0(extdatadir,'/scm_database/'))}
  dest_path=paste0(extdatadir,'/scm_database/',template,'.zip')
  
  #Check if cells data are present and download
  if (!dir.exists(paste0(extdatadir,'/scm_database/',template)))
  {
    missingobj=1
    
    prompt = utils::menu(c("Yes", "No"), title=paste0(
      "\nThe database directory was not detected inside VertexWiseR's installed package directory (",paste0(extdatadir,'/scm_database/',template,'/'), 
      "). It is needed to smooth, plot and analyse subcortical surfaces from SubCortexMesh. It can be downloaded from the github VertexWiseR directory.\n\nDo you want the scm_database database for ",template," (~",datsize," MB) to be downloaded now?"))
    if (prompt==1) {
      
      #function to check if url exists
      #courtesy of Schwarz, March 11, 2020, CC BY-SA 4.0:
      #https://stackoverflow.com/a/60627969
      valid_url <- function(url_in,t=2){
        con <- url(url_in)
        check <- suppressWarnings(try(open.connection(con,open="rt",timeout=t),silent=TRUE)[1])
        suppressWarnings(try(close.connection(con),silent=TRUE))
        ifelse(is.null(check),TRUE,FALSE)}
      
      #downloader 
      #Check if URL works and avoid returning error but only print message as requested by CRAN:
      url=paste0("https://raw.githubusercontent.com/CogBrainHealthLab/VertexWiseR/main/inst/extdata/scm_database/",template,".zip")
      if(valid_url(url)) {
        download.file(url=url,
                      destfile = dest_path)
      } else { 
        stop(paste0("The scm_database.zip failed to be downloaded from the github VertexWiseR directory. Please check your internet connection. Alternatively, you may visit https://github.com/CogBrainHealthLab/VertexWiseR/tree/main/inst/extdata/ and download the object manually.")) 
      }
      
      #unzipper
      out_dir=paste0(extdatadir,'/scm_database/')
      #unzip or untar as unzip struggles with long paths
      zip_result <- suppressWarnings(try(unzip(zipfile = dest_path, exdir = out_dir), silent = TRUE))
      if (is.null(zip_result) | inherits(zip_result, "try-error")==TRUE) {untar(tarfile = dest_path, exdir = out_dir)}
      unlink(dest_path) #clear zip file after extraction
      
      if (!dir.exists(out_dir))
      {
        stop(paste0("The scm_database.zip failed to be unzipped via unzip() or untar(). You might want to try manually."))   
      }
      
      #if user refuses, stops if required, just returns a message if optional at this stage
    } else {
      stop("\nSubCortexMesh subcortical data can only be smoothed, plotted and analysed with VertexWiseR's template database.\n") 
    } 
    
  } #silent if database exists
}

############################################################################################################################
############################################################################################################################

#' @title Subcortical plotting parameters
#'
#' @description For each subcortical region-of-interest (ROI), loads the corresponding points and cells data, and sets up default zoom, size for the plotter if not specified in plot_surf().
#' @param surf_data  A numeric vector (length of V) or a matrix (N rows x V columns), where N is the number of subplots, and V is the number of vertices. It can be the output of SCMvextract(), or masks outputted by vertex-wise analyses function that match the number of vertices of applicable ROIs.
#' @param size A combined pair of numeric vector indicating the image dimensions (width and height in pixels).
#' @param zoom A numeric value for adjusting the level of zoom on the figures. Default is 1.25 for whole-brain surface and 1.20 for hippocampal surface.
#' @param template A string indicating the subcortical segmentation scheme. Current options are 'fsaverage' and 'fslfirst'.
#' @param twocbars A boolean object stating whether the plot will include two colour bars as per plot_overlay_surf() or not. It matters for the size to encompass them properly for some regions. 
#' @returns Cell data for the left and right hemispheres (only one for the Brain-Stem), and the size and zoom parameters.
#' @noRd

scm_plot_parameters=function(surf_data,size,zoom,template,twocbars=FALSE)
{
  n_vert=max(dim(t(surf_data)))
  
  ###FSAVERAGE ASEG
  if (template=='fsaverage'){
    #accumbens nuclei
    if (n_vert==2044)
    {
      lh_vert=1022; rh_vert=1022 
      if(is.null(size) & twocbars==FALSE) {size=c(600,400)}
      if(is.null(size) & twocbars==TRUE) {size=c(800,400)}
      if(is.null(zoom)) {zoom=1.4}
    } 
    #amygdalae
    else if (n_vert==3430) 
    {
      lh_vert=1638; rh_vert=1792 
      if(is.null(size) & twocbars==FALSE) {size=c(600,400)}
      if(is.null(size) & twocbars==TRUE) {size=c(800,400)}
      if(is.null(zoom)) {zoom=1.4}
    }
    #caudate
    else if (n_vert==6940)
    {
      lh_vert=3440; rh_vert=3500     
      if(is.null(size) & twocbars==FALSE) {size=c(500,500)}
      if(is.null(size) & twocbars==TRUE) {size=c(500,400)}
      if(is.null(zoom)) {zoom=2.5}
    }
    #cerebella
    else if (n_vert==39214) 
    {
      lh_vert=19559; rh_vert=19664
      if(is.null(size) & twocbars==FALSE) {size=c(400,400)}
      if(is.null(size) & twocbars==TRUE) {size=c(600,400)}
      if(is.null(zoom)) {zoom=1.7}
    }
    #hippocampi
    else if (n_vert==8132) 
    {
      lh_vert=4046; rh_vert=4086 
      if(is.null(size) & twocbars==FALSE) {size=c(400,400)}
      if(is.null(size) & twocbars==TRUE) {size=c(600,400)}
      if(is.null(zoom)) {zoom=1.7}
    }
    #pallidum
    else if (n_vert==3200) 
    {
      lh_vert=1600; rh_vert=1600 
      if(is.null(size) & twocbars==FALSE) {size=c(400,400)}
      if(is.null(size) & twocbars==TRUE) {size=c(700,400)}
      if(is.null(zoom)) {zoom=1.7}
    }
    #putamen
    else if (n_vert==8394) 
    {
      lh_vert=4268; rh_vert=4126 
      if(is.null(size) & twocbars==FALSE) {size=c(400,400)}
      if(is.null(size) & twocbars==TRUE) {size=c(700,400)}
      if(is.null(zoom)) {zoom=1.7}
    }
    #thalamus
    else if (n_vert==7768) 
    {
      lh_vert=3936; rh_vert=3832 
      if(is.null(size) & twocbars==FALSE) {size=c(600,400)}
      if(is.null(size) & twocbars==TRUE) {size=c(600,400)}
      if(is.null(zoom)) {zoom=1.9}
    }
    #ventral diencephala
    else if (n_vert==7144) 
    {
      lh_vert=3550; rh_vert=3594 
      if(is.null(size) & twocbars==FALSE) {size=c(400,400)}
      if(is.null(size) & twocbars==TRUE) {size=c(600,400)}
      if(is.null(zoom)) {zoom=1.9}
    }
    #brain stem
    else if (n_vert==9452) 
    {
      lh_vert=9452; rh_vert=9452 
      if(is.null(size) & twocbars==FALSE) {size=c(400,400)}
      if(is.null(size) & twocbars==TRUE) {size=c(500,400)}
      if(is.null(zoom)) {zoom=1.9}
    }
    #all aseg
    else if (n_vert==95718) 
    {
      lh_vert=95718; rh_vert=95718 
      if(is.null(size)) {size=c(900,400)}
      if(is.null(zoom)) {zoom=1.55}
    }
    else
    {
      stop('the surf_data object given does not match the number of vertices of any subcortical template surfaces.')
    }
  }
  ###FSLFIRST
  else if (template=='fslfirst')
  {
    #accumbens nuclei
    if (n_vert==2026)
    {
      lh_vert=1104; rh_vert=922 
      if(is.null(size) & twocbars==FALSE) {size=c(400,300)}
      if(is.null(size) & twocbars==TRUE) {size=c(800,400)}
      if(is.null(zoom)) {zoom=1.4}
    } 
    #amygdalae
    else if (n_vert==3592) 
    {
      lh_vert=1834; rh_vert=1758 
      if(is.null(size) & twocbars==FALSE) {size=c(500,400)}
      if(is.null(size) & twocbars==TRUE) {size=c(800,400)}
      if(is.null(zoom)) {zoom=1.4}
    }
    #caudate
    else if (n_vert==7570)
    {
      lh_vert=3600; rh_vert=3970     
      if(is.null(size) & twocbars==FALSE) {size=c(500,500)}
      if(is.null(size) & twocbars==TRUE) {size=c(500,400)}
      if(is.null(zoom)) {zoom=2.5}
    }
    #cerebella
    else if (n_vert==31466) 
    {
      lh_vert=15806; rh_vert=15660
      if(is.null(size) & twocbars==FALSE) {size=c(400,400)}
      if(is.null(size) & twocbars==TRUE) {size=c(600,400)}
      if(is.null(zoom)) {zoom=1.7}
    }
    #hippocampi
    else if (n_vert==8244) 
    {
      lh_vert=4044; rh_vert=4200 
      if(is.null(size) & twocbars==FALSE) {size=c(400,400)}
      if(is.null(size) & twocbars==TRUE) {size=c(600,400)}
      if(is.null(zoom)) {zoom=1.7}
    }
    #pallidum
    else if (n_vert==3548) 
    {
      lh_vert=1778; rh_vert=1770 
      if(is.null(size) & twocbars==FALSE) {size=c(400,400)}
      if(is.null(size) & twocbars==TRUE) {size=c(700,400)}
      if(is.null(zoom)) {zoom=1.7}
    }
    #putamen
    else if (n_vert==7908) 
    {
      lh_vert=3978; rh_vert=3930 
      if(is.null(size) & twocbars==FALSE) {size=c(400,400)}
      if(is.null(size) & twocbars==TRUE) {size=c(700,400)}
      if(is.null(zoom)) {zoom=1.7}
    }
    #thalamus
    else if (n_vert==8542) 
    {
      lh_vert=4316; rh_vert=4226 
      if(is.null(size) & twocbars==FALSE) {size=c(400,500)}
      if(is.null(size) & twocbars==TRUE) {size=c(600,400)}
      if(is.null(zoom)) {zoom=1.9}
    }
    #brain stem
    else if (n_vert==9516) 
    {
      lh_vert=9516; rh_vert=9516 
      if(is.null(size) & twocbars==FALSE) {size=c(400,400)}
      if(is.null(size) & twocbars==TRUE) {size=c(500,400)}
      if(is.null(zoom)) {zoom=1.9}
    }
    #all aseg
    else if (n_vert==82412) 
    {
      lh_vert=82412; rh_vert=82412 
      if(is.null(size)) {size=c(900,400)}
      if(is.null(zoom)) {zoom=1.55}
    }
    else
    {
      stop('the surf_data object given does not match the number of vertices of any subcortical template surfaces.')
    }
  } else {stop(paste0(template,' is not a recognized template. Options are "fsaverage" or "fslfirst".'))}
  
  
  #get cell data
  celldat=scm_database_fetcher(n_vert,'points_cells',template)
  lh_celldat_path=readRDS(celldat$lh_celldat_path)
  rh_celldat_path=readRDS(celldat$rh_celldat_path)
  
  return(list(lh_celldata=lh_celldat_path, lh_vert=lh_vert,
              rh_celldata=rh_celldat_path, rh_vert=rh_vert,
              size=size, zoom=zoom))
}


############################################################################################################################
############################################################################################################################

#' @title SubCortexMesh database fetcher
#'
#' @description For each subcortical region-of-interest (ROI), loads the corresponding edge list
#' @param n_vert A single number or integer, representing the number of vertices contained in the template surface data.
#' @param data_object A string object stating the data object to be extracted from the scm_database ('edgelist', 'ROImap', 'template', or 'points_cells').
#' @param template A string indicating the subcortical segmentation scheme. Current options are 'fsaverage' and 'fslfirst'.
#' @details #' Number of vertices in the (bilateral) matrix for each region-of-interest:
#'  - fsaverage/fslfirst: ROI
#'  - 2044/2026: accumbens area
#'  - 3430/3592: amygdala
#'  - 6940/7570: caudate nuclei
#'  - 39214/31466: cerebellum
#'  - 8132/8244: hippocampus
#'  - 3200/3548: pallidum
#'  - 8394/7908: putamen
#'  - 7768/8542: thalamus
#'  - 7144/NA: ventral diencephalon
#'  - 9452/9516: brain stem
#'  - 95718/82412: all ROIs merged
#'  
#' @returns If data_object=='edgelist': A N x 2 matrix object listing each vertex of the template and the vertices adjacent to it (making an edge together).
#' If data_object=='ROImap': A matrix object with N vertices from a template and each parcellation number the vertices correspond to in 6 atlases (6 columns).
#' If data_object=='template': A string object indicating the path to the ROI's surface template in the scm_database for reticulate to retrieve it.
#' If data_object=='points_cells': A string object indicating the path to the ROI's points and cells data in the scm_database.
#' @noRd

scm_database_fetcher=function(n_vert, data_object, template){
  
  if (missing(template)) {
    stop('A template argument must be provided ("fsaverage" or "fslfirst").')
  }
  
  #preset template dimensions
  if (template=='fsaverage'){
    #fsaverage
    if (n_vert==2044) {roilabel='accumbens'}
    else if (n_vert==3430) {roilabel='amygdala'}
    else if (n_vert==6940) {roilabel='caudate'}
    else if (n_vert==39214) {roilabel='cerebellum'}
    else if (n_vert==8132) {roilabel='hippocampus'}
    else if (n_vert==3200) {roilabel='pallidum'}
    else if (n_vert==8394) {roilabel='putamen'}
    else if (n_vert==7768) {roilabel='thalamus'}
    else if (n_vert==7144) {roilabel='ventraldc'}
    else if (n_vert==9452) {roilabel='brainstem'}
    else if (n_vert==95718) {roilabel='allaseg'}
  }
  else if (template=='fslfirst'){
    #fslfirst
    if (n_vert==2026) {roilabel='accumbens'}
    else if (n_vert==3592) {roilabel='amygdala'}
    else if (n_vert==7570) {roilabel='caudate'}
    else if (n_vert==31466) {roilabel='cerebellum'}
    else if (n_vert==8244) {roilabel='hippocampus'}
    else if (n_vert==3548) {roilabel='pallidum'}
    else if (n_vert==7908) {roilabel='putamen'}
    else if (n_vert==8542) {roilabel='thalamus'}
    else if (n_vert==9516) {roilabel='brainstem'}
    else if (n_vert==82412) {roilabel='allfslfirst'}
  } else {stop(paste0(template,' is not a recognized template. Options are "fsaverage" or "fslfirst".'))}
  #the "database" path's existence should be checked by scm_database_check() already
  if(data_object=='edgelist')
  {
    database=system.file(paste0('extdata/scm_database/',template,'/roi_edgelists'), package='VertexWiseR')
    load(paste0(database,'/edgelist_',roilabel,'.rdata'))
    edgelist <- get(paste0('edgelist_',roilabel))@data
    return(edgelist)
  }
  
  if(data_object=='ROImap')
  {
    database=system.file(paste0('extdata/scm_database/',template,'/roi_ROImaps'), package='VertexWiseR')
    load(paste0(database,'/ROImap_',roilabel,'.rdata'))
    ROImap <- get(paste0('ROImap_',roilabel))
    return(ROImap)
  }
  
  if(data_object=='template')
  { 
    database=system.file(paste0('extdata/scm_database/',template,'/roi_templates'), package='VertexWiseR')
    templatepath=paste0(database, '/', roilabel,'.vtk')
    return(templatepath)
  }
  
  if(data_object=='points_cells')
  {
    database=system.file('extdata/scm_database/',template,'/roi_points_cells', 
                         package='VertexWiseR')
    if (roilabel=='brainstem') 
    {lh_celldat_path=paste0(database,"/",roilabel,"_points_cells.rds")
    rh_celldat_path=paste0(database,"/",roilabel,"_points_cells.rds")
    } else if (roilabel=='allaseg' | roilabel=='allfslfirst')
    {lh_celldat_path=paste0(database,"/",roilabel,"_points_cells.rds")
    rh_celldat_path=paste0(database,"/",roilabel,"_points_cells.rds")
    } else {
      lh_celldat_path=paste0(database,"/left-",roilabel,"_points_cells.rds")
      rh_celldat_path=paste0(database,"/right-",roilabel,"_points_cells.rds")
    }
    return(list(lh_celldat_path=lh_celldat_path, 
                rh_celldat_path=rh_celldat_path))
  }
}