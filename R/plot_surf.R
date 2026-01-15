#' @title Surface plotter
#'
#' @description Plots surface data in a grid with one or multiple rows in a .png file
#'
#' @param surf_data  A numeric vector (length of V) or a matrix (N rows x V columns), where N is the number of subplots, and V is the number of vertices. It can be the output from SURFvextract(), CAT12vextract(), FSLRvextract(), HIPvextract() as well as masks or vertex-wise results outputted by analyses functions. Alternatively, atlas ROI values as supported by atlas_to_surf() may be given.
#' @param filename A string object containing the desired name of the output .png. Default is 'plot.png' in the R temporary directory (tempdir()).Only filenames with a .png extension are allowed.
#' @param title A string object for setting the title in the plot. Default is none. For titles that too long to be fully displayed within the plot, we recommend splitting them into multiple lines by inserting "\\n".
#' @param surface A string object containing the name of the type of cortical surface background rendered. Possible options include "white", "smoothwm","pial" and "inflated" (default). The surface parameter is ignored for hippocampal surface data.
#' @param cmap A string object specifying the name of an existing colormap or a vector of hexadecimal color codes to be used as a custom colormap. The names of existing colormaps are listed in the \href{https://matplotlib.org/stable/gallery/color/colormap_reference.html}{'Matplotlib' plotting library}. 
#' 
#' Default cmap is set to `"Reds"` for positive values, `"Blues_r"` for negative values and `"RdBu"` when both positive and negative values exist. 
#' @param limits A combined pair of numeric vector composed of the lower and upper color scale limits of the plot. If the limits are specified, the same limits will be applied to all subplots. When left unspecified, the same symmetrical limits c(-max(abs(surf_dat),max(abs(surf_dat))) will be used for all subplots. If set to NULL, each subplot will have its own limits corresponding to their min and max values
#' @param colorbar A logical object stating whether to include a color bar in the plot or not (default is TRUE).
#' @param size A combined pair of numeric vector indicating the image dimensions (width and height in pixels). Default is c(1920,400) for whole-brain surface and c(400,200) for hippocampal surface.
#' @param zoom A numeric value for adjusting the level of zoom on the figures. Default is 1.25 for whole-brain surface and 1.20 for hippocampal surface.
#' @param transparent_bg A logical object to determine if the background of the image is set to transparent (Default is FALSE).
#' @param show_nan A logical object to determine if the NaN vertices are to be plotted (Default is TRUE).
#' @param alpha A numeric object between 0 and 1 to determine the opacity of the non-empty vertices. Note that this is not a true opacity setting, it will blend the colour into that of the NaN vertices (white if show_nan is FALSE)
#' @param show.plot.window A logical object to determine if the generated plot is to be shown within RStudio's plot window
#' @param VWR_check A boolean object specifying whether to check and validate system requirements. Default is TRUE.
#'
#' @returns Outputs the plot as a .png image
#' @examples
#' results = runif(20484,min=0, max=1);
#' plot_surf(surf_data = results, 
#' filename=paste0(tempdir(),"/output.png"),
#' title = 'Cortical thickness', 
#' surface = 'inflated', cmap = 'Blues',
#' VWR_check=FALSE)
#' @importFrom reticulate tuple import np_array source_python
#' @importFrom grDevices col2rgb
#' @importFrom png readPNG
#' @importFrom grid grid.raster
#' @export
######################################################################################################################################################
######################################################################################################################################################
plot_surf=function(surf_data, filename, title="",surface="inflated",cmap,limits, colorbar=TRUE, size, zoom, transparent_bg=FALSE, show_nan=T, alpha=1, show.plot.window=TRUE,VWR_check=TRUE)
{
  #Check required python dependencies. If files missing:
  #Will prompt the user to get them in interactive session 
  #Will stop if it's a non-interactive session 
  if (VWR_check == TRUE){
    message("Checking for VertexWiseR system requirements ...")
    check = VWRfirstrun(n_vert=max(dim(t(surf_data))))
    if (!is.null(check)) {return(check)}
  } else if(interactive()==FALSE) { return(message('Non-interactive sessions need requirement checks'))}
  
  if (missing("filename")) {
    message('No filename argument was given. The plot will be saved as "plot.png" in R temporary directory (tempdir()).\n')
    filename=paste0(tempdir(),'/plot.png')
  }
  
  #format title for single row
  if(is.null(nrow(surf_data)))
  {
    title=list('left'=list(title))
    rows=1
    surf_data=as.numeric(surf_data)
  } else 
  {
    if(nrow(surf_data)==1)
    {
      title=list('left'=list(title))
      rows=1
      surf_data=as.numeric(surf_data)
    }
    else 
    {
    rows=nrow(surf_data)
    }
  }
  
  #in a multi-row data scenario: insert a dummy title if title is missing  or repeat the title nrow times
  if(rows>1) 
  {
    if(missing("title")) {title=rep(NULL,rows)}
    else if (!missing("title")) {title=rep(title,rows)}
  }
  
  #check length of vector
  n_vert=length(surf_data)
  #if surface template is inputted
  if(n_vert%%20484==0) {template="fsaverage5"}
  else if (n_vert%%64984==0) {template="fslr32k"} 
  else if (n_vert%%81924==0) {template="fsaverage6"} 
  else if (n_vert%%14524==0) {template="CIT168"} 
  #if atlas object is inputted
  else if (max(dim(t(surf_data))) == 10) {template="CIT168";
  surf_data=atlas_to_surf(surf_data, template)} 
  else if (max(dim(t(surf_data))) %in% c(70,148,360,100,200,400)) {template="fsaverage5";
  surf_data=atlas_to_surf(surf_data, template)} 
  else{stop("surf_data vector should only contain 20484 (fsaverage5), 81924 (fsaverage6), 64984 (fslr32k) or 14524 (hippocampal vertices) columns. If you intended to plot an atlas' parcels, please refer to ?atlas_to_surf() for information about accepted atlases.")
  }
  
  #if cmap is missing, select cmaps depending on whether the image contains positive only or negative only values
  if(missing("cmap"))
  {
    if(range(surf_data,na.rm = TRUE)[1]>=0)  {cmap="Reds"}
    else if (range(surf_data,na.rm = TRUE)[2]<=0)  {cmap="Blues_r"}
    else  {cmap="RdBu_r"}  
  }
  
  #custom cmap— if a vector of hex color codes is specified
  if(length(cmap)>1)
  {
    #need to manually specify class as colors because the output of some color palette functions are not automatically assigned the class of colors
    class(cmap)="colors"
    matplotlib=reticulate::import("matplotlib", delay_load = TRUE)
    
    custom_colors=t(grDevices::col2rgb(cmap)/255) # convert hex color codes to RGB codes, then divide by 255 to convert to RGBA codes
    
    #save as python cmap object
    mymap = matplotlib$colors$LinearSegmentedColormap$from_list('my_colormap', custom_colors)
    matplotlib$colormaps$unregister(name = "custom_cmap")
    matplotlib$colormaps$register(cmap = mymap,name="custom_cmap")
    cmap="custom_cmap"  
  }
  
  #show nan or not
  if (show_nan==TRUE)
  {nan_color=reticulate::tuple(0.7, 0.7, 0.7, 1)} else 
  {nan_color=reticulate::tuple(0.7, 0.7, 0.7, 0)}
  
  #set opacity
  if (alpha<1 & alpha>=0)
  {
    #establish default NaN colour 
    if (show_nan==TRUE)
    {nan_rgb <- c(0.7, 0.7, 0.7)} else
    {nan_rgb <- c(1, 1, 1)}
    
    #read cmap
    matplotlib <- reticulate::import("matplotlib", delay_load = TRUE)
    np <- reticulate::import("numpy", delay_load = TRUE)
    cmap_obj <- matplotlib$cm$get_cmap(cmap)
    # read rgb colors across the colormap
    rgb_vals <- cmap_obj(np$linspace(0, 1, 256L))[, 1:3]

    #blend cmap toward NaN colour (not true opacity but works the same)
    blended <- (1 - alpha) * matrix(nan_rgb, nrow = nrow(rgb_vals), ncol = 3, byrow = TRUE) +
      alpha * rgb_vals
    
    #save as python cmap object
    mymap = matplotlib$colors$LinearSegmentedColormap$from_list('blended_map', blended)
    matplotlib$colormaps$unregister(name = "blended_map")
    matplotlib$colormaps$register(cmap = mymap,name="blended_map")
    cmap="blended_map"  
  }
  
  #setting color scale limits
  maxlimit=max(abs(range(surf_data,na.rm = TRUE)))
  if(missing("limits")) 
  {
    if(range(surf_data,na.rm = TRUE)[1]>=0) {limits=reticulate::tuple(0,range(surf_data,na.rm = TRUE)[2])} ##if image contains all positive values
    else if(range(surf_data,na.rm = TRUE)[2]<=0) {limits=reticulate::tuple(range(surf_data,na.rm = TRUE)[1],0)} ##if image contains all negative values
    else {limits=reticulate::tuple(-maxlimit,maxlimit)} ##symmetrical limits will be used if image contains both positive and negative values
  } else {
    ##user specified limits
    if(!is.null(limits))
    {
      limits=reticulate::tuple(limits[1],limits[2])  
    }   
  }
  
  if(template %in% c('fsaverage5','fsaverage6','fslr32k'))
  {
    ##cortical surface fplots
    #import python libraries
    brainstat.datasets=reticulate::import("brainstat.datasets", delay_load = TRUE)  
    brainspace.plotting=reticulate::import("brainspace.plotting", delay_load = TRUE)
    
    #For brainstat data, it will look either in default $HOME path or 
    #custom if it's been set
    # If custom installation paths have been defined by the user, source
    # them from the package directory:
    Renvironpath=paste0(tools::R_user_dir(package='VertexWiseR'),'/.Renviron')
    if (file.exists(Renvironpath)) {readRenviron(Renvironpath)}
    
    if (Sys.getenv('BRAINSTAT_DATA')=="")
    { 
      brainstat_data_path=fs::path_home()
    } 
    else if (!Sys.getenv('BRAINSTAT_DATA')=="") 
    {
      brainstat_data_path=Sys.getenv('BRAINSTAT_DATA')
    }
    
    #loading fsaverage surface
    left=brainstat.datasets$fetch_template_surface(template, join=FALSE, layer=surface,data_dir = paste0(brainstat_data_path,'/brainstat_data/surface_data/'))[1]
    right=brainstat.datasets$fetch_template_surface(template, join=FALSE, layer=surface,data_dir = paste0(brainstat_data_path,'/brainstat_data/surface_data/'))[2]
    
    #default cortical size and zoom parametes
    if(missing("size")) { size=c(1920,rows*400)}
    if(missing("zoom")) { zoom=1.25 }
    
    surf_plot=brainspace.plotting$plot_hemispheres(
      left[[1]], 
      right[[1]],  
      array_name=reticulate::np_array(surf_data),
      cmap=cmap,
      size=reticulate::tuple(as.integer(size)),
      nan_color=nan_color,
      return_plotter=TRUE,
      background=reticulate::tuple(as.integer(c(1,1,1))),
      zoom=zoom,
      color_range=limits,
      label_text=title,
      interactive=FALSE, 
      color_bar=colorbar,  
      transparent_bg=transparent_bg)  ##disabling interactive mode because this causes RStudio to hang
    
    
  #######Hippocampal template
  } else if (template=='CIT168')
  {
    #Solves the "no visible binding for global variable" issue
    . <- surfplot_canonical_foldunfold  <- NULL 
    internalenv <- new.env()
    assign("surfplot_canonical_foldunfold", surfplot_canonical_foldunfold, envir = internalenv)
    
    ##hippocampal plots
    #import python libraries
    reticulate::source_python(paste0(system.file(package='VertexWiseR'),'/python/hipp_plot.py'))
    
    #default hippocampal size and zoom parametes
    if(missing("size")) { size=c(400,200)}
    if(missing("zoom")) { zoom=1.2 }
    
    #reshaping surf_data into a 7262 x 2 x N array
    if(is.null(nrow(surf_data)))  {surf_data=cbind(surf_data[1:7262],surf_data[7263:14524])} #if N=1
    else  
    {
      if(nrow(surf_data)==1)
      {surf_data=cbind(surf_data[1:7262],surf_data[7263:14524])} 
      else
      {
      surf_data.3d=array(NA,c(7262,2,nrow(surf_data))) #if N>1
      for (row in 1:nrow(surf_data))  {surf_data.3d[,,row]=cbind(surf_data[row,1:7262],surf_data[row,7263:14524])}
      surf_data=surf_data.3d
      }
    }
    
    surf_plot=surfplot_canonical_foldunfold(
      surf_data,
      hipdat =get('hip_points_cells'),
      color_bar=colorbar,
      share="row",
      nan_color=nan_color,
      size=as.integer(size), 
      zoom=zoom,
      cmap=cmap,
      color_range=limits,
      label_text=title, 
      return_plotter=TRUE,
      interactive=FALSE) ##disabling interactive mode because this causes RStudio to hang
  }
  #output plot as a .png image
  surf_plot$screenshot(filename=filename,transparent_bg = transparent_bg)
  if(show.plot.window==TRUE)
  {
    grid::grid.newpage() #resets window if another plot was made
    img=png::readPNG(source =filename)
    grid::grid.raster(img)
  }
}
