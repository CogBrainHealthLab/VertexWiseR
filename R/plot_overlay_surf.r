#' @title Surface overlay plotter
#'
#' @description Plots surface data in a grid with one or multiple rows in a .png file
#' @param model_output A list object outputted by RFT_vertex_analysis() or TFCE_threshold(). The 'tstat_map' will automatically be treated as background map, and the 'thresholded_tstat_map' as overlay map. See surf_data_1 or surf_data_2 to assign any map manually.
#' @param surf_data_1 A numeric vector (length of V), where V is the number of vertices. It can be one row from the output from SURFvextract(), FSLRvextract(), HIPvextract() as well as masks or vertex-wise results outputted by analyses functions. This is the background surface.
#' @param surf_data_2 Same as surf_data_1. This is the overlay surface.
#' @param cmap_1 A string object specifying the name of an existing colormap or a vector of hexadecimal color codes to be used as a custom colormap. The names of existing colormaps are listed in the \href{https://matplotlib.org/stable/gallery/color/colormap_reference.html}{'Matplotlib' plotting library}. 
#' 
#' Default cmap is set to `"Reds"` for positive values, `"Blues_r"` for negative values and `"RdBu"` when both positive and negative values exist. This is meant for the background surface. 
#' @param cmap_2 Same as cmap_1. This is meant for the overlay surface. Default is the same as cmaps_1, adapted to limits of surf_data_2.
#' @param limits_1 A combined pair of numeric vector composed of the lower and upper color scale limits of surf_data_1. When left unspecified, the  symmetrical limits c(-max(abs(surf_data_1),max(abs(surf_data_1))) will be used. If set to NULL, the limits will correspond to the min and max values of surf_data_1. This is meant for the background surface.
#' @param limits_2 Same  as limits_1, or also the string object "same". Default is symmetrical limits relative to surf_data_2, while 'same' will apply the exact same limits as for limits_1. This is meant for the overlay surface.
#' @param colorbar_1 A logical object stating whether to include the color bar for the background layer in the plot or not (default is TRUE).
#' @param colorbar_2 A logical object stating whether to include the color bar for the overlay layer in the plot or not (default is TRUE).
#' @param alpha_1 A numeric object between 0 and 1 to determine the opacity of the non-empty vertices. Note that this is not a true opacity setting, it will blend the colour into that of the NaN vertices (white if show_nan is FALSE). This is meant for the background surface.
#' @param alpha_2 Same as alpha_2. This is meant for the overlay surface. Default is the same as alpha_1.
#' @param overlay_boundaries A logical object stating whether to plot black contour of the overlay layer.
#' @param show_nan A logical object to determine if the NaN vertices are to be plotted (Default is TRUE). This is meant for the background surface. The overlay surface will always omit NaN vertices to make the background visible.
#' @param filename A string object containing the desired name of the output .png. Default is 'combined_plots.png' in the R temporary directory (tempdir()).Only filenames with a .png extension are allowed.
#' @param title A string object for setting the title in the plot. Default is none. For titles that too long to be fully displayed within the plot, we recommend splitting them into multiple lines by inserting "\\n".
#' @param surface A string object containing the name of the type of cortical surface background rendered. Possible options include "white", "smoothwm","pial" and "inflated" (default). The surface parameter is ignored for hippocampal surface data.
#' @param size A combined pair of numeric vector indicating the image dimensions (width and height in pixels). Default is c(1920,400) for whole-brain surface and c(400,200) for hippocampal surface.
#' @param zoom A numeric value for adjusting the level of zoom on the figures. Default is 1.25 for whole-brain surface and 1.20 for hippocampal surface.
#' @param transparent_bg A logical object to determine if the background of the image is set to transparent (Default is FALSE).
#' @param show.plot.window A logical object to determine if the generated plot is to be shown within RStudio's plot window
#' @param VWR_check A boolean object specifying whether to check and validate system requirements. Default is TRUE.
#'
#' @returns Outputs the plot as a .png image
#' @examples
#' #simulate t-map
#' background=runif(20484,NA);
#' ROImap_fs5 <- get('ROImap_fs5')
#' ROImap <- list(ROImap_fs5@data,ROImap_fs5@atlases)
#' neg_parts <- c(1:100)
#' idx_neg <- which(ROImap[[1]][,1] %in% neg_parts)
#' background[idx_neg]=runif(length(idx_neg),-1,0)
#' pos_parts <- c(18:46)
#' idx_pos <- which(ROImap[[1]][,1] %in% pos_parts)
#' background[idx_pos]=runif(length(idx_neg),0,1)
#' #simulate clusters (randomly picked fsaverage parcels)
#' idx_sig=which(ROImap[[1]][,1]==29|
#'                 ROImap[[1]][,1]==19|
#'                                 ROImap[[1]][,1]==45)
#'                                 overlay = runif(20484,NA);
#'                                 overlay[idx_sig]=1
#'                                 
#'plot_overlay_surf(surf_data_1=background,
#'                 surf_data_2=overlay,
#'                 cmap_1='RdBu_r', cmap_2='Reds',
#'                 alpha_1=0.4, alpha_2=1,
#'                 filename=paste0(tempdir(),"/simulated_plot.png"), 
#'                title="Significant effects",
#'                overlay_boundaries=TRUE,
#'                VWR_check=FALSE)
#' @importFrom reticulate tuple import np_array source_python
#' @importFrom grDevices col2rgb
#' @importFrom png readPNG writePNG
#' @importFrom grid grid.raster
#' @importFrom stringr str_extract 
#' @export
######################################################################################################################################################
######################################################################################################################################################
plot_overlay_surf=function(model_output=NULL,
                           surf_data_1=NULL, surf_data_2=NULL, 
                           cmap_1, cmap_2,
                           limits_1, limits_2,
                           alpha_1=1, alpha_2=1,
                           colorbar_1=TRUE,
                           colorbar_2=TRUE,
                           show_nan=TRUE, 
                           filename, 
                           title="",
                           surface="inflated", 
                           overlay_boundaries=FALSE,
                           size, 
                           zoom, 
                           transparent_bg=FALSE, 
                           show.plot.window=TRUE,
                           VWR_check=TRUE)
{
  
  
  #fetch thresholded and unthresholded tstat map automatically
  if(!missing(model_output))
  {
    if('tstat_map' %in% names(model_output))
    {
      surf_data_1=model_output$tstat_map
    } else {stop('No tstat_map was found in your model_output variable. Make sure you are inputting the output of RFT_vertex_analysis() or TFCE_threshold(), or use the surf_data_1 and surf_data_2 arguments instead of the model_output argument.')}
    if('thresholded_tstat_map' %in% names(model_output))
    {
      surf_data_2=model_output$thresholded_tstat_map
    } else {stop('No thresholded_tstat_map was found in your model_output variable. Make sure you are inputting the output of RFT_vertex_analysis() or TFCE_threshold(), or use the surf_data_1 and surf_data_2 arguments instead of the model_output argument.')}
  } else {
    
    # If both model_output and one of surf_data_1/2 are NULL, throw error
    if (is.null(surf_data_1) | is.null(surf_data_2)) {
      stop("If no model_output is given, both surf_data_1 and surf_data_2 arguments should be provided.")
    }
  }
  
  #Check required python dependencies. If files missing:
  #Will prompt the user to get them in interactive session 
  #Will stop if it's a non-interactive session 
  if (VWR_check == TRUE){
    message("Checking for VertexWiseR system requirements ...")
    check = VWRfirstrun(n_vert=max(dim(t(surf_data_1))))
    if (!is.null(check)) {return(check)}
  } else if(interactive()==FALSE) { return(message('Non-interactive sessions need requirement checks'))}
  
  if (missing("filename")) {
    message('No filename argument was given. The plot will be saved as "overlay_plot.png" in R temporary directory (tempdir()).\n')
    filename=paste0(tempdir(),'/combined_plots.png')
  }
  
  #Both surf_data MUST have the same size
  if(is.null(nrow(surf_data_2)))
  {
    if (length(surf_data_1)!=length(surf_data_2)) 
    {stop('Surface vectors should have the exact same length.')}
  } else {
    if(dim(surf_data_1)[1]!=dim(surf_data_2)[1]|
       dim(surf_data_1)[2]!=dim(surf_data_2)[2])
    {stop('Surface matrices should have the exact same dimensions.')}
  }
  
  #format title for single row
  if(is.null(nrow(surf_data_1)))
  {
    title <- reticulate::dict(left = list(title)); surf_data=as.numeric(surf_data_1)
  } else {
    if(nrow(surf_data_1)){stop('Only one set of surfaces (background and overlay) can be plotted at a time. surf_data_1 and surf_data_2 should be individual vectors and not multiple-row matrices.')}
  }
  
  #check length of vector (only check surf_data_1 since they both must have the same dimensions anyways)
  n_vert=length(surf_data_1)
  #if surface template is inputted
  if(n_vert%%20484==0) {template="fsaverage5"
  } else if (n_vert%%64984==0) {template="fslr32k"
  } else if (n_vert%%81924==0) {template="fsaverage6"
  } else if (n_vert%%14524==0) {template="CIT168"
  #if atlas object is inputted
  } else if (max(dim(t(surf_data_1))) == 10) {template="CIT168";
  surf_data_1=atlas_to_surf(surf_data_1, template)
  } else if (max(dim(t(surf_data_1))) == 70) {template="fsaverage5";
  surf_data_1=atlas_to_surf(surf_data_1, template)
  } else if (max(dim(t(surf_data_1))) == 148) {template="fsaverage5";
  surf_data_1=atlas_to_surf(surf_data_1, template)
  } else if (max(dim(t(surf_data_1))) == 360) {template="fsaverage5";
  surf_data_1=atlas_to_surf(surf_data_1, template)
  } else if (max(dim(t(surf_data_1))) == 100) {template="fsaverage5";
  surf_data_1=atlas_to_surf(surf_data_1, template)
  } else if (max(dim(t(surf_data_1))) == 200) {template="fsaverage5";
  surf_data_1=atlas_to_surf(surf_data_1, template)
  } else if (max(dim(t(surf_data_1))) == 400) {template="fsaverage5";
  surf_data_1=atlas_to_surf(surf_data_1, template)
  } else{stop("surf_data vectors should only contain 20484 (fsaverage5), 81924 (fsaverage6), 64984 (fslr32k) or 14524 (hippocampal vertices) columns. If you intended to plot an atlas' parcels, please refer to ?atlas_to_surf() for information about accepted atlases.")
  }
  
  #NA-ify vertices if alpha=0 before computing color maps
  if (alpha_1==0){surf_data_1=rep(NaN,length(surf_data_1))
  } else if (alpha_2==0){surf_data_2=rep(NaN,length(surf_data_2))}
  
  #if cmap is missing, select cmaps depending on whether the image contains positive only or negative only values
  if(missing("cmap_1"))
  {
    if(range(surf_data_1,na.rm = TRUE)[1]>=0)  {cmap_1="Reds"
    } else if (range(surf_data_1,na.rm = TRUE)[2]<=0)  {cmap_1="Blues_r"
    } else  {cmap_1="RdBu_r"}  
  }
  if(missing("cmap_2"))
  {
    if(range(surf_data_2,na.rm = TRUE)[1]>=0)  {cmap_2="Reds"
    } else if (range(surf_data_2,na.rm = TRUE)[2]<=0)  {cmap_2="Blues_r"
    } else  {cmap_2="RdBu_r"}
  }
  
  #custom cmap— if a vector of hex color codes is specified
  #wrapped in a function to apply for both cmap_1 and cmap_2
  custom_map=function(cmap)
  {
    # Get the number of the variable passed (1 for "cmap_1" etc.)
    cmap_name <- deparse(substitute(cmap))
    suffix <- stringr::str_extract(cmap_name, "\\d+")
    # Construct unique colormap name
    custom_name <- paste0("custom_map_", suffix)
    
    #need to manually specify class as colors because the output of some color palette functions are not automatically assigned the class of colors
    class(cmap)="colors"
    matplotlib=reticulate::import("matplotlib", delay_load = TRUE)
    
    custom_colors=t(grDevices::col2rgb(cmap)/255) # convert hex color codes to RGB codes, then divide by 255 to convert to RGBA codes
    
    #save as python cmap object
    mymap = matplotlib$colors$LinearSegmentedColormap$from_list('my_colormap', custom_colors)
    matplotlib$colormaps$unregister(name = custom_name)
    matplotlib$colormaps$register(cmap = mymap,name=custom_name)
    cmap=custom_name 
    return(cmap)
  }
  if(length(cmap_1)>1){cmap_1=custom_map(cmap_1)}
  if(length(cmap_2)>1){cmap_2=custom_map(cmap_2)}
  
#set opacity
 #wrapped in a function to apply for both alpha_1 and alpha_2
 alpha_function=function(alpha,show_nan){
    
    # Get the number of the variable passed (1 for "alpha_1" etc.)
    alpha_name <- deparse(substitute(alpha))
    suffix <- stringr::str_extract(alpha_name, "\\d+")
  
    if (alpha<1 & alpha>=0)
    {
      #establish default NaN colour 
      if (show_nan==TRUE) {nan_rgb <- c(0.7, 0.7, 0.7)} else{nan_rgb <- c(1, 1, 1)}
      
      #read cmap
      matplotlib <- reticulate::import("matplotlib", delay_load = TRUE)
      np <- reticulate::import("numpy")
      cmap_obj <- matplotlib$cm$get_cmap(get(paste0('cmap_',suffix)))
      # read rgb colors across the colormap
      rgb_vals <- cmap_obj(np$linspace(0, 1, 256L))[, 1:3]
      
      #blend cmap toward NaN colour (not true opacity but works the same)
      blended <- (1 - alpha) * matrix(nan_rgb, nrow = nrow(rgb_vals), ncol = 3, byrow = TRUE) + alpha * rgb_vals
      
      #save as python cmap object
      mymap = matplotlib$colors$LinearSegmentedColormap$from_list('blended_map', blended)
      matplotlib$colormaps$unregister(name = "blended_map")
      matplotlib$colormaps$register(assign(paste0("cmap_", suffix), mymap),
                                    name="blended_map")
      assign(paste0("cmap_", suffix), "blended_map")
    }
    
    return(get(paste0('cmap_',suffix)))
  }
  cmap_1=alpha_function(alpha_1,show_nan)
  cmap_2=alpha_function(alpha_2,show_nan)
  
  #setting color scale limits
  #wrapped in a function to apply for both limits_1 and limits_2
  limit_function=function(limit_no)
  {
    suffix <- limit_no
    
    maxlimit=max(abs(range(get(paste0('surf_data_',suffix)),na.rm = TRUE)))
    if(range(get(paste0('surf_data_',suffix)),na.rm = TRUE)[1]>=0) 
    {limits=reticulate::tuple(0,range(get(paste0('surf_data_',suffix)),na.rm = TRUE)[2])} ##if image contains all positive values
    else if(range(get(paste0('surf_data_',suffix)),na.rm = TRUE)[2]<=0) {limits=reticulate::tuple(range(get(paste0('surf_data_',suffix)),na.rm = TRUE)[1],0)} ##if image contains all negative values
    else {limits=reticulate::tuple(-maxlimit,maxlimit)} ##symmetrical limits will be used if image contains both positive and negative values 
    return(as.numeric(reticulate::py_to_r(limits))) #still need numeric for plotter
  }
  
  
  if(missing("limits_1")) {limits_1=limit_function(1)} 
  else #user specified limits
  {
    if(!is.null(limits_1)) 
    {limits=as.numeric(c(limits_1[1],limits_1[2]))} 
  }
  if(missing("limits_2")) {limits_2=limit_function(2)} 
  else #user specified limits
  {
    if(limits_2=='same'){limits_2=limits_1}
    else if(!is.null(limits_2)) 
    {limits=as.numeric(c(limits_2[1],limits_2[2]))} 
  }
  
  if(template %in% c('fsaverage5','fsaverage6','fslr32k'))
  {
    ##cortical surface plots
    #import python libraries
    brainstat.datasets=reticulate::import("brainstat.datasets", delay_load = TRUE)  
    brainspace.plotting=reticulate::import("brainspace.plotting", delay_load = TRUE)  
    #need PTuple function for overlays
    utils=reticulate::import("brainspace.plotting.utils", delay_load = TRUE) 
    
    #For brainstat data, it will look either in default $HOME path or 
    #custom if it's been set
    # If custom installation paths have been defined by the user, source
    # them from the package directory:
    Renvironpath=paste0(tools::R_user_dir(package='VertexWiseR'),'/.Renviron')
    if (file.exists(Renvironpath)) {readRenviron(Renvironpath)}
    
    if (Sys.getenv('BRAINSTAT_DATA')=="")
    { 
      brainstat_data_path=fs::path_home()
    } else if (!Sys.getenv('BRAINSTAT_DATA')=="") 
    {
      brainstat_data_path=Sys.getenv('BRAINSTAT_DATA')
    }
    
    #loading fsaverage surface
    left=brainstat.datasets$fetch_template_surface(template, join=FALSE, layer=surface,data_dir = paste0(brainstat_data_path,'/brainstat_data/surface_data/'))[1]
    right=brainstat.datasets$fetch_template_surface(template, join=FALSE, layer=surface,data_dir = paste0(brainstat_data_path,'/brainstat_data/surface_data/'))[2]
    
    #default cortical size and zoom parametes
    if(missing("size")) { size=c(1700,400)}
    if(missing("zoom")) { zoom=1.25 }

    #slice surface hemisphere-wise
    n_lh = left[[1]]$GetNumberOfPoints()
    n_rh = right[[1]]$GetNumberOfPoints()
    surf_data_1_lh = surf_data_1[1:n_lh]
    surf_data_1_rh = surf_data_1[(n_lh + 1):(n_lh + n_rh)]
    surf_data_2_lh = surf_data_2[1:n_lh]
    surf_data_2_rh = surf_data_2[(n_lh + 1):(n_lh + n_rh)]
    
    #Convert your overlay to VTK and name surface
    #and attach it to the corresponding surface template
    vtk_np = reticulate::import("vtk.util.numpy_support", convert = FALSE)
    vtk_array_lh_1 = vtk_np$numpy_to_vtk(reticulate::np_array(surf_data_1_lh), deep = TRUE)
    vtk_array_lh_1$SetName("lh_data_1")
    left[[1]]$GetPointData()$AddArray(vtk_array_lh_1)
    
    vtk_array_rh_1 = vtk_np$numpy_to_vtk(reticulate::np_array(surf_data_1_rh), deep = TRUE)
    vtk_array_rh_1$SetName("rh_data_1")
    right[[1]]$GetPointData()$AddArray(vtk_array_rh_1)
    
    vtk_array_lh_2 = vtk_np$numpy_to_vtk(reticulate::np_array(surf_data_2_lh, dtype = "float"), deep = TRUE)
    vtk_array_lh_2$SetName("lh_data_2")
    left[[1]]$GetPointData()$AddArray(vtk_array_lh_2)
  
    vtk_array_rh_2 = vtk_np$numpy_to_vtk(reticulate::np_array(surf_data_2_rh, dtype = "float"), deep = TRUE)
    vtk_array_rh_2$SetName("rh_data_2")
    right[[1]]$GetPointData()$AddArray(vtk_array_rh_2)
    
    surf_obj=reticulate::dict(lh = left[[1]],rh = right[[1]])
    
    ###
    #If NaN are plotted, add a coloured NaN layer under both layers if show_nan==T, because the default NaN has to be removed for the background to show behind the overlay's own NaN 
    if (show_nan==TRUE)
    { nan_mesh_lh=rep(0.43, length(surf_data_1_lh))
      nan_mesh_rh=rep(0.43, length(surf_data_1_rh))
    } else 
    { nan_mesh_lh=rep(NaN, length(surf_data_1_lh))
      nan_mesh_rh=rep(NaN, length(surf_data_1_rh))
    }
    
    # Create VTK arrays for NaN too
    vtk_array_lh_nan <- vtk_np$numpy_to_vtk(reticulate::np_array(nan_mesh_lh), deep = TRUE); 
    vtk_array_lh_nan$SetName("lh_data_nan")
    vtk_array_rh_nan <- vtk_np$numpy_to_vtk(reticulate::np_array(nan_mesh_rh), deep = TRUE); 
    vtk_array_rh_nan$SetName("rh_data_nan")
    # Add nan plot to surf_obj
    left[[1]]$GetPointData()$AddArray(vtk_array_lh_nan)
    right[[1]]$GetPointData()$AddArray(vtk_array_rh_nan)
    surf_obj=reticulate::dict(lh = left[[1]],rh = right[[1]])
    
    ###
    #If user wants boundaries plotted around overlay:
    mesh_ops <- reticulate::import("brainspace.mesh.array_operations", convert = FALSE)
    binary_mask_lh <- ifelse(is.na(surf_data_2_lh), 0L, 1L)
    binary_mask_rh <- ifelse(is.na(surf_data_2_rh), 0L, 1L)
    boundary_lh <- mesh_ops$get_labeling_border(left[[1]], reticulate::np_array(binary_mask_lh))$astype("float")
    boundary_rh <- mesh_ops$get_labeling_border(right[[1]], reticulate::np_array(binary_mask_rh))$astype("float")

    #NaN for non boundaries if wanted, and all NaN if to be hidden
    if (overlay_boundaries==TRUE)
    {
      #mask out non boundary vertices
      boundary_lh = reticulate::py_to_r(boundary_lh)
      boundary_rh = reticulate::py_to_r(boundary_rh) 
      boundary_lh[boundary_lh==0] = NaN
      boundary_rh[boundary_rh==0] = NaN
    } else {
      boundary_lh=rep(NaN, length(surf_data_1_rh))
      boundary_rh=rep(NaN, length(surf_data_1_rh))
    }

    
    # Create VTK arrays for boundaries too
    vtk_array_lh_bd <- vtk_np$numpy_to_vtk(reticulate::np_array(boundary_lh), deep = TRUE); 
    vtk_array_lh_bd$SetName("lh_data_bd")
    vtk_array_rh_bd <- vtk_np$numpy_to_vtk(reticulate::np_array(boundary_rh), deep = TRUE); 
    vtk_array_rh_bd$SetName("rh_data_bd")
    # Add nan plot to surf_obj
    left[[1]]$GetPointData()$AddArray(vtk_array_lh_bd)
    right[[1]]$GetPointData()$AddArray(vtk_array_rh_bd)
    surf_obj=reticulate::dict(lh = left[[1]],rh = right[[1]])
    
    #############################################################
    #Produce the plot 
    surf_plot = brainspace.plotting$plot_surf(
      surfs = surf_obj,
      array_name = list(list(
        utils$PTuple("lh_data_nan","lh_data_1", "lh_data_2","lh_data_bd"),
        utils$PTuple("lh_data_nan","lh_data_1", "lh_data_2","lh_data_bd"),
        utils$PTuple("rh_data_nan","rh_data_1", "rh_data_2","rh_data_bd"),
        utils$PTuple("rh_data_nan","rh_data_1", "rh_data_2","rh_data_bd")
      )
      ), #Referencing the attached scalar by name
      layout = list(list("lh","lh","rh","rh")),
      view = list(list("lateral", "medial","lateral","medial")),
      cmap = list(list(
        utils$PTuple('Greys',cmap_1,cmap_2,'Greys'), 
        utils$PTuple('Greys',cmap_1,cmap_2,'Greys'), 
        utils$PTuple('Greys',cmap_1,cmap_2,'Greys'), 
        utils$PTuple('Greys',cmap_1,cmap_2,'Greys')
      )),
      size = as.integer(size),
      nan_color = reticulate::tuple(0.7, 0.7, 0.7, 0),
      background = reticulate::tuple(as.integer(c(1, 1, 1))),
      zoom = zoom,
      color_range = list(list(
        utils$PTuple(c(0,1),limits_1,limits_2,c(0,1)), 
        utils$PTuple(c(0,1),limits_1,limits_2,c(0,1)), 
        utils$PTuple(c(0,1),limits_1,limits_2,c(0,1)), 
        utils$PTuple(c(0,1),limits_1,limits_2,c(0,1))
      )),
      label_text = title,
      return_plotter=TRUE,
      interactive = FALSE,
      transparent_bg = transparent_bg
    )
    
    #color bars need to be added on top of the original plot, since brainspace's plot_surf won't render them properly
    #plot_hemispheres does, so we can plot the bar(s) with it, take the bar(s), and append it to the main image
    if (colorbar_1==TRUE | colorbar_2==TRUE)
    {
      dummy_mesh=as.numeric(rep(NaN, length(surf_data_1)))
      
      background_plot=brainspace.plotting$plot_hemispheres(
        left[[1]], right[[1]],  
        array_name=reticulate::np_array(dummy_mesh),
        cmap=cmap_1,
        size=reticulate::tuple(as.integer(size)),
        nan_color=reticulate::tuple(0.7, 0.7, 0.7, 0), 
        background=reticulate::tuple(as.integer(c(1,1,1))),
        zoom=zoom,
        color_range=list(list(limits_1)),
        return_plotter = TRUE,
        interactive=FALSE, 
        color_bar=TRUE,  
        transparent_bg=transparent_bg)  
      bg_bar <- background_plot$screenshot(filename=paste0(tempdir(),"/background_cbar.png"),transparent_bg = transparent_bg)
        
        overlay_plot=brainspace.plotting$plot_hemispheres(
          left[[1]], right[[1]],  
          array_name=reticulate::np_array(dummy_mesh),
          cmap=cmap_2,
          size=reticulate::tuple(as.integer(size)),
          nan_color=reticulate::tuple(0.7, 0.7, 0.7, 0),
          background=reticulate::tuple(as.integer(c(1,1,1))),
          zoom=zoom,
          color_range=list(list(limits_2)),
          return_plotter = TRUE,
          interactive=FALSE, 
          color_bar=TRUE,  
          transparent_bg=transparent_bg)  
        ol_bar <- overlay_plot$screenshot(filename=paste0(tempdir(),"/overlay_cbar.png"), transparent_bg = transparent_bg)
        
        #Get colorbar position/size depending on plot size (approximate % of color bar depending on pic size)
        crop_scalar_bar <- function(img_path, start_frac = 0.90, end_frac = 0.947) {
          img <- png::readPNG(img_path)
          #get width based on proportion of image width
          start_col <- round(dim(img)[2] * start_frac)
          end_col   <- round(dim(img)[2] * end_frac)
          #crops based on dimensions 
          cropped_img <- img[, start_col:end_col, ]
          #overwrite
          png::writePNG(cropped_img, img_path)
        }
        crop_scalar_bar(paste0(tempdir(),"/background_cbar.png"))
        crop_scalar_bar(paste0(tempdir(),"/overlay_cbar.png"))
        
        #now append the cbar images
        # Load the images
        #Save with no color bar
        img_main=surf_plot$screenshot(filename=paste0(tempdir(),"/no_bar.png"),transparent_bg = transparent_bg)
        img_main <- png::readPNG(paste0(tempdir(),"/no_bar.png"))
        bg_bar  <- png::readPNG(paste0(tempdir(),"/background_cbar.png"))
        ol_bar  <- png::readPNG(paste0(tempdir(),"/overlay_cbar.png"))
        
        # Ensure same height and channels
        h <- dim(img_main)[1]; c <- dim(img_main)[3]
        w1 <- dim(img_main)[2]; w2 <- dim(bg_bar)[2]; w3 <- dim(ol_bar)[2]
        #expand dimensions depending on which color bars will be present
        w_total = w1
        if (colorbar_1==T) w_total <- w_total + w2
        if (colorbar_2==T) w_total <- w_total + w3
        # Create combined image
        combined_img <- array(0, dim = c(h, w_total, c))
        # Fill slices conditionally
        pos <- 1 #index of where color bar start in the plot
        combined_img[, pos:(pos + w1 - 1), ] <- img_main
        pos <- pos + w1
        if (colorbar_1==T) {
          combined_img[, pos:(pos + w2 - 1), ] <- bg_bar
          pos <- pos + w2}
        if (colorbar_2==T) {
          combined_img[, pos:(pos + w3 - 1), ] <- ol_bar}
        # Save
        png::writePNG(combined_img, filename) 
    } else
    {
      #Save with no color bar
      img_main=surf_plot$screenshot(filename=filename,transparent_bg = transparent_bg)
    }
    
    #display plot
    if(show.plot.window==TRUE)
    {
      img=png::readPNG(source=filename)
      grid::grid.raster(img)
    }
    
  } else { stop('This function only works with fsaverage5,fsaverage6, and fslr32k templates')}
} 



