plot_surf=function(surf_data, filename, title="",surface="inflated",cmap,limits, colorbar=TRUE, size, zoom, show.plot.window=TRUE,VWR_check=TRUE)
{
  #Check required python dependencies. If files missing:
  #Will prompt the user to get them in interactive session 
  #Will stop if it's a non-interactive session 
  if (VWR_check == TRUE){
    message("Checking for VertexWiseR system requirements ...")
    check = VWRfirstrun(n_vert=max(dim(t(surf_data))))
    if (!is.null(check)) {return(check)} else {message("\u2713 \n")}
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
  } else {rows=nrow(surf_data)}
  
  #in a multi-row data scenario: insert a dummy title if title is missing  or repeat the title nrow times
  if(rows>1) 
  {
    if(missing("title")) {title=rep(NULL,rows)}
    else if (missing("title")) {title=rep(title,rows)}
  }
  
  #check length of vector
  n_vert=length(surf_data)
  if(n_vert%%20484==0) {template="fsaverage5"}
  else if (n_vert%%64984==0) {template="fslr32k"} 
  else if (n_vert%%81924==0) {template="fsaverage6"} 
  else if (n_vert%%14524!=0) {stop("surf_data vector should only contain 20484 (fsaverage5), 81924 (fsaverage6), 64984 (fslr32k) or 14524 (hippocampal vertices) columns")}
  
  #if cmap is missing, select cmaps depending on whether the image contains positive only or negative only values
  if(missing("cmap"))
  {
    if(range(surf_data,na.rm = TRUE)[1]>=0)  {cmap="Reds"}
    else if (range(surf_data,na.rm = TRUE)[2]<=0)  {cmap="Blues_r"}
    else  {cmap="RdBu_r"}  
  }
  
  #custom cmap— if a vector of hex color codes is specified
  if(inherits(cmap,"colors")==TRUE)
  {
    matplotlib=reticulate::import("matplotlib")
    
    custom_colors=t(col2rgb(cmap)/255) # convert hex color codes to RGB codes, then divide by 255 to convert to RGBA codes
    
    #save as python cmap object
    mymap = matplotlib$colors$LinearSegmentedColormap$from_list('my_colormap', custom_colors)
    matplotlib$colormaps$unregister(name = "custom_cmap")
    matplotlib$colormaps$register(cmap = mymap,name="custom_cmap")
    cmap="custom_cmap"  
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
  
  if(n_vert%%14524!=0)
  {
    ##cortical surface fplots
    #import python libraries
    brainstat.datasets=reticulate::import("brainstat.datasets", delay_load = TRUE)  
    brainspace.plotting=reticulate::import("brainspace.plotting", delay_load = TRUE)  
    
    #loading fsaverage surface
    left=brainstat.datasets$fetch_template_surface(template, join=FALSE, layer=surface)[1]
    right=brainstat.datasets$fetch_template_surface(template, join=FALSE, layer=surface)[2]
    
    #default cortical size and zoom parametes
    if(missing("size")) { size=c(1920,rows*400)}
    if(missing("zoom")) { zoom=1.25 }
    
    surf_plot=brainspace.plotting$plot_hemispheres(left[[1]], right[[1]],  array_name=reticulate::np_array(surf_data),cmap=cmap, 
                                                   size=reticulate::tuple(as.integer(size)),nan_color=reticulate::tuple(0.7, 0.7, 0.7, 1),
                                                   return_plotter=TRUE,background=reticulate::tuple(as.integer(c(1,1,1))),zoom=zoom,color_range=limits,
                                                   label_text=title,interactive=FALSE, color_bar=colorbar,  transparent_bg=FALSE)  ##disabling interactive mode because this causes RStudio to hang
  } else
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
      surf_data.3d=array(NA,c(7262,2,nrow(surf_data))) #if N>1
      for (row in 1:nrow(surf_data))  {surf_data.3d[,,row]=cbind(surf_data[row,1:7262],surf_data[row,7263:14524])}
      surf_data=surf_data.3d
    }
    
    surf_plot=surfplot_canonical_foldunfold(surf_data,hipdat =get('hip_points_cells'),color_bar=colorbar,share="row",nan_color=reticulate::tuple(0.7, 0.7, 0.7, 1),size=as.integer(size), zoom=zoom,
                                            cmap=cmap,color_range=limits,label_text=title, return_plotter=TRUE,interactive=FALSE) ##disabling interactive mode because this causes RStudio to hang
  }
  #output plot as a .png image
  surf_plot$screenshot(filename=filename,transparent_bg = FALSE)
  if(show.plot.window==T)
  {
    img=png::readPNG(source =filename)
    grid::grid.raster(img)
  }
}
