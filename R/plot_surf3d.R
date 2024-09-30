#' @title 3D Surface plotter
#'
#' @description Plots 3D cortical surfaces
#'
#' @param surf_data  A numeric vector (length of V) 
#' @param surf_color  color of the cortical surface. Set to `'grey'` by default
#' @param cmap A string vector containing 2 to 4 color names/codes specifying the colors to be used for the color scale. If none are specified, appropriate colors will be automatically selected according to `range(surf_data)`
#' @param limits A combined pair of numeric vector composed of the lower and upper color scale limits of the plot. When left unspecified, the symmetrical limits `c(-max(abs(surf_dat),max(abs(surf_dat)))` will be used. 
#' @param atlas atlas used for identifying region labels. 1=Desikan, 2=Destrieux-148, 3=Glasser-360, 4=Schaefer-100, 5=Schaefer-200, 6=Schaefer-400. Set to `1` by default. This argument is ignored for hippocampal surfaces.
#' @param hemi A string specifying the hemisphere to plot. Possible values are `l` (left), `r` (right) or `b` (both).
#' @param orientation_labels A boolean object specifying if orientation labels are to be displayed. Set to `TRUE` by default
#' @param VWR_check A boolean object specifying whether to check and validate system requirements. Default is TRUE.
#'
#' @returns a plot_ly object
#' @examples
#' surf_data = runif(20484);
#' plot_surf3d(surf_data = surf_data)
#' @importFrom plotly plot_ly add_trace layout
#' @export
######################################################################################################################################################
######################################################################################################################################################
plot_surf3d=function(surf_data, surf_color="grey",cmap,limits, atlas=1, hemi="b",orientation_labels=TRUE,VWR_check=TRUE)
{
  #Check required python dependencies. If files missing:
  #Will prompt the user to get them in interactive session 
  #Will stop if it's a non-interactive session 
    if (VWR_check == TRUE){
      message("Checking for VertexWiseR system requirements ...")
      check = VWRfirstrun(n_vert=max(dim(t(surf_data))))
      if (!is.null(check)) {return(check)}
    } else if(interactive()==FALSE) { return(message('Non-interactive sessions need requirement checks'))}
    
  ##check selected hemispheres and write orientation labels
    if(hemi=="l") {xlab=c("Left", "Medial")}
    else if(hemi=="r") {xlab=c("Medial", "Right")}
    else if(hemi=="b") {xlab=c("Left", "Right")}
    else {stop("hemi can only take values of 'l' (left), 'r' (right) or 'b' (both)")}
    
  ##color scale
    if(missing("cmap"))
    {
      if(range(surf_data,na.rm = TRUE)[1]>=0)  {cmap=c("#A51122","#F5FACD")}
      else if (range(surf_data,na.rm = TRUE)[2]<=0)  {cmap=c("#324DA0","#E7F1D5")}
      else  {cmap=c("#E7F1D5","#324DA0","#A51122","#F5FACD")}  
    }
    
    # enabling custom color scales
    if(length(cmap)==2) {cmap=list(list(0,cmap[1]), list(1,cmap[2]))} 
    else if(length(cmap)==3){cmap=list(list(0,cmap[1]), list(0.5,cmap[2]),list(1,cmap[3]))}
    else if(length(cmap)==4){cmap=list(list(0,cmap[1]), list(0.5,cmap[2]),list(0.51,cmap[3]),list(1,cmap[4]))}
    
  
  ## selecting template and ROI map depending on no. vertices
    n_vert=length(surf_data)
    if(n_vert==20484)
    {
      tri=t(get_faces("fsaverage5"))
      coords=t(get_MNIcoords("fsaverage5"))
      ROImap_fs5 <- get('ROImap_fs5')
      ROImap <- list(ROImap_fs5@data,ROImap_fs5@atlases)
    } else if (n_vert==64984)
    {
      tri=t(VertexWiseRget_faces("fslr32k"))
      coords=t(get_MNIcoords("fslr32k"))
      ROImap_fslr32k <- get('ROImap_fslr32k')
      ROImap <- list(ROImap_fslr32k@data,ROImap_fslr32k@atlases)
    } else if (n_vert==81924)
    {
      tri=t(get_faces("fsaverage6"))
      coords=t(get_MNIcoords("fsaverage6"))
      ROImap_fs6 <- get('ROImap_fs6')
      ROImap <- list(ROImap_fs6@data,ROImap_fs6@atlases)
    } else if (n_vert==14524)
    {
      LH.hip.mni=get("hip_points_cells")[[1]]
      RH.hip.mni=LH.hip.mni
      RH.hip.mni[,1]=LH.hip.mni[,1]*-1
      coords=rbind(LH.hip.mni,RH.hip.mni)
      LH.tri=get("hip_points_cells")[[2]]+1
      tri=rbind(LH.tri,LH.tri+7262)
      
      ROImap_hip <- get('ROImap_hip')
      ROImap <- list(ROImap_hip@data,ROImap_hip@atlases)
      atlas=1
    } else 
    {stop("data vector should only contain 20484 (fsaverage5), 81924 (fsaverage6), 64984 (fslr32k) or 14524 (hippocampal vertices) columns")}

  ##Averaging vertices values to obtain face values
    face.stat=rowMeans(cbind(surf_data[tri[,1]],surf_data[tri[,2]],surf_data[tri[,3]]))
    non0.idx=which(face.stat>0)
  
  ##setting color scale limits
    maxlimit=max(abs(range(face.stat,na.rm = TRUE)))
    if(missing(limits)) 
    {
      limits.range=range(face.stat,na.rm = T)
      if(limits.range[1]>=0) {limits=c(0,limits.range[2])} ##if image contains all positive values
      else if(limits.range[2]<=0) {limits=c(limits.range[1],0)} ##if image contains all negative values
      else if(limits.range[1]<0 & limits.range[2]>0){limits=c(-maxlimit,maxlimit)} ##symmetrical limits will be used if image contains both positive and negative values
    } else {limits=c(limits[1],limits[2])    }
    
  ##splitting cortical data in to LH and RH if necessary 
    mid.idx=n_vert/2
    mid.tri=NROW(tri)/2
    
    if(hemi=="l")
    {
      surf_data=surf_data[1:mid.idx]
      coords=coords[1:mid.idx,]
      tri=tri[1:mid.tri,]
      face.stat=face.stat[1:mid.tri]
      face.stat.non0.idx=which(abs(face.stat)>0)
      
      ROI.idx=ROImap[[1]][,atlas][1:mid.idx]
      ROI.idx[ROI.idx==0]=max(ROI.idx)+1
      ROI.text=ROImap[[2]][,atlas][ROI.idx]
      ROI.text=ROI.text[1:mid.idx]
      
    } else if(hemi=="r")
    {
      surf_data=surf_data[(mid.idx+1):n_vert]
      coords=coords[(mid.idx+1):n_vert,]
      tri=tri[(mid.tri+1):NROW(tri),]-mid.idx
      face.stat=face.stat[(mid.tri+1):length(tri)]
      face.stat.non0.idx=which(abs(face.stat)>0)
      
      
      ROI.idx=ROImap[[1]][,atlas]
      ROI.idx[ROI.idx==0]=max(ROI.idx)+1
      ROI.text=ROImap[[2]][,atlas][ROI.idx]
      ROI.text=ROI.text[(mid.idx+1):n_vert]
    } else if(hemi=="b")
    {
      face.stat.non0.idx=which(abs(face.stat)>0)
      ROI.idx=ROImap[[1]][,atlas]
      ROI.idx[ROI.idx==0]=max(ROI.idx)+1
      ROI.text=ROImap[[2]][,atlas][ROI.idx]
    }
  
  ##create blank cortical surface
    fig=plotly::plot_ly(type = 'mesh3d',
                        x = coords[,1],
                        y = coords[,2],
                        z = coords[,3],
                        i = tri[, 1] - 1,  # plotly uses 0-based indexing, so subtract 1
                        j = tri[, 2] - 1,
                        k = tri[, 3] - 1,facecolor=surf_color)
  
  ##overlay statistical map on cortical surface
    fig=add_trace(fig,type = 'mesh3d',
                  x = coords[,1],
                  y = coords[,2],
                  z = coords[,3],
                  i = tri[face.stat.non0.idx,][, 1] - 1,  # plotly uses 0-based indexing, so subtract 1
                  j = tri[face.stat.non0.idx,][, 2] - 1,
                  k = tri[face.stat.non0.idx,][, 3] - 1,
                  intensitymode="cell",intensity=face.stat[face.stat.non0.idx],
                  colorscale = cmap,
                  cauto = F,
                  cmin = limits[1],
                  cmax = limits[2])

  ##add mouse-over text                 
    fig=add_trace(fig,text=ROI.text,intensitymode="vertex", intensity=0, opacity=0,showscale= F,
                  x = coords[,1],
                  y = coords[,2],
                  z = coords[,3],
                  i = tri[, 1] - 1,  # plotly uses 0-based indexing, so subtract 1
                  j = tri[, 2] - 1,
                  k = tri[, 3] - 1,
                  customdata=surf_data,
                  hovertemplate=paste("Region label: %{text}<br>",
                                      "MNI coords: %{x:.1f},%{y:.1f},%{z:.1f}<br>",
                                      "statistic:%{customdata:.2f}<extra></extra>"))
    
  ##axis parameters
  p=plotly::layout(fig,
                   hoverlabel = list(align = "left"),
                   scene = list(camera=list(eye = list(x = 0, y = 1.5, z = 1.5)),
                                xaxis = list(showgrid = T,showticklabels=T,showspikes=F,zeroline=F, title=""),
                                yaxis = list(showgrid = T,showticklabels=T,showspikes=F,zeroline=F, title=""),
                                zaxis = list(showgrid = T,showticklabels=T,showspikes=F,zeroline=F, title="")))
  
  ##add optional orientation labels
  if(orientation_labels==T)
  {
    axx = list(ticketmode = 'array',ticktext = xlab,tickvals = range(coords[,1]))
    axy = list(ticketmode = 'array',ticktext = c("Posterior","Anterior"),tickvals = range(coords[,2]))
    axz = list(ticketmode = 'array',ticktext = c("Inferior","Superior"),tickvals = range(coords[,3]))
    
    p = layout(p,scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  }
  return(p)
}

############################################################################################################################
############################################################################################################################
##getting the triangle faces of brainstat templates
get_faces=function(template)
{
  #Load brainstat tools
  brainstat.datasets=reticulate::import("brainstat.datasets")  
  brainspace.mesh.mesh_elements=reticulate::import("brainspace.mesh.mesh_elements")
  
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
  
  #Loads template surfaces
  surf.template=brainstat.datasets$fetch_template_surface(template=template, join=TRUE, data_dir=data_dir)
  
  #Returns MNI coordinates
  return(t(brainspace.mesh.mesh_elements$get_cells(surf.template))+1) #python indices start from 0, hence+1
}
