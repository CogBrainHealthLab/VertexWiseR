## OTHER VERTEX-WISE FUNCTIONS
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY
############################################################################################################################
############################################################################################################################
## various checks
check.inputs=function(packages,CT_data, all_predictors, IV_of_interest)
{
  #check if required packages are installed
  if(!missing(packages))
    {
      packages=as.character(packages)
      new.packages = packages[!(packages %in% installed.packages()[,"Package"])]
      if(length(new.packages)) 
      {
        cat(paste("The following package(s) are required and will be installed:\n",new.packages,"\n"))
        install.packages(new.packages)
      }  
    }
  #check if nrow is consistent for all_predictors and CT_data
  if(NROW(CT_data)!=NROW(all_predictors))  {stop(paste("The number of rows for CT_data (",NROW(CT_data),") and all_predictors (",NROW(all_predictors),") are not the same",sep=""))}
  
  #incomplete data check
  idxF=which(complete.cases(all_predictors)==F)
  if(length(idxF)>0)
  {
    cat(paste("all_predictors contains",length(idxF),"subjects with incomplete data. Subjects with incomplete data will be excluded in the current analysis"))
    all_predictors=all_predictors[-idxF,]
    IV_of_interest=IV_of_interest[-idxF]
    CT_data=CT_data[-idxF,]
  }
  
  #check IV_of_interest
  for(colno in 1:(NCOL(all_predictors)+1))
  {
    if(colno==(NCOL(all_predictors)+1))  {stop("IV_of_interest is not contained within all_predictors")}
    
    if(class(IV_of_interest) != "integer" & class(IV_of_interest) != "numeric") 
    {
      if(identical(IV_of_interest,all_predictors[,colno]))  {break} 
    } else 
    {
      if(identical(as.numeric(IV_of_interest),as.numeric(all_predictors[,colno])))  {break}
    }
  }

  #check categorical variable
  for (column in 1:NCOL(all_predictors))
  {
    if(class(all_predictors[,column]) != "integer" & class(all_predictors[,column]) != "numeric")
    {
      if(length(unique(all_predictors[,column]))==2)
      {
        cat(paste("The binary variable '",colnames(all_predictors)[column],"' will be recoded with ",unique(all_predictors[,column])[1],"=0 and ",unique(all_predictors[,column])[2],"=1 for the analysis\n",sep=""))
        
        recode=rep(0,NROW(all_predictors))
        recode[all_predictors[,column]==unique(all_predictors[,column])[2]]=1
        all_predictors[,column]=recode
        IV_of_interest=all_predictors[,colno]
      } else if(length(unique(all_predictors[,column]))>2)    {cat(paste("The categorical variable '",colnames(all_predictors)[column],"' contains more than 2 levels, please code it into binarized dummy variables",sep=""))}
    }      
  }
  
  #check length of CT data and smoothing
  n_vert=ncol(CT_data)
  if(n_vert==20484)
  {
    template="fsaverage5"
    load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap_fs5.rdata?raw=TRUE"))
  }
  else if (n_vert==81924)
  {
    template="fsaverage6"
    load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap_fs6.rdata?raw=TRUE"))
  }
  else {stop("CT_data should only contain 20484 (fsaverage5) or 81924 (fsaverage6) columns")}
  
  #check for collinearity
  cormat=cor(all_predictors,use = "pairwise.complete.obs")
  cormat.0=cormat
  cormat.0[cormat.0==1]=NA
  if(max(abs(cormat.0),na.rm = T) >0.5)
  {
    warning(paste("correlations among variables in all_predictors are observed to be as high as ",round(max(abs(cormat.0),na.rm = T),2),", suggesting potential collinearity among predictors.\nAnalysis will continue...",sep=""))
  }
  #save objects to global environment
  listobj=list(template,ROImap)
  names(listobj)=c("template","ROImap")
  list2env(listobj,envir=globalenv())
}

############################################################################################################################
############################################################################################################################
## Efficient way to extract t statistics from linear regression models to speed up the permutation process
## adapted from https://stackoverflow.com/questions/15820623/obtain-t-statistic-for-regression-coefficients-of-an-mlm-object-returned-by-l
extract.t=function(mod,row)
{
  p = mod$rank
  rdf = mod$df.residual
  Qr = mod$qr
  p1 = 1L:p
  r = mod$residuals
  rss = colSums(r^2)
  resvar = rss/rdf
  R = chol2inv(Qr$qr[p1, p1, drop = FALSE])  
  se = (sqrt(diag(R) %*% t(resvar)))[row,]
  est = mod$coefficients[row,]
  tval = est/se                          
}
############################################################################################################################
############################################################################################################################
##find clusters using edgelist
getClusters=function(data)
{ 
  n_vert=length(data)
  
  #listing out non-zero vertices
  vert=which(data!=0)

  #matching non-zero vertices with adjacency matrices to obtain list of edges connecting between the non-zero vertices
  edgelist0=edgelist[!is.na(match(edgelist[,1],vert)),]
  if(length(edgelist0)>2)  {edgelist1=edgelist0[!is.na(match(edgelist0[,2],vert)),]} 
  else if (length(edgelist0)==2)  ##if only a single edge was identified, edgelist will no longer be a Nx2 matrix, hence need to reshape it into a matrix
    { 
    edgelist0=matrix(edgelist0,ncol=2,nrow=1)
    edgelist1=edgelist0[!is.na(match(edgelist0[,2],vert)),]
    } else {edgelist1=0}
  remove(data,vert,edgelist0)
  
  if(length(edgelist1)>2) #if at least 2 edges are identified
  {
    #extracting cluster-related info from list of non-zero edges
    com=igraph::components(igraph::graph.data.frame(edgelist1, directed = F))
    clust.size=com$csize
    
    #cluster mappings
    clust.map=rep(NA,n_vert)
    clust.map[as.numeric(names(com$membership))]=com$membership
  
  } else if(length(edgelist1)==2) #bypass cluster extraction procedure if only 1 edge is identified
  {
    clust.size=2
    clust.map=rep(NA,n_vert)
    clust.map[edgelist1]=1
  } else #bypass cluster extraction procedure if no edges are identified
  {
    clust.map="noclusters"
    clust.size="noclusters"
  }
  return(list(clust.map,clust.size))
}
############################################################################################################################
############################################################################################################################
## To extract atlas ROI values from fsaverage5 vertex-wise data and vice-versa
fs5_to_atlas=function(data,atlas) ## atlas: 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360, 5=Destrieux-148
{  
  #check length of vector
  if(length(data)%%20484!=0) {stop("Length of data is not a multiple of 20484")}
  
  #load atlas mapping data
  load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap.rdata?raw=TRUE"))
  
  #init variables
  nregions=max(ROImap[[1]][,atlas])
  data[is.na(data)]=0

  #mapping fsaverage5 space vertice to atlas regions if data is a 1x20484 vector
  if(length(data)==20484) 
  {
    data=matrix(data,ncol=20484,nrow=1)  
    ROI=rep(NA,nregions)
    for (region in 1:nregions)  {ROI[region]=mean(data[which(ROImap[[1]][,atlas]==region)])} #vertices are averaged within the atlas ROI
  } else 
  {
  #mapping fsaverage5 space vertice to atlas regions if data is a Nx20484 matrix
    ROI=matrix(NA, nrow=NROW(data), ncol=nregions)
    for (region in 1:nregions)  {ROI[,region]=rowMeans(data[,which(ROImap[[1]][,atlas]==region)])} #vertices are averaged within the atlas ROI
  }
  return(ROI)
}

atlas_to_fs5=function(data,atlas) ## atlas: 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360, 5=Destrieux-148
  {
    #load atlas mapping data
    load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap.rdata?raw=TRUE"))
  
    #init variables
    nregions=max(ROImap[[1]][,atlas])
    fs5_dat=rep(NA,20484)
  
    #mapping atlas label to fsaverag5 space
    for (region in 1:nregions)  {fs5_dat[which(ROImap[[1]][,atlas]==region)]=data[region]}
    return(as.numeric(fs5_dat))
  }
############################################################################################################################
############################################################################################################################
#convert between fsaverage5 and fsaverage6 spacing
fs5_to_fs6=function(data)
{
  #check length of vector
  if(length(data)%%20484!=0) {stop("Length of data is not a multiple of 20484")}
  
  #load atlas mapping data
  load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/fs6_to_fs5.rdata?raw=TRUE"))
  #mapping fsaverage5 to fsaverage6 space if data is a Nx20484 matrix
  if(length(data)==20484) {data.fs6=data[fs6_to_fs5]} 
  #mapping fsaverage5 to fsaverage6 space if data is a Nx20484 matrix
  else {data.fs6=data[,fs6_to_fs5]}
  return(data.fs6)
}

fs6_to_fs5=function(data)
{
  #check length of vector
  if(length(data)%%81924!=0) {stop("Length of data is not a multiple of 81924")}
  
  #load atlas mapping data
  load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/fs6_to_fs5.rdata?raw=TRUE"))
  
  if(length(data)==81924) #mapping fsaverage6 to fsaverage5 space if data is a Nx20484 matrix
  {
    data=matrix(data,ncol=81924,nrow=1)  
    data.fs5=matrix(NA,ncol=20484,nrow=1)
    
    for (vert in 1:20484)  {data.fs5[vert]=mean(data[fs6_to_fs5==vert],na.rm = T)} 
  } else #mapping fsaverage6 to fsaverage5 space if data is a Nx20484 matrix
  {
    data.fs5=matrix(NA,ncol=20484,nrow=NROW(data))
    for (vert in 1:20484)  {data.fs5[,vert]=rowMeans(data[,fs6_to_fs5==vert],na.rm = T)} 
  }
  return(data.fs5)
}
############################################################################################################################
############################################################################################################################
##smoothing fsaverage5 and fsaverage6 data
smooth=function(data, FWHM)
{
  ##import python libraries
  reticulate::source_python("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/smooth.py?raw=TRUE")

  ##fsaverage space specific parameters
  if(ncol(data)==20484) ##fsaverage5 parameters
    {
      load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/edgelistfs5.rdata?raw=TRUE"))
      vert_mm=3.5
    } else if (ncol(data)==81924) ##fsaverage6 parameters
    {
      load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/edgelistfs6.rdata?raw=TRUE"))
      vert_mm=2
    } 
    else {stop("data vector should only contain 20484 (fsaverage5) or 81924 (fsaverage6) columns")}

  ##smoothing
  smooth=mesh_smooth(Y=data, edg=edgelist, FWHM = FWHM/vert_mm)
  smooth[is.na(smooth)]=0
  return(smooth)  
}
############################################################################################################################
############################################################################################################################
##CT surface plots
plotCT=function(data, filename,title="",surface="inflated",cmap,fs_path, range=NULL , colorbar=T)
{
  #check length of vector
  n_vert=length(data)
  if(n_vert==20484) {template="fsaverage5"}
  else if (n_vert==81924) {template="fsaverage6"} 
  else {stop("data vector should only contain 20484 (fsaverage5) or 81924 (fsaverage6) columns")}

  #legacy input message
  if(!missing("fs_path")){cat("The fs_path parameter and the fsaverage5 files are no longer needed in the updated plotCT function\n")}

  #import python libraries
  brainstat.datasets=reticulate::import("brainstat.datasets")  
  brainspace.plotting=reticulate::import("brainspace.plotting")  

  #loading fsaverage surface
  left=brainstat.datasets$fetch_template_surface(template, join=F, layer=surface)[1]
  right=brainstat.datasets$fetch_template_surface(template, join=F, layer=surface)[2]

  #setting color maps
  if(missing("cmap"))
  {
    if(range(data,na.rm = T)[1]>=0)
      {
      cmap="Reds"
      range=NULL
      }
    else if (range(data,na.rm = T)[2]<=0)
      {cmap="Blues_r"
      range=NULL
      }
    else
      {
      cmap="RdBu_r"
      range="sym"
      }  
  }

  #plot object
  CTplot=brainspace.plotting$plot_hemispheres(left[[1]], right[[1]],  array_name=reticulate::np_array(data),cmap=cmap, 
                                              size=reticulate::tuple(as.integer(c(1920,400))),nan_color=reticulate::tuple(c(0.7, 0.7, 0.7, 1)),
                                              return_plotter=T,background=reticulate::tuple(as.integer(c(1,1,1))),zoom=1.25,color_range=range,
                                              label_text=list('left'=list(title)),interactive=F, color_bar=colorbar,  transparent_bg=FALSE)
  #output plot as a .png image
  CTplot$screenshot(filename=filename,transparent_bg = F)
}
############################################################################################################################
############################################################################################################################
##CT image decoding
decode_img=function(img,contrast="positive")
{
  ##checks
    #check length of vector
    n_vert=length(img)
    if(n_vert==20484) {template="fsaverage5"}
    else if (n_vert==81924) {stop("decoding of fsaverage6-space image is current not implemented, please resample the image to fsaverage5 space")} 
    else {stop("img vector should only contain 20484 (fsaverage5)")}

    #check contrast
    if(contrast != "positive" & contrast != "negative")  {stop("contrast has to be either positive or negative")} 
  
  ##import python libraries
  interpolate=reticulate::import("brainstat.mesh.interpolate")
  discrete=reticulate::import("nimare.decode")
  nimare.dataset=reticulate::import("nimare.dataset")
  
  ##selecting contrasts
  if(contrast=="positive")
  {
    img[is.na(img)]=0
    img[img<0]=0
    img[img>0]=1
  } else if (contrast=="negative")
  {
    img[is.na(img)]=0
    img[img>0]=0
    img[img<0]=1
  }

  ##convert img vector to nii image
  stat_labels=reticulate::r_to_py(img)
  stat_nii = interpolate$`_surf2vol`(template, stat_labels)
  
  ##download neurosynth database if necessary 
  if(file.exists("neurosynth_dataset.pkl.gz")==F)
  {
    cat("\nneurosynth_dataset.pkl.gz is not detected in the current working directory. The neurosynth database will be downloaded\n")
    download.file(url="https://raw.githubusercontent.com/CogBrainHealthLab/VertexWiseR/main/data/neurosynth_dataset.pkl.gz",destfile = "neurosynth_dataset.pkl.gz")
  } 
  ##running the decoding procedure
  neurosynth_dset = nimare.dataset$Dataset$load("neurosynth_dataset.pkl.gz")
  cat("Correlating input image with images in the neurosynth database. This may take a while\n")
  decoder = discrete$ROIAssociationDecoder(stat_nii)
  decoder$fit(neurosynth_dset)

  ##compiling the results
  decoder_df = data.matrix(decoder$transform())
  row.names(decoder_df)=gsub(pattern = "terms_abstract_tfidf__",x=row.names(decoder_df), replacement = "")
  result=data.frame(row.names(decoder_df),round(as.numeric(decoder_df),3))
  colnames(result)=c("keyword","r")
  result=result[order(-result$r),]
  
  return(result)
}  
############################################################################################################################
############################################################################################################################
