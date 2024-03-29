## FUNCTION FOR VERTEX-WISE ANALYSIS WITH MIXED EFFECTS
## ADAPTED FROM brainstat python library (https://brainstat.readthedocs.io/en/master/_modules/brainstat/stats/SLM.html#SLM)
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

############################################################################################################################
############################################################################################################################
##vertex wise analysis with mixed effects
vertex_analysis=function(all_predictors,IV_of_interest, random, CT_data, p=0.05, atlas=1, smooth_FWHM)  ## atlas: 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360, 5=Destrieux-148
{
  if(class(IV_of_interest)=="integer") {IV_of_interest=as.numeric(IV_of_interest)}
  
  ##load other vertex-wise functions
  source("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/otherfunc.r?raw=TRUE")
  
  ##checks
  #check if required packages are installed
  packages="reticulate"
  new.packages = packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) 
  {
    cat(paste("The following package(s) are required and will be installed:\n",new.packages,"\n"))
    install.packages(new.packages)
  }  
  
  #check IV_of_interest
  if(NCOL(all_predictors)>1)
  {
    for(colno in 1:(NCOL(all_predictors)+1))
    {
      if(colno==(NCOL(all_predictors)+1))  {warning("IV_of_interest is not contained within all_predictors")}
      
      if(class(IV_of_interest)=="character") 
      {
        if(identical(IV_of_interest,data.matrix(all_predictors)[,colno]))  {break} 
      } else 
      {
        if(identical(suppressWarnings(as.numeric(IV_of_interest)),suppressWarnings(as.numeric(all_predictors[,colno]))))  {break}
      }
    }
  }  else
  {
    if(class(IV_of_interest)=="character") 
    {
      if(identical(IV_of_interest,all_predictors))  {colno=1} 
      else  {warning("IV_of_interest is not contained within all_predictors")}
    } else
    {
      if(identical(as.numeric(IV_of_interest),as.numeric(all_predictors)))  {colno=1}
      else  {warning("IV_of_interest is not contained within all_predictors")}
    }
  }
  
  #check if nrow is consistent for all_predictors and CT_data
  if(NROW(CT_data)!=NROW(all_predictors))  {stop(paste("The number of rows for CT_data (",NROW(CT_data),") and all_predictors (",NROW(all_predictors),") are not the same",sep=""))}
  
  #incomplete data check
  idxF=which(complete.cases(all_predictors)==F)
  if(length(idxF)>0)
  {
    cat(paste("all_predictors contains",length(idxF),"subjects with incomplete data. Subjects with incomplete data will be excluded in the current analysis\n"))
    all_predictors=all_predictors[-idxF,]
    IV_of_interest=IV_of_interest[-idxF]
    CT_data=CT_data[-idxF,]
  }
  
  #check categorical and recode variable
  if(NCOL(all_predictors)>1)
  {
    for (column in 1:NCOL(all_predictors))
    {
      if(class(all_predictors[,column])=="character") 
      {
        if(length(unique(all_predictors[,column]))==2)
        {
          cat(paste("The binary variable '",colnames(all_predictors)[column],"' will be recoded with ",unique(all_predictors[,column])[1],"=0 and ",unique(all_predictors[,column])[2],"=1 for the analysis\n",sep=""))
          
          recode=rep(0,NROW(all_predictors))
          recode[all_predictors[,column]==unique(all_predictors[,column])[2]]=1
          all_predictors[,column]=recode
          IV_of_interest=all_predictors[,colno]
        } else if(length(unique(all_predictors[,column]))>2)    {stop(paste("The categorical variable '",colnames(all_predictors)[column],"' contains more than 2 levels, please code it into binarized dummy variables",sep=""))}
      }      
    }
  } else
  {
    if(class(all_predictors)=="character") 
    {
      if(length(unique(all_predictors))==2)
      {
        cat(paste("The binary variable '",colnames(all_predictors),"' will be recoded such that ",unique(all_predictors)[1],"=0 and ",unique(all_predictors)[2],"=1 for the analysis\n",sep=""))
        
        recode=rep(0,NROW(all_predictors))
        recode[all_predictors==unique(all_predictors)[2]]=1
        all_predictors=recode
        IV_of_interest=all_predictors
      } else if(length(unique(all_predictors))>2)    {stop(paste("The categorical variable '",colnames(all_predictors),"' contains more than 2 levels, please code it into binarized dummy variables",sep=""))}
    }      
  }
  
  #check length of CT data and load the appropriate fsaverage files
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
  else if (n_vert==14524)
  {
    if(file.exists("hip_template.fs")==F)
    {
      cat("\nhip_template.fs is not detected in the current working directory. The hippocampus surface template will be downloaded\n")
      download.file(url="https://raw.githubusercontent.com/CogBrainHealthLab/VertexWiseR/main/data/hip_template.fs",destfile ="hip_template.fs",mode = "wb")
    } 
    brainspace.mesh.mesh_io=reticulate::import("brainspace.mesh.mesh_io")
    template=brainspace.mesh.mesh_io$read_surface("hip_template.fs")
    load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap_hip.rdata?raw=TRUE"))
  }
  else {stop("data vector should only contain 20484 (fsaverage5), 81924 (fsaverage6) or 14524 (hippocampal vertices) columns")}
  
  ##smoothing
  n_vert=NCOL(CT_data)
  if(missing("smooth_FWHM"))
  {
    if(n_vert==20484) 
    {
      cat("CT_data will be smoothed using the default 10mm FWHM kernel for fsaverage5 images\n")
      CT_data=smooth(CT_data, FWHM=10)
    } else if(n_vert==81924) 
    {
      cat("CT_data will be smoothed using the default 5mm FWHM kernel for fsaverage6 images")
      CT_data=smooth(CT_data, FWHM=5)
    } else if(n_vert==14524) 
    {
      cat("CT_data will be smoothed using the default 5mm FWHM kernel for hippocampal maps\n")
      CT_data=smooth(CT_data, FWHM=5)
    }
  } else if(smooth_FWHM>0) 
  {
    cat(paste("CT_data will be smoothed using a ",smooth_FWHM,"mm FWHM kernel", sep=""))
    CT_data=smooth(CT_data, FWHM=smooth_FWHM)
  }
  CT_data[is.na(CT_data)]=0
  
  ##import python libaries
  brainstat.stats.terms=reticulate::import("brainstat.stats.terms")
  brainstat.stats.SLM=reticulate::import("brainstat.stats.SLM")
  brainstat.datasets=reticulate::import("brainstat.datasets")  
  
  ##fitting model
  #preparing mask for model
  mask=array(rep(T,NCOL(CT_data)))
  maskNA=which(colSums(CT_data != 0) == 0)
  mask[which(colSums(CT_data != 0) == 0)]=F
  
  #fit model
  if(missing("random")) {model0=brainstat.stats.terms$FixedEffect(all_predictors, "_check_categorical" = F)}
  else {model0=brainstat.stats.terms$MixedEffect(ran = random,fix = all_predictors,"_check_categorical" = F)}
  model=brainstat.stats.SLM$SLM(model = model0,
                                contrast=IV_of_interest,
                                surf = template,
                                mask=mask,
                                correction=c("fdr", "rft"),
                                cluster_threshold=p)
  model$fit(CT_data)
  
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
    pos_clusterIDmap=rep(0, NCOL(CT_data))
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
    neg_clusterIDmap=rep(0, NCOL(CT_data))
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
  
  ##generating thresholded t-stat vector (for plotting)
  tstat[intersect(which(neg_clusterIDmap==0),which(pos_clusterIDmap==0))]=NA
  tstat[is.na(tstat)]=0
  tstat[maskNA]=NA
  
  ##generating positive and negative masks
  posmask=array(rep(0,NCOL(CT_data)))
  posmask[which(tstat>0)]=1
  
  negmask=array(rep(0,NCOL(CT_data)))
  negmask[which(tstat<0)]=1
  
  #listing objects to return
  returnobj=list(cluster_results,tstat,posmask,negmask,pos_clusterIDmap,neg_clusterIDmap)
  names(returnobj)=c("cluster_level_results","thresholded_tstat_map","pos_mask","neg_mask","pos_clusterIDmap","neg_clusterIDmap")
  return(returnobj)
}
############################################################################################################################
############################################################################################################################

##example
# source("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/vertwise.r?raw=TRUE")
# 
# model=vertex_analysis(fixed_predictors = dat_beh[,c(2:4)], IV_of_interest = dat_beh$AGE_AT_SCAN, random = dat_beh$SUB_ID,CT_data = dat_CT,p = 0.01, atlas=1)
# model$cluster_level_results
