#' Region-of-Interest mapping object
#' 
#' A class for surface vertices mapping on atlas labels 
#' 
#' @slot matrix  A matrix object with N vertices from a template and each parcellation number the vertices correspond in 6 atlases (6 columns).
#' @slot atlases Each available of the 6 available atlases and their corresponding labels (1=aparc, 2=Destrieux-148, 3=Glasser-360, 4=Schaefer-100, 5=Schaefer-200, 6=Schaefer-400). 
#' @slot name The name of the template surface
#' 
#' @export

setClass("ROImap",
         slots = list(data = "matrix", 
                      atlases = 'data.frame',
                      name = 'character')
        )

ROI_map <- function(data,atlases) 
{
  if (nrow(data)==20484) {templatename='fsaverage5'}
  if (nrow(data)==81924) {templatename='fsaverage6'}
  if (nrow(data)==14524) {templatename='CITI168'}
  
  new("ROImap", data = data, atlases = atlases, name = templatename)
}

setMethod("show", "ROImap", function(object) {
  cat("[Region-of-Interest mapping object]\n")
  cat("Surface template:",object@name, "\n")
  
  if (length(colnames(object@atlases))>1) {
  cat("Matrix:",nrow(object@data), "vertices *", ncol(object@data), "matching label number in each atlas\n")
  cat("Atlases: "); cat(colnames(object@atlases))
  } else
  {cat("Matrix:",nrow(object@data), "vertices *", ncol(object@data), "matching label number in the atlas\n")
  cat("Atlas: "); cat(colnames(object@atlases))
  }
  
})

#' MNI surface map object
#' 
#' A class for surface coordinates in MNI space
#' 
#' @slot matrix  A matrix object with N vertices columns x  3 rows corresponding to each vertex's X,Y,Z coordinates in MNI space.
#' @slot name The name of the template surface
#' 
#' @export

setClass("MNIsurface",
         slots = list(data = "matrix", 
                      name = 'character')
         )

MNI_surface <- function(data) {
  if (ncol(data)==20484) {templatename='fsaverage5'}
  if (ncol(data)==81924) {templatename='fsaverage6'}
  if (ncol(data)==14524) {templatename='CITI168'}
  
  new("MNIsurface", data = data, name = templatename)
}

setMethod("show", "MNIsurface", function(object) {
  cat("[MNI surface map object]\n")
  cat("Surface template:",object@name,"\n") 
  cat(ncol(object@data), "vertices *", nrow(object@data), "X,Y,Z MNI coordinates\n")
})

#' Edge list object
#' 
#' A class for adjacent vertex correspondance 
#' 
#' @slot matrix   A N x 2 matrix object listing each vertex of the template and the vertices adjacent to it (making an edge together).
#' @slot name The name of the template surface
#' 
#' @export

setClass("edgelist",
slots = list(data = "matrix", 
             name = 'character')
        )

edge_list <- function(data) {
  if (nrow(data)==55876) {templatename='fsaverage5'}
  if (nrow(data)==224310) {templatename='fsaverage6'}
  if (nrow(data)==43054) {templatename='CITI168'}
  
  new("edgelist", data = data, name = templatename)
}

setMethod("show", "edgelist", function(object) {
  cat("[Edge listing object]\n")
  cat("Surface template:",object@name,"\n") 
  cat(ncol(object@data), "index numbers *", nrow(object@data), "pairs of adjacent vertices (edges)\n")
})
