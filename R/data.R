#' Atlas parcellations of fsaverage5
#'
#' A list containing two data frames, 1) listing fsaverage5 vertices and each parcellation number they correspond to, and 2) listing each available atlas and their corresponding labels (1=aparc, 2=Destrieux-148, 3=Glasser-360, 4=Schaefer-100, 5=Schaefer-200, 6=Schaefer-400). 
#'
#' @format ## `ROImap_fs5`
#' A list object with two data frame objects: () 
#' \describe{
#'   \item{vertices}{data frame with 20484 rows (vertices), 6 columns (atlases)}
#'   \item{atlases}{data frame with 400 rows (labels, not all are filled depending on atlas), 6 columns (atlases)}
#' }
"ROImap_fs5"

#' Atlas parcellations of fsaverage6
#'
#' A list containing two data frames, 1) listing fsaverage6 vertices and each atlas parcellation number they correspond to, and 2) listing each available atlas and their corresponding labels (1=aparc, 2=Destrieux-148, 3=Glasser-360, 4=Schaefer-100, 5=Schaefer-200, 6=Schaefer-400). 
#'
#' @format ## `ROImap_fs6`
#' A list object with two data frame objects: () 
#' \describe{
#'   \item{vertices}{data frame with 81924 rows (vertices), 6 columns (atlases)}
#'   \item{atlases}{data frame with 400 rows (labels, not all are filled depending on atlas), 6 columns (atlases)}
#' }
"ROImap_fs6"

#' Atlas parcellations of the hippocampus (CITI168)
#'
#' A list containing 1) a matrix listing CITI168 vertices and each parcellation number they correspond to, and  2) a data frame listing the  hippocampal atlas labels. 
#'
#' @format ## `ROImap_hip`
#' A list object with two data frame objects: () 
#' \describe{
#'   \item{vertices}{array of 14524 numeric vectors (vertices)}
#'   \item{atlases}{data frame with 10 rows listing names of left and right hippocampal subfields}
#' }
"ROImap_hip"

#' fsaverage5 surface in MNI space
#'
#' A matrix with 20484 columns corresponding to the fsaverage5 vertices and 3 rows corresponding to each vertex's X,Y,Z coordinates in MNI space
#'
#' @format ## `MNImap_fs5`
#' A 3x20494 matrix object
#' \describe{
#'   \item{coordinates}{20484 rows (vertices), 3 columns (X,Y,Z coordinates)}
#' }
"MNImap_fs5"


#' fsaverage6 surface in MNI space
#'
#' A matrix with 81924 columns corresponding to the fsaverage6 vertices and 3 rows corresponding to each vertex's X,Y,Z coordinates in MNI space
#'
#' @format ## `MNImap_fs6`
#' A 3x81924 matrix object
#' \describe{
#'   \item{coordinates}{81924 rows (vertices), 3 columns (X,Y,Z coordinates)}
#' }
"MNImap_fs6"


#' Hippocampal surface in MNI space
#'
#' A matrix with 14524 columns corresponding to the hippocampal template vertices and 3 rows corresponding to each vertex's X,Y,Z coordinates in MNI space
#'
#' @format ## `MNImap_hip`
#' A 3x14524 matrix object
#' \describe{
#'   \item{coordinates}{14524 rows (vertices), 3 columns (X,Y,Z coordinates)}
#' }
"MNImap_hip"

#' List of edges for the fsaverage5 template
#'
#' A Nx2 matrix object listing each vertex of the fsaverage5 template and the vertices adjacent to it (making an edge together).
#'
#' @format ## `edgelist_fs5`
#' \describe{
#'   \item{Nx2 matrix object}{Matrix with two columns and N rows corresponding to the unique edges in the fsaverage5 surface}
#' }
"edgelist_fs5"

#' List of edges for the fsaverage6 template
#'
#' A Nx2 matrix object listing each vertex of the fsaverage5 template and the vertices adjacent to it (making an edge together).
#'
#' @format ## `edgelist_fs6`
#' \describe{
#'   \item{Nx2 matrix object}{Matrix with two columns and N rows corresponding to the unique edges in the fsaverage6 surface}
#' }
"edgelist_fs6"

#' List of edges for the hippocampal template
#'
#' A Nx2 matrix object listing each vertex of the hippocampal template and the vertices adjacent to it (making an edge together).
#'
#' @format ## `edgelist_hip`
#' \describe{
#'   \item{Nx2 matrix object}{Matrix with two columns and N rows corresponding to the unique edges in the fsaverage5 surface}
#' }
"edgelist_hip"

#' points and cells data required to build the hippocampus surface template
#'
#' @format ## `hip_points_cells` 
#' A list object with two data frame objects: () 
#' \describe{
#'   \item{vertices}{data frame with 7262 rows (vertices), 3 columns (MNI coordinates X, y, Z)}
#'   \item{vertices}{data frame with 14266 rows (vertices), 3 columns (vertices of all unique triangles}
#'   \item{vertices}{data frame with 7262 rows (vertices), 3 columns (MNI coordinates X, y, Z for unfolded hippocampal surface)}
#' }
#' @docType data
"hip_points_cells"


#' fsaverage6 template object for nearest neighbor conversion in fs6_to_fs5()
#'
#' @format ## `fs6_to_fs5_map` 
#' An array of 81924 integers () 
#' \describe{
#'   \item{vertices}{81924 integers corresponding to each fsaverage6 vertex}
#' }
#' @docType data
"fs6_to_fs5_map"
