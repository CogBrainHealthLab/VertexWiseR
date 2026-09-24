#' @title Decode surface data
#'
#' @description Correlates the significant clusters of an earlier vertex-wise analysis with a database of task-based fMRI and voxel-based morphometric statistical maps and associate them with relevant key words. Decoding currently works with surfaces in fsaverage5 space only."
#'
#' @details An internal image decoding is used, reproducing the \href{https://nimare.readthedocs.io/en/stable/index.html}{'NiMARE'} python module's ROIAssociationDecoder in base R. The function also downloads the \href{https://github.com/neurosynth/neurosynth-data}{'Neurosynth' database} in the package's inst/extdata directory (converted to .RDS, ~7.2 MB) for the analysis.
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
  if(file.exists(system.file('extdata','neurosynth_dataset.rds', package='VertexWiseR'))==TRUE)
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
    
    ## Load the internal Neurosynth dataset
    # This dataset was converted form a NiMARE pickle to a .rds file so R can handle it without python. This made use of reticulate"s importer to load the dataset with nimare$Dataset$load, and the py_to_r() function to extract coordinates, annotations and sources as R arrays.
    neurosynth_dset <- readRDS(system.file("extdata", "neurosynth_dataset.rds", package = "VertexWiseR"))
    
    ##running the decoding procedure
    message("\u2713 \nCorrelating input image with images in the neurosynth database. This may take a while ... ")
    decoder <- ROIAssociationDecoder(nibabel_roi_to_R(stat_nii))
    decoder$fit(neurosynth_dset)
    
    ##compiling the results
    decoder_df = data.matrix(decoder$transform())
    row.names(decoder_df)=gsub(pattern = "terms_abstract_tfidf__",x=row.names(decoder_df), replacement = "")
    result=data.frame(row.names(decoder_df),round(as.numeric(decoder_df),3))
    colnames(result)=c("keyword","r")
    result=result[order(-result$r),]
    message("\u2713 \n")
    return(result)
  } else {
    stop("The neurosynth database was not found in the package data. Please try running `VWRfirstrun(requirement='neurosynth')`.")
  }
}  


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

#R-native replacement for the NiMARE 0.5.4 decoder using python for BrainStat's surface wrapper
# Bridge for the existing BrainStat _surf2vol result, not part of R decoding.
# Explicit py_to_r conversion is important for NumPy arrays and Python tuples.
nibabel_roi_to_R <- function(img) {
  cv <- function(x) if (inherits(x, "python.builtin.object"))
    reticulate::py_to_r(x) else x
  list(data = cv(img$get_fdata()), affine = cv(img$affine),
       zooms = as.numeric(unlist(cv(img$header$get_zooms())))[1:3])
}

.vwr_validate_roi <- function(masker) {
  if (!is.list(masker) || is.null(masker$data) || is.null(masker$affine))
    stop("masker must be list(data=3D array, affine=4x4 matrix, zooms=c(dx,dy,dz)).")
  d <- masker$data
  A <- masker$affine
  if (!(is.numeric(d) || is.logical(d)) || length(dim(d)) != 3L ||
      any(dim(d) < 1L) || any(!is.finite(d))) stop("ROI must be a finite 3D array.")
  if (!is.matrix(A) || !is.numeric(A) || !identical(dim(A), c(4L, 4L)) ||
      any(!is.finite(A)) || any(abs(A[4, ] - c(0, 0, 0, 1)) > 1e-8))
    stop("Invalid 4x4 voxel-to-world affine.")
  inv <- tryCatch(solve(A), error = function(e) stop("Singular ROI affine."))
  z <- masker$zooms
  if (is.null(z)) z <- sqrt(colSums(A[1:3, 1:3, drop = FALSE]^2))
  if (!is.numeric(z) || length(z) != 3L || any(!is.finite(z)) || any(z <= 0))
    stop("zooms must contain three positive voxel sizes in mm.")
  keep <- which(d != 0)
  if (!length(keep)) stop("The ROI is empty.")
  list(data = d != 0, affine = A, inverse = inv, zooms = z, keep = keep)
}

#' Construct an R ROI association decoder
#' @param masker List with data (3D array), affine (4x4 voxel-to-world matrix),
#'   and zooms (three header voxel sizes in mm; default affine column norms).
#'   Array element (1,1,1) has zero-based voxel coordinate (0,0,0).
#' @param kernel_transformer Only "MKDAKernel" is implemented.
#' @param feature_group Optional annotation prefix, with or without trailing __.
#' @param features Optional annotation names; unprefixed when group is supplied.
#' @param kernel__r Sphere radius, in mm; default 10.
#' @param kernel__value Value inside a sphere; default 1.
#' @return Environment exposing $fit(dataset), $transform(), $roi_values_,
#'   $ids_, and $features_. $transform() returns a data.frame with column r.
#' @noRd

ROIAssociationDecoder <- function(masker, kernel_transformer = "MKDAKernel",
                                  feature_group = NULL, features = NULL,
                                  kernel__r = 10, kernel__value = 1) {
  if (!identical(kernel_transformer, "MKDAKernel"))
    stop("This implementation supports MKDAKernel only.")
  if (length(kernel__r) != 1L || !is.finite(kernel__r) || kernel__r < 0 ||
      length(kernel__value) != 1L || !is.finite(kernel__value))
    stop("Invalid kernel radius/value.")
  roi <- .vwr_validate_roi(masker)
  shape <- dim(roi$data)
  # Same discrete sphere as NiMARE: header voxel sizes, inclusive radius.
  extent <- floor(kernel__r / roi$zooms)
  offsets <- as.matrix(expand.grid(lapply(extent, function(e) seq(-e, e))))
  offsets <- offsets[rowSums(sweep(offsets, 2L, roi$zooms, "*")^2) <=
                       kernel__r^2, , drop = FALSE]
  # Map full-volume linear positions to compact ROI positions (R column order).
  lookup <- integer(length(roi$data))
  lookup[roi$keep] <- seq_along(roi$keep)
  roi_ijk <- arrayInd(roi$keep, shape) - 1L
  lower <- apply(roi_ijk, 2L, min) - extent
  upper <- apply(roi_ijk, 2L, max) + extent
  self <- new.env(parent = emptyenv())
  fitted_annotations <- NULL
  self$fit <- function(dataset, drop_invalid = TRUE) {
    fitted_annotations <<- NULL
    self$roi_values_ <- self$ids_ <- self$features_ <- NULL
    a <- dataset$annotations
    c <- dataset$coordinates
    ids <- sort(unique(c$id), method = "radix")
    invalid <- union(setdiff(a$id, ids), setdiff(ids, a$id))
    if (length(invalid) && !drop_invalid)
      stop("Some IDs are missing coordinates or annotations.")
    ids <- ids[ids %in% a$id]
    if (length(ids) < 2L) stop("Need at least two studies with coordinates and annotations.")
    a <- a[match(ids, a$id), , drop = FALSE]
    c <- c[c$id %in% ids, , drop = FALSE]
    f <- setdiff(names(a), c("id", "study_id", "contrast_id"))
    if (!is.null(feature_group)) {
      if (length(feature_group) != 1L || is.na(feature_group)) stop("Invalid feature_group.")
      prefix <- paste0(sub("_+$", "", feature_group), "__")
      f <- f[startsWith(f, prefix)]
      if (!is.null(features)) {
        wanted <- paste0(prefix, features)
        if (any(!wanted %in% f)) stop("Unknown features in selected group.")
        f <- wanted
      }
    } else if (!is.null(features)) {
      if (any(!features %in% f)) stop("Unknown annotation features.")
      f <- features
    }
    if (!length(f) || anyDuplicated(f)) stop("No features or duplicate features selected.")
    if (!all(vapply(a[f], is.numeric, logical(1)))) stop("Features must be numeric.")
    Y <- as.matrix(a[f])
    if (any(!is.finite(Y))) stop("Non-finite annotation weights are not supported.")
    f <- f[colSums(Y > 0) > 0]
    if (!length(f)) stop("No features have positive weights in any retained study.")
    Y <- Y[, f, drop = FALSE]
    # NiMARE mm2vox casts to int: truncate toward zero, NOT round().
    ijk <- trunc((cbind(as.matrix(c[c("x", "y", "z")]), 1) %*%
                    t(roi$inverse))[, 1:3, drop = FALSE])
    groups <- split(seq_len(nrow(c)), factor(c$id, levels = ids))
    values <- numeric(length(ids))
    for (s in seq_along(ids)) {
      pts <- unique(ijk[groups[[s]], , drop = FALSE])
      eligible <- rowSums(sweep(pts, 2L, lower, ">=")) == 3L &
        rowSums(sweep(pts, 2L, upper, "<=")) == 3L
      pts <- pts[eligible, , drop = FALSE]
      covered <- logical(length(roi$keep))
      for (j in seq_len(nrow(pts))) {
        q <- sweep(offsets, 2L, pts[j, ], "+")
        inside <- rowSums(q >= 0) == 3L &
          rowSums(sweep(q, 2L, shape, "<")) == 3L
        q <- q[inside, , drop = FALSE]
        idx <- 1 + q[, 1] + shape[1] * (q[, 2] + shape[2] * q[, 3])
        hits <- lookup[idx]
        covered[hits[hits > 0L]] <- TRUE
      }
      values[s] <- kernel__value * sum(covered) / length(covered)
    }
    self$ids_ <- ids
    self$features_ <- f
    self$roi_values_ <- stats::setNames(values, ids)
    fitted_annotations <<- Y
    invisible(self)
  }
  self$transform <- function() {
    if (is.null(fitted_annotations)) stop("Call $fit(dataset) before $transform().")
    # Correlation across STUDIES, one result per feature. Constant vectors
    # return NaN (undefined), as in NiMARE's Pearson implementation.
    x <- unname(self$roi_values_)
    x <- x - mean(x)
    Y <- sweep(fitted_annotations, 2L, colMeans(fitted_annotations), "-")
    den <- sqrt(sum(x^2) * colSums(Y^2))
    r <- as.numeric(crossprod(x, Y)) / den
    r[den == 0] <- NaN
    finite <- is.finite(r)
    r[finite] <- pmax(-1, pmin(1, r[finite]))
    data.frame(r = r, row.names = self$features_, check.names = FALSE)
  }
  class(self) <- "VWR_ROIAssociationDecoder"
  self
}