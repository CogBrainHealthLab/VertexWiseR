# ASEGvextract

Extracts surface-based measures from ASeg subcortical segmentations in
each subject. These must be specifically computed from a FreeSurfer
output subjects directory using the python toolbox "SubCortexMesh".
ASEGvextract() reads through the metrics .vtk files outputted by
SubCortexMesh (surface_metrics/ directory), extracts the vertex-wise
values, and stores them as a .rds file, merging left and right
hemispheres when applicable, for each (selected) subcortical region.

## Usage

``` r
ASEGvextract(
  sdirpath = "./",
  outputdir,
  measure = "thickness",
  roilabel,
  subj_ID = TRUE,
  silent = FALSE,
  VWR_check = TRUE
)
```

## Arguments

- sdirpath:

  A string object containing the path to the 'SubCortexMesh' surface
  metrics directory, containing surface-based metrics (.vtk files).
  Default is the current working directory ("./").

- outputdir:

  A string object containing the path of the directory where all ROI RDS
  files will be stored. Default is 'subcortices' in the R temporary
  directory (tempdir()).

- measure:

  A string object containing the name of the measure of interest.
  Options include thickness, and surfarea. Default is thickness.

- roilabel:

  A string object or vector of string objects containing the name(s) of
  the subcortical region(s) to extract (both hemispheres automatically
  included when applicable). E.g. "thalamus" or c("thalamus",
  "caudate"). Default is all regions in sdirpath.

- subj_ID:

  A logical object to determine whether to return a list object
  containing both subject ID and data matrix.

- silent:

  A logical object to determine whether messages will be silenced. Set
  to 'FALSE' by default

- VWR_check:

  A boolean object specifying whether to check and validate system
  requirements. Default is TRUE.

## Value

A directory containing - for each bilateral subcortical region
separately - a .RDS files, each with a list containing 1. the list of
subject IDs (first element) and 2. a surface data matrix object (second
element), or only a data matrix object. Each matrix has N subjects x M
vertices dimensions and can be used readily by VertexWiseR statistical
analysis functions. Each row corresponds to a subject (in the order they
are listed in the folder) and contains the left to right hemispheres'
vertex-wise values (if applicable).

## Details

The function looks for standard subject folder names (starting "sub-")
and surface file names as defined by SubCortexMesh (starting
"left-","right-", or "brain-stem"). ASEGvextract() requires reticulate
and python to read the .vtk files.

Number of vertices in the (bilateral) matrix for each
region-of-interest:

- 2044: accumbens area

- 3430: amygdala

- 6940: caudate nuclei

- 39214: cerebellum

- 8132: hippocampus

- 3200: pallidum

- 8394: putamen

- 7768: thalamus

- 7144: ventral diencephalon

- 9452: brain stem

- 95718: all aseg ROIs

## Examples

``` r
ASEGvextract(sdirpath = "subcortexmesh_output_metrics", 
outputdir=paste0(tempdir(), "\\subcortices"), measure="surfarea") 
#> Checking for VertexWiseR system requirements ... 
#> Subcortical meshes could not be found in the set sdirpath.
```
