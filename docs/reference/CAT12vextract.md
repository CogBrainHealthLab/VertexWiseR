# CAT12vextract

Extracts vertex-wise surface-based measures for each subject from a
[CAT12](https://neuro-jena.github.io/cat12-help/) preprocessed
directory, resampled to a 32k mesh, and stores them all as a single .RDS
file.

## Usage

``` r
CAT12vextract(
  sdirpath = "./",
  filename,
  measure = "thickness",
  subj_ID = TRUE,
  silent = FALSE,
  VWR_check = TRUE
)
```

## Arguments

- sdirpath:

  A string object containing the path to the CAT12 subjects preprocessed
  directory. Default is the current working directory ("./").

- filename:

  A string object containing the desired name of the output RDS file.
  Default is 'CAT12_measure.rds' in the R temporary directory
  (tempdir()).

- measure:

  A string object containing the name of the measure of interest.
  Options are 'thickness', 'depth', 'fractaldimension', 'gyrification',
  and 'toroGI20mm'. Default is 'thickness.'

- subj_ID:

  A logical object to determine whether to return a list object
  containing both subject ID and data matrix. Subject IDs are assumed to
  be the top directory names in the sdirpath.

- silent:

  A logical object to determine whether messages will be silenced. Set
  to 'FALSE' by default

- VWR_check:

  A boolean object specifying whether to check and validate system
  requirements. Default is TRUE.

## Value

A .RDSfile with a list containing 1. the list of subject IDs (first
element) and 2. a surface data matrix object (second element), or only a
data matrix object. The matrix has N subjects x M vertices dimensions
and can be readily used by VertexWiseR statistical analysis functions.
Each row corresponds to a subject (in the order they are listed in the
folder) and contains the left to right hemispheres' vertex-wise values.

## Details

The function searches inside the CAT12 preprocessed for 32k meshes
(.gii) with the user-selected measure, extracts the data from these
files, and organizes the left and right hemisphere vertex data for each
subject as rows in a N x 64984 data matrix within a .rds object. Python
and reticulate are required as the [NiBabel](https://nipy.org/nibabel/)
package is used to import .gii files outputted by CAT12.

## Examples

``` r
CAT12vextract(sdirpath="./", 
filename='thickness.rds', 
measure='thickness', 
subj_ID = TRUE, 
VWR_check=FALSE)
#> Non-interactive sessions need requirement checks
```
