# HIPvextract

Extracts hippocampal vertex-wise surface-based measures for each subject
in the 'HippUnfold' subjects directory, and stores it as a single .RDS
file.

## Usage

``` r
HIPvextract(sdirpath = "./", filename, measure = "thickness", subj_ID = TRUE)
```

## Arguments

- sdirpath:

  A string object containing the path to the 'HippUnfold' subjects
  directory. Default is the current working directory ("./").

- filename:

  A string object containing the desired name of the output RDS file.
  Default is 'hip_measure.rds' in the R temporary directory (tempdir()).

- measure:

  A string object containing the name of the measure of interest.
  Options are 'thickness','curvature','gyrification' and 'surfarea' (For
  more information see [the 'HippUnfold'
  documentation](https://hippunfold.readthedocs.io/en/latest/outputs/output_files.html#surface-metrics)).
  Default is thickness.

- subj_ID:

  A logical object stating whether to return a list object containing
  both subject ID and data matrix.

## Value

A .RDSfile with a list containing 1. the list of subject IDs (first
element) and 2. a surface data matrix object (second element), or only a
data matrix object. The matrix has N subjects x M vertices dimensions
and can be readily used by VertexWiseR statistical analysis functions.
Each row corresponds to a subject (in the order they are listed in the
folder) and contains the left to right hemispheres' hippocampal
vertex-wise values.

## Details

The function searches for the hippocampal surface data by listing out
files with certain suffixes, extract the data from these files, and
organize the left and right hippocampal vertex data for each subject as
rows in a N x 14524 data matrix within a .rds object.

## Examples

``` r
HIPvextract(sdirpath = "./", filename = paste0(tempdir(),"/hip_data.RDS"), measure = "thickness") 
#> HippUnfold data could not be found in the set sdirpath
```
