# FSLRvextract

Extracts vertex-wise surface-based measures for each subject from
sMRIprep, fMRIprep, ASLprep, HCP processing or XCP-D preprocessed
directories, and stores it as a single .RDS file.

## Usage

``` r
FSLRvextract(
  sdirpath = "./",
  filename,
  dscalar,
  subj_ID = TRUE,
  silent = FALSE,
  VWR_check = TRUE
)
```

## Arguments

- sdirpath:

  A string object containing the path to the sMRIprep, fMRIprep,
  ASLprep, HCP preprocessed directory. Default is the current working
  directory ("./").

- filename:

  A string object containing the desired name of the output RDS file.
  Default is 'fslr32k.rds' in the R temporary directory (tempdir()).

- dscalar:

  A string object containing the filename suffix of the dscalar file.
  These dscalar files are named differently depending on the
  preprocessing pipeline used. Examples of filename suffixes are shown
  below

  - `.thickness_MSMAll.32k_fs_LR.dscalar.nii` (HCP MSMAll pipeline)

  - `.sulc_MSMAll.32k_fs_LR.dscalar.nii` (HCP MSMAll pipeline)

  - `.thickness.32k_fs_LR.dscalar.nii` (HCP legacy pipeline)

  - `.sulc.32k_fs_LR.dscalar.nii` (HCP legacy pipeline)

  - `_space-fsLR_den-91k_thickness.dscalar.nii`
    (sMRIprep/fMRIprep/ASLprep; using the `--cifti-output 91k` flag)

  - `_space-fsLR_den-91k_curv.dscalar.nii` (sMRIprep/fMRIprep/ASLprep;
    using the `--cifti-output 91k` flag)

  - `_space-fsLR_den-91k_sulc.dscalar.nii` (sMRIprep/fMRIprep/ASLprep;
    using the `--cifti-output 91k` flag)

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

A .RDSfile with a list containing 1. the list of subject IDs (first
element) and 2. a surface data matrix object (second element), or only a
data matrix object. The matrix has N subjects x M vertices dimensions
and can be readily used by VertexWiseR statistical analysis functions.
Each row corresponds to a subject (in the order they are listed in the
folder) and contains the left to right hemispheres' vertex-wise values.

## Details

The function searches in the preprocessed directory by listing out files
with certain suffixes, extracts the data from these files, and organizes
the left and right hemisphere vertex data for each subject as rows in a
N x 64984 data matrix within a .rds object.

## Examples

``` r
dat_fslr32k=FSLRvextract(sdirpath="./", 
filename="dat_fslr32k.rds",
dscalar=".thickness_MSMAll.32k_fs_LR.dscalar.nii", 
subj_ID = TRUE, 
silent=FALSE,
VWR_check=FALSE)
#> Non-interactive sessions need requirement checks
```
