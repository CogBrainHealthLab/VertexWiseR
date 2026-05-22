# Surface to volume

Converts surface data to volumetric data (.nii file)

## Usage

``` r
surf_to_vol(surf_data, filename, VWR_check = TRUE)
```

## Arguments

- surf_data:

  A numeric vector or object containing the surface data, either in
  fsaverage5 (1 x 20484 vertices) or fsLR32k (1 x 64984 vertices) space.
  It can only be one row of vertices (not a cohort surface data matrix).

- filename:

  A string object containing the desired name of the output .nii file
  (default is 'output.nii' in the R temporary directory (tempdir())).

- VWR_check:

  A boolean object specifying whether to check and validate system
  requirements. Default is TRUE.

## Value

A .nii volume file

## Examples

``` r
CTv = runif(20484,min=0, max=100);
surf_to_vol(CTv, filename = paste0(tempdir(),'/volume.nii'), VWR_check=FALSE)
#> Non-interactive sessions need requirement checks
```
