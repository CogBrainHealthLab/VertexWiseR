# Decode surface data

Correlates the significant clusters of an earlier vertex-wise analysis
with a database of task-based fMRI and voxel-based morphometric
statistical maps and associate them with relevant key words. Decoding
currently works with surfaces in fsaverage5 space only."

## Usage

``` r
decode_surf_data(surf_data, contrast = "positive", VWR_check = TRUE)
```

## Arguments

- surf_data:

  A numeric vector or object containing the surface data, in fsaverage5
  (1 x 20484 vertices) or fsLR32k (1 x 64984 vertices) space. It can
  only be one row of vertices (not a cohort surface data matrix).

- contrast:

  A string object indicating whether to decode the positive or negative
  mask ('positive' or 'negative')

- VWR_check:

  A boolean object specifying whether to check and validate system
  requirements. Default is TRUE.

## Value

A data.frame object listing the keywords and their Pearson's R values

## Details

An internal image decoding is used, reproducing the
['NiMARE'](https://nimare.readthedocs.io/en/stable/index.html) python
module's ROIAssociationDecoder in base R. The function also downloads
the ['Neurosynth'
database](https://github.com/neurosynth/neurosynth-data) in the
package's inst/extdata directory (converted to .RDS, ~7.2 MB) for the
analysis.

## Examples

``` r
CTv = rbinom(20484, 1, 0.001) 
decoding = decode_surf_data(CTv, 'positive', VWR_check=FALSE);
#> Converting and interpolating the surface data ... 
#> ✓ 
#> Correlating input image with images in the neurosynth database. This may take a while ... 
#> ✓ 
head(decoding)
#>         keyword     r
#> 195    encoding 0.052
#> 538   retrieval 0.051
#> 311    language 0.046
#> 510 recognition 0.044
#> 202    episodic 0.043
#> 381       motor 0.041
```
