# fsaverage5 to fslr

Remaps vertex-wise surface data in fsaverage5 space to fslr space using
the nearest neighbor approach

## Usage

``` r
fs5_to_fslr(surf_data)
```

## Arguments

- surf_data:

  A N x V matrix object containing the surface data (N row for each
  subject, V for each vertex), in fsaverage5 (20484 vertices) space. See
  also SURFvextract() output format.

## Value

A matrix object containing vertex-wise surface data mapped in fslr space

## See also

[`fslr_to_fs5`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/fslr_to_fs5.md)

## Examples

``` r
CTv = runif(20484,min=0, max=100);
CTv_fslr = fs5_to_fslr(CTv);
#> Error in fs5_to_fslr(CTv): could not find function "fs5_to_fslr"
```
