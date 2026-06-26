# fslr to fsaverage5

Remaps vertex-wise surface data in fslr space to fsaverage5 space using
the nearest neighbor approach

## Usage

``` r
fslr_to_fs5(surf_data)
```

## Arguments

- surf_data:

  A N x V matrix object containing the surface data (N row for each
  subject, V for each vertex), in fslr (64984 vertices) space. See also
  SURFvextract() output format.

## Value

A matrix object containing vertex-wise surface data mapped in fsaverage5
space

## See also

[`fs5_to_fslr`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/fs5_to_fslr.md)

## Examples

``` r
surf_data = runif(64984,min=0, max=100);
fs5_data=fslr_to_fs5(surf_data)
#> Error in fslr_to_fs5(surf_data): could not find function "fslr_to_fs5"
```
