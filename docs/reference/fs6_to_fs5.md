# fsaverage6 to fsaverage5

Remaps vertex-wise surface data in fsaverage6 space to fsaverage5 space
using the nearest neighbor approach

## Usage

``` r
fs6_to_fs5(surf_data)
```

## Arguments

- surf_data:

  A N x V matrix object containing the surface data (N row for each
  subject, V for each vertex), in fsaverage6 (81924 vertices) space. See
  also SURFvextract() output format.

## Value

A matrix object containing vertex-wise surface data mapped in fsaverage5
space

## See also

[`fs5_to_fs6`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/fs5_to_fs6.md)

## Examples

``` r
surf_data = runif(81924,min=0, max=100);
fs5_data=fs6_to_fs5(surf_data)
```
