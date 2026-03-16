# fsaverage5 to fsaverage6

Remaps vertex-wise surface data in fsaverage5 space to fsaverage6 space
using the nearest neighbor approach

## Usage

``` r
fs5_to_fs6(surf_data)
```

## Arguments

- surf_data:

  A N x V matrix object containing the surface data (N row for each
  subject, V for each vertex), in fsaverage5 (20484 vertices) space. See
  also SURFvextract() output format.

## Value

A matrix object containing vertex-wise surface data mapped in fsaverage6
space

## See also

[`fs6_to_fs5`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/fs6_to_fs5.md)

## Examples

``` r
CTv = runif(20484,min=0, max=100);
CTv_fs6 = fs5_to_fs6(CTv);
```
