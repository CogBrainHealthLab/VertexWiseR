# Surface to atlas

Returns the mean or sum of vertex-wise surface data for each ROI of a
selected atlas

## Usage

``` r
surf_to_atlas(surf_data, atlas, mode = "mean")
```

## Arguments

- surf_data:

  A N x V matrix object containing the surface data (N row for each
  subject, V for each vertex), in fsaverage5 (20484 vertices),
  fsaverage6 (81924 vertices), fslr32k (64984 vertices) or hippocampal
  (14524 vertices) space. See also Hipvextract(), SURFvextract() or
  FSLRvextract output formats.

- atlas:

  A numeric integer object corresponding to the atlas of interest.
  1=Desikan, 2=Destrieux-148, 3=Glasser-360, 4=Schaefer-100,
  5=Schaefer-200, 6=Schaefer-400. Set to `1` by default. This argument
  is ignored for hippocampal surfaces.

- mode:

  A string indicating whether to extract the sum ('sum') or the average
  ('mean') of the ROI vertices values. Default is 'mean'.

## Value

A matrix object with ROI as column and corresponding average vertex-wise
values as row

## Details

The function currently works with the aparc/Desikan-Killiany-70,
Destrieux-148, Glasser-360, Schaefer-100, Schaefer-200, Schaefer-400
atlases. ROI to vertex mapping data were obtained from the ['ENIGMA
toolbox'](https://github.com/MICA-MNI/ENIGMA/tree/master/enigmatoolbox/datasets/parcellations)
; data for Destrieux came from ['Nilearn' 's
nilearn.datasets.fetch_atlas_surf_destrieux](https://github.com/nilearn/nilearn/blob/a366d22e426b07166e6f8ce1b7ac6eb732c88155/nilearn/datasets/atlas.py)

For hippocampal data, the function currently works with the "bigbrain"
10-parcels atlas integrated in 'HippUnfold.' See also
[doi:10.1016/j.neuroimage.2019.116328](https://doi.org/10.1016/j.neuroimage.2019.116328)
.

## See also

[`atlas_to_surf`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/atlas_to_surf.md)

## Examples

``` r
CTv = runif(20484,min=0, max=100)
parcel_data = surf_to_atlas(CTv, 1)
```
