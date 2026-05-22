# Atlas to surface

Maps average parcellation surface values (e.g. produced with the
surf_to_atlas() function) to the fsaverage5, fsaverage6 or fslr32k space

## Usage

``` r
atlas_to_surf(parcel_data, template)
```

## Arguments

- parcel_data:

  A matrix or vector object containing average surface measures for each
  region of interest, see the surf_to_atlas() output format.

- template:

  A string object stating the surface space on which to map the data
  ('fsaverage5', 'fsaverage6', 'fslr32k', 'CIT168' (hippocampal)).

## Value

A matrix or vector object containing vertex-wise surface data mapped in
fsaverage5, fsaverage6, fslr32k, or CIT168 space

## Details

The function currently supports the Desikan-Killiany-70, Schaefer-100,
Schaefer-200, Schaefer-400, Glasser-360, or Destrieux-148 atlases for
cortical surfaces, and the 'bigbrain' 10-parcels atlas for hippocampal
surfaces. ROI to vertex mapping data for 1 to 4 were obtained from the
['ENIGMA
toolbox'](https://github.com/MICA-MNI/ENIGMA/tree/master/enigmatoolbox/datasets/parcellations)
; and data for 5 from ['Nilearn' 's
nilearn.datasets.fetch_atlas_surf_destrieux](https://github.com/nilearn/nilearn/blob/a366d22e426b07166e6f8ce1b7ac6eb732c88155/nilearn/datasets/atlas.py)
. atlas_to_surf() will automatically detect the atlas based on the
number of columns.

## See also

[`surf_to_atlas`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/surf_to_atlas.md)

## Examples

``` r
parcel_data = t(runif(100,min=0, max=100));
surf_data = atlas_to_surf(parcel_data, template='fsaverage5');
```
