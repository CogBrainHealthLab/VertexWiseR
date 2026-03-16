# 3D Surface plotter

Plots 3D cortical surfaces

## Usage

``` r
plot_surf3d(
  surf_data,
  surf_color = "grey",
  cmap,
  limits,
  atlas = 1,
  hemi = "b",
  medial_gap = 0,
  orientation_labels = TRUE,
  plot_grid = TRUE,
  transparent_bg = FALSE,
  smooth_mesh = 0,
  VWR_check = TRUE
)
```

## Arguments

- surf_data:

  A numeric vector (length of V), where V is the number of vertices. It
  can be the output from SURFvextract(), CAT12vextract(),
  FSLRvextract(), HIPvextract(), ASEGvextract().

- surf_color:

  color of the cortical surface for NA values. Set to `'grey'` by
  default. A RGBA string can also be given, with A as the opacity (0 to
  1), e.g. for transparent grey: surf_color="rgba(100,100,100,0.5)".

- cmap:

  A string vector containing 2 to 4 color names/codes specifying the
  colors to be used for the color scale; or a single string object with
  the name of a color map listed in
  [`RColorBrewer::display.brewer.all()`](https://rdrr.io/pkg/RColorBrewer/man/ColorBrewer.html).
  If none are specified, appropriate colors will be automatically
  selected according to `range(surf_data)`

- limits:

  A combined pair of numeric vector composed of the lower and upper
  color scale limits of the plot. When left unspecified, the symmetrical
  limits `c(-max(abs(surf_dat),max(abs(surf_dat)))` will be used.

- atlas:

  atlas used for identifying region labels. 1=Desikan, 2=Destrieux-148,
  3=Glasser-360, 4=Schaefer-100, 5=Schaefer-200, 6=Schaefer-400. Set to
  `1` by default. This argument is ignored for HippUnfold hippocampal
  surfaces and SubCortexMesh subcortical surface.

- hemi:

  A string specifying the hemisphere to plot. Possible values are `l`
  (left), `r` (right) or `b` (both).

- medial_gap:

  A numeric value specifying the amount of gap (in MNI coordinate units)
  to separate the left and right hemispheres. Set to `0` (no gap between
  hemispheres) by default. In order to view the medial surfaces clearly,
  it is recommended that this value is set to `20`. This argument is
  ignored if `hemi!='b'`

- orientation_labels:

  A boolean object specifying if orientation labels are to be displayed.
  Set to `TRUE` by default

- plot_grid:

  A boolean object specifying whether to plot the orientation grid or
  not (default is `TRUE`).

- transparent_bg:

  A boolean object specifying whether to get a transparent background
  upon saving the image (default is `FALSE`, white background).

- smooth_mesh:

  A numeric object stating the number of iterations of (Laplacian)
  smoothing to make the surface appear smoother. Not applicable for
  HippUnfold hippocampal data. Default is 0.

- VWR_check:

  A boolean object specifying whether to check and validate system
  requirements. Default is TRUE.

## Value

a plot_ly object

## Examples

``` r
surf_data = runif(20484);
plot_surf3d(surf_data = surf_data, VWR_check=FALSE)
#> Non-interactive sessions need requirement checks
```
