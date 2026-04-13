# Surface plotter

Plots surface data in a grid with one or multiple rows in a .png file

## Usage

``` r
plot_surf(
  surf_data,
  filename,
  title = "",
  surface = "inflated",
  cmap,
  limits,
  colorbar = TRUE,
  size,
  zoom,
  transparent_bg = FALSE,
  show_nan = T,
  alpha = 1,
  smooth_mesh = 0,
  show.plot.window = TRUE,
  VWR_check = TRUE
)
```

## Arguments

- surf_data:

  A numeric vector (length of V) or a matrix (N rows x V columns), where
  N is the number of subplots, and V is the number of vertices. It can
  be the output from SURFvextract(), CAT12vextract(), FSLRvextract(),
  HIPvextract(), SCMvextract(), as well as masks or vertex-wise results
  outputted by analyses functions. Alternatively, atlas ROI values as
  supported by atlas_to_surf() may be given.

- filename:

  A string object containing the desired name of the output .png.
  Default is 'plot.png' in the R temporary directory (tempdir()).Only
  filenames with a .png extension are allowed.

- title:

  A string object for setting the title in the plot. Default is none.
  For titles that too long to be fully displayed within the plot, we
  recommend splitting them into multiple lines by inserting "\n".

- surface:

  A string object containing the name of the type of cortical surface
  background rendered (only applicable on fsaverage5 and fsaverage6).
  Possible options include "white", "smoothwm","pial" and "inflated"
  (default).

- cmap:

  A string object specifying the name of an existing colormap or a
  vector of hexadecimal color codes to be used as a custom colormap. The
  names of existing colormaps are listed in the ['Matplotlib' plotting
  library](https://matplotlib.org/stable/gallery/color/colormap_reference.html).

  Default cmap is set to `"Reds"` for positive values, `"Blues_r"` for
  negative values and `"RdBu"` when both positive and negative values
  exist.

- limits:

  A combined pair of numeric vector composed of the lower and upper
  color scale limits of the plot. If the limits are specified, the same
  limits will be applied to all subplots. When left unspecified, the
  same symmetrical limits c(-max(abs(surf_dat),max(abs(surf_dat))) will
  be used for all subplots. If set to NULL, each subplot will have its
  own limits corresponding to their min and max values

- colorbar:

  A logical object stating whether to include a color bar in the plot or
  not (default is TRUE).

- size:

  A combined pair of numeric vector indicating the image dimensions
  (width and height in pixels). Default is c(1920,400) for whole-brain
  surface and c(400,200) for HippUnfold hippocampal surface.

- zoom:

  A numeric value for adjusting the level of zoom on the figures.
  Default is 1.25 for whole-brain surface and 1.20 for HippUnfold
  hippocampal surface.

- transparent_bg:

  A logical object to determine if the background of the image is set to
  transparent (Default is FALSE).

- show_nan:

  A logical object to determine if the NaN vertices are to be plotted
  (Default is TRUE).

- alpha:

  A numeric object between 0 and 1 to determine the opacity of the
  non-empty vertices. Note that this is not a true opacity setting, it
  will blend the colour into that of the NaN vertices (white if show_nan
  is FALSE)

- smooth_mesh:

  A numeric object stating the number of iterations of smoothing to make
  the surface appear smoother. Not applicable for HippUnfold hippocampal
  data. Default is 0.

- show.plot.window:

  A logical object to determine if the generated plot is to be shown
  within RStudio's plot window

- VWR_check:

  A boolean object specifying whether to check and validate system
  requirements. Default is TRUE.

## Value

Outputs the plot as a .png image

## Examples

``` r
results = runif(20484,min=0, max=1);
plot_surf(surf_data = results, 
filename=paste0(tempdir(),"/output.png"),
title = 'Cortical thickness', 
surface = 'inflated', cmap = 'Blues',
VWR_check=FALSE)
#> Non-interactive sessions need requirement checks
```
