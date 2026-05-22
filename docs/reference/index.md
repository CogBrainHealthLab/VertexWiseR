# Package index

## Demo & tutorials

- [`getting_started()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/getting_started.md)
  : Getting started with tutorials

## Requirements check

- [`VWRfirstrun()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/VWRfirstrun.md)
  : VertexWiseR system requirements installation

## Surface data extraction

- [`HIPvextract()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/HIPvextract.md)
  : HIPvextract
- [`SURFvextract()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/SURFvextract.md)
  : SURFvextract
- [`FSLRvextract()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/FSLRvextract.md)
  : FSLRvextract
- [`CAT12vextract()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/CAT12vextract.md)
  : CAT12vextract
- [`DTSERIESvextract()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/DTSERIESvextract.md)
  : DTSERIESvextract
- [`SCMvextract()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/SCMvextract.md)
  : SCMvextract

## Smoothing

- [`smooth_surf()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/smooth_surf.md)
  : Smooth surface

## Main analysis functions

### Vertex-Wise analysis

- [`RFT_vertex_analysis()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/RFT_vertex_analysis.md)
  : Vertex-wise analysis with random field theory cluster correction
- [`TFCE_vertex_analysis()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/TFCE_vertex_analysis.md)
  : Vertex-wise analysis with threshold-free cluster enhancement (fixed
  effect)
- [`TFCE_vertex_analysis_mixed()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/TFCE_vertex_analysis_mixed.md)
  : Vertex-wise analysis with threshold-free cluster enhancement (mixed
  effect)
- [`TFCE_threshold()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/TFCE_threshold.md)
  : Thresholding TFCE output

### Plotting

- [`plot_surf()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/plot_surf.md)
  : Surface plotter
- [`plot_surf3d()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/plot_surf3d.md)
  : 3D Surface plotter
- [`plot_overlay_surf()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/plot_overlay_surf.md)
  : Surface overlay plotter

## Other functions

### Decoding surface data

- [`decode_surf_data()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/decode_surf_data.md)
  : Decode surface data

### Conversion functions

- [`fs5_to_fs6()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/fs5_to_fs6.md)
  : fsaverage5 to fsaverage6
- [`fs6_to_fs5()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/fs6_to_fs5.md)
  : fsaverage6 to fsaverage5
- [`surf_to_vol()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/surf_to_vol.md)
  : Surface to volume
- [`atlas_to_surf()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/atlas_to_surf.md)
  : Atlas to surface
- [`surf_to_atlas()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/surf_to_atlas.md)
  : Surface to atlas

### Descriptive statistics extraction

- [`fs_stats()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/fs_stats.md)
  : fs_stats()

## Data objects

### Object classes

- [`MNIsurface-class`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/MNIsurface-class.md)
  : MNI surface map object
- [`edgelist-class`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/edgelist-class.md)
  : Edge list object
- [`ROImap-class`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/ROImap-class.md)
  : Region-of-Interest mapping object

### Edges

- [`edgelist_hip`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/edgelist_hip.md)
  : List of edges for the hippocampal template

### Template maps

- [`fs6_to_fs5_map`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/fs6_to_fs5_map.md)
  : fsaverage6 template object for nearest neighbor conversion in
  fs6_to_fs5()
- [`hip_points_cells`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/hip_points_cells.md)
  : points and cells data required to build the hippocampus surface
  template

### Parcellation maps

- [`ROImap_fs5`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/ROImap_fs5.md)
  : Atlas parcellations of fsaverage5
- [`ROImap_fs6`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/ROImap_fs6.md)
  : Atlas parcellations of fsaverage6
- [`ROImap_fslr32k`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/ROImap_fslr32k.md)
  : Atlas parcellations of FS_LR32k
- [`ROImap_hip`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/ROImap_hip.md)
  : Atlas parcellations of the hippocampus (CITI168)
