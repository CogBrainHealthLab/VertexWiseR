# Thresholding TFCE output

Threshold TFCE maps from the TFCE_vertex_analysis() output and
identifies significant clusters at the desired threshold.

## Usage

``` r
TFCE_threshold(TFCEoutput, p = 0.05, atlas = 1, k = 20, VWR_check = TRUE)
```

## Arguments

- TFCEoutput:

  An object containing the output from TFCE_vertex_analysis()

- p:

  A numeric object specifying the p-value to threshold the results
  (Default is 0.05)

- atlas:

  An integer corresponding to the atlas of interest. 1=Desikan,
  2=Destrieux-148, 3=Glasser-360, 4=Schaefer-100, 5=Schaefer-200,
  6=Schaefer-400. Set to `1` by default. This argument is ignored for
  hippocampal surfaces. For applicable SubCortexMesh surfaces, 1=default
  base ROI, 2=anatomical atlas.

- k:

  Cluster-forming threshold (Default is 20)

- VWR_check:

  A boolean object specifying whether to check and validate system
  requirements. Default is TRUE.

## Value

A list object containing the cluster level results, unthresholded t-stat
map, thresholded t-stat map, and positive, negative and bidirectional
cluster maps.

## Examples

``` r
model1_TFCE=readRDS(file = url(paste0("https://github.com/CogBrainHealthLab",
"/VertexWiseR/blob/main/inst/demo_data/model1_TFCE.rds?raw=TRUE")))

TFCEanalysis_output=TFCE_threshold(model1_TFCE, p=0.05, atlas=1,
VWR_check=FALSE)
TFCEanalysis_output$cluster_level_results
#> $`Positive contrast`
#> [1] "No significant clusters"
#> 
#> $`Negative contrasts`
#>   clusid nverts      P   X     Y     Z tstat              region
#> 1      1   8098 <0.001  47   4.0 -16.6 12.64 rh-superiortemporal
#> 2      2   7728 <0.001 -34 -25.7  16.2 14.23           lh-insula
#> 
```
