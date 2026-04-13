# Example analyses with VertexWiseR - Example 3

------------------------------------------------------------------------

This article introduces vertex-wise statistical modelling of subcortical
surfaces available since update v1.5.0. The subcortical surfaces it is
currently compatible with are either based on the FSL FIRST (Patenaude
et al. 2011) or on the FreeSurfer “ASeg” (Fischl 2012) subcortical
volumetric segmentations. The volumes are converted to 3D meshes and
surface-based metrics (thickness, surface area, curvature) produced
using the python module
[SubCortexMesh](https://github.com/chabld/SubCortexMesh). VertexWiseR
can extract the surface measures obtained with the latter package.

## Requirements

Analysis and manipulation of the subcortical surfaces require additional
utility data to be downloaded and stored in the R package’s external
data (`system.file('extdata',package='VertexWiseR')`). A prompt will
automatically be triggered to assist downloading it whenever subcortical
surfaces are inputted (users can trigger it themselves by typing
`VertexWiseR:::scm_database_check(template='fsaverage')`). The size
depends on the template used (~19.9 MB for FreeSurfer ‘fsaverage’, ~17.3
MB for FSL ‘fslfirst’).

## Extracting ASeg subcortical surfaces

SubCortexMesh’s surface metric output in template space can be extracted
by VertexWiseR as in the other [surface extraction
functions](https://cogbrainhealthlab.github.io/VertexWiseR/articles/VertexWiseR_surface_extraction.html),
generating compact (bilateral) matrices for all applicable subcortical
surfaces. Here is an example of code which was run to obtain the demo
data:

``` r
aseg_CTv = SCMvextract(sdirpath = 'surface_metrics/', 
                        outputdir = "aseg_matrices/",
                        measure = c('thickness','curvature'), 
                        roilabel = c('thalamus','caudate'), 
                        template='fsaverage',
                        subj_ID = TRUE)
```

## Example analysis on the thalamus and caudate nuclei

The analysis will use surface data already extracted in R from a
preprocessed subjects directory, which we make available so you do not
need to preprocess a sample yourself. To obtain it, we had thickness and
curvature metrics data from a SubCortexMesh output directory based on
FreeSurfer segmentations of the SUDMEX_CONN dataset (Garza-Villarreal et
al. 2017).

The demo data (~216 MB) can be downloaded from the package’s github
repository with the following function:

``` r
#This will save the demo_data directory in a temporary directory (tempdir(), but you can change it to your own path)
demodata=VertexWiseR:::dl_demo(path=tempdir(), quiet=TRUE)
```

As in the SCMvextract() code above, metrics from the bilateral thalami
and caudate nuclei were extracted. The surface data can be loaded that
way:

``` r
thalamus_thickness = readRDS(file=paste0(demodata,"/SUDMEX_CONN_thalamus_thickness.rds"))
thalamus_curvature = readRDS(file=paste0(demodata,"/SUDMEX_CONN_thalamus_curvature.rds"))

caudate_thickness = readRDS(file=paste0(demodata,"/SUDMEX_CONN_caudate_thickness.rds"))
caudate_curvature = readRDS(file=paste0(demodata,"/SUDMEX_CONN_caudate_curvature.rds"))
```

To load the corresponding behavioural data (all subjects
indiscriminately):

``` r
beh_data = readRDS(paste0(demodata,"/SUDMEX_CONN_behdata.rds"))
```

Here, we are interested in the effect of group (cocaine use disorder
(CUD) versus control) on the thalami and caudate nuclei (replicating a
similar analysis that was done by (Xu, Xu, and Guo 2023)). We can run a
vertex-wise analysis model with random field theory-based cluster
correction, testing for the effect of group on each region and each
metric:

``` r
for (surf_data in c('thalamus_thickness', 'thalamus_curvature',
                    'caudate_thickness', 'caudate_curvature'))
{   
model=RFT_vertex_analysis(model = as.factor(beh_data$group),
                          contrast = as.factor(beh_data$group),
                          surf_data = get(surf_data))
assign(paste0('model_',surf_data),model)
}
```

``` r
model_thalamus_thickness$cluster_level_results
```

    ## $`Positive contrast`
    ## [1] "No significant clusters"
    ## 
    ## $`Negative contrast`
    ##   clusid nverts      P     X     Y     Z tstat         region
    ## 1      1    213 <0.001 106.1 132.1 163.4 -3.05 Right Thalamus
    ## 2      2    186 <0.001 153.8 128.7 163.5 -3.62  Left Thalamus
    ## 3      3    162  0.001 132.3 133.0 148.4 -3.55  Left Thalamus
    ## 4      4    127  0.002 115.3 118.5 161.0 -3.35 Right Thalamus
    ## 5      5    171  0.002 136.0 112.5 141.9 -3.20  Left Thalamus

``` r
model_thalamus_curvature$cluster_level_results
```

    ## $`Positive contrast`
    ##   clusid nverts      P     X     Y     Z tstat        region
    ## 1      1    119 <0.001 139.1 113.9 135.3  3.45 Left Thalamus
    ## 2      2     62  0.005 134.0 130.5 154.0  4.13 Left Thalamus
    ## 3      3     45  0.014 153.0 133.0 157.5  3.60 Left Thalamus
    ## 
    ## $`Negative contrast`
    ##   clusid nverts      P     X     Y     Z tstat         region
    ## 1      1    119 <0.001 143.5 116.7 133.0 -3.57  Left Thalamus
    ## 2      2     72  0.024 122.5 116.9 154.3 -4.19 Right Thalamus
    ## 3      3     63  0.025 135.0 126.6 159.2 -3.42  Left Thalamus
    ## 4      4     85  0.031 128.0 122.5 142.7 -3.96 Right Thalamus

``` r
model_caudate_thickness$cluster_level_results
```

    ## $`Positive contrast`
    ## [1] "No significant clusters"
    ## 
    ## $`Negative contrast`
    ##   clusid nverts      P     X     Y     Z tstat        region
    ## 1      1     77 <0.001 112.9 133.0 116.7 -3.44 Right Caudate
    ## 2      2     74  0.004 107.3 120.9 116.1 -2.93 Right Caudate
    ## 3      3     74  0.046 119.5 128.7 110.0 -3.61 Right Caudate

``` r
model_caudate_curvature$cluster_level_results
```

    ## $`Positive contrast`
    ## [1] "No significant clusters"
    ## 
    ## $`Negative contrast`
    ##   clusid nverts    P     X     Y     Z tstat        region
    ## 1      1     72 0.01 121.3 132.2 116.7 -3.63 Right Caudate

We can plot them together in one summary plot by binding the maps:

``` r
thalamusmaps=rbind(model_thalamus_thickness$thresholded_tstat_map,
                  model_thalamus_curvature$thresholded_tstat_map)

caudatemaps= rbind(model_caudate_thickness$thresholded_tstat_map,
              model_caudate_curvature$thresholded_tstat_map)
```

``` r
thalamusplot=plot_surf(thalamusmaps, 
          filename='SUDMEX_CONN_thalamus_tstatmaps.png',
          title=c('Thickness','Curvature'),
          smooth_mesh=20,
          show.plot.window=TRUE)
```

![Significant clusters after RFT correction on the
thalamus](VertexWiseR_Example_3_files/figure-html/unnamed-chunk-8-1.png)

``` r
caudateplot=plot_surf(caudatemaps, 
          filename='SUDMEX_CONN_caudate_tstatmaps.png',
          title=c('Thickness','Curvature'), 
          smooth_mesh=20,
          show.plot.window=TRUE)
```

![Significant clusters after RFT correction on the caudate
nuclei](VertexWiseR_Example_3_files/figure-html/unnamed-chunk-9-1.png)

The outcome shows that compared to controls, CUD patients had decreases
in thickness of the bilateral thalami, and various shapes alterations in
curvature.

Note regarding anatomical orientations: for the amygdala, brain-stem and
cerebellar plots, the upper part of the image is the superior part of
the brain and lower the inferior. While for the rest, more horizontal
structures, upper means anterior and lower means posterior, where the
second and third panels are medial views (looking from the center toward
the side/exterior of the brain) and vice versa for the first and fourth
panel. You may use plot_surf3d() to get a more explicit orientation.

## Analysing all subcortices together

SubCortexMesh has a function to merge all subcortices, assuming all of
them are available, into one single object for each metric. For example,
the “allaseg” can be outputted in SubCortexMesh’s surface_metrics/
directory, and extracted via SCMvextract():

``` r
allaseg_thickness_data = SCMvextract(sdirpath = 'surface_metrics/', 
                                     outputdir = "sudmex_conn_surf_data_scm/", 
                                     measure = 'thickness', 
                                     roilabel = 'allaseg', 
                                     template='fsaverage',
                                     
                                     subj_ID = TRUE)
```

As an example, we provide the “allaseg” thickness for the same dataset
and run the same analysis:

``` r
allaseg_thickness = readRDS(file=paste0(demodata,"/SUDMEX_CONN_allaseg_thickness.rds"))
allaseg_model=RFT_vertex_analysis(
  model = as.factor(beh_data$group),
  contrast = as.factor(beh_data$group),
  surf_data = allaseg_thickness)
```

``` r
allaseg_model$cluster_level_results
```

    ## $`Positive contrast`
    ##   clusid nverts      P     X     Y     Z tstat                  region
    ## 1      1    363 <0.001 121.1 137.6 125.2  3.62         right-ventraldc
    ## 2      2    157 <0.001 144.8 136.5 137.8  3.31          left-ventraldc
    ## 3      3    312  0.002 113.0 135.4 155.4  2.76              brain-stem
    ## 4      4    176  0.005 147.9 131.4 158.3  3.14          left-ventraldc
    ## 5      5    111  0.006  92.6 137.9 136.2  2.89           right-putamen
    ## 6      6    160  0.008 107.3 131.5 159.0  3.18         right-ventraldc
    ## 7      7   1344  0.015 150.0 179.1 167.6  2.91  left-cerebellum-cortex
    ## 8      8     65  0.015 100.6 133.1 136.3  3.10           right-putamen
    ## 9      9    620  0.028  99.8 170.5 167.9  2.94 right-cerebellum-cortex
    ## 
    ## $`Negative contrast`
    ##    clusid nverts      P     X     Y     Z tstat                  region
    ## 1       1   2062 <0.001 128.2 134.8 177.7 -3.22  left-cerebellum-cortex
    ## 2       2    223 <0.001 125.7 135.9 136.1 -4.36         right-ventraldc
    ## 3       3   1757 <0.001 104.5 148.0 171.7 -3.10 right-cerebellum-cortex
    ## 4       4    213 <0.001 106.1 132.1 163.4 -3.05          right-thalamus
    ## 5       5    186  0.001 153.8 128.7 163.5 -3.62           left-thalamus
    ## 6       6     77  0.001 112.9 133.0 116.7 -3.44           right-caudate
    ## 7       7     99  0.005 151.4 134.2 150.0 -3.44          left-ventraldc
    ## 8       8    162  0.007 132.3 133.0 148.4 -3.55           left-thalamus
    ## 9       9    120  0.008 150.5 160.0 186.3 -3.58  left-cerebellum-cortex
    ## 10     10   1605  0.014 171.1 174.7 204.5 -4.07  left-cerebellum-cortex
    ## 11     11     92  0.014 100.9 136.2 145.3 -3.73         right-ventraldc
    ## 12     12    127  0.016 115.3 118.5 161.0 -3.35          right-thalamus
    ## 13     13    171  0.022 136.0 112.5 141.9 -3.20           left-thalamus
    ## 14     14     74  0.029 107.3 120.9 116.1 -2.93           right-caudate
    ## 15     15     96  0.031 106.3 147.3 151.8 -3.38       right-hippocampus
    ## 16     16     68  0.035 106.2 134.8 141.4 -3.26         right-ventraldc

We can plot the subcortices together the same way:

``` r
allasegplot=plot_surf(
  surf_data=allaseg_model$thresholded_tstat_map,
  filename='SUDMEX_CONN_allaseg_tstatmaps.png',
  title='Thickness',
  smooth_mesh=20,
  show.plot.window=TRUE)
```

![Significant clusters after RFT correction across all ASeg subcortices,
2D plot](VertexWiseR_Example_3_files/figure-html/unnamed-chunk-13-1.png)

Since some regions may be harder to see in the 2d plot, plot_surf3d()
function has a special interactive slider for the “all merged” surface
in the GUI that allows spacing out each regions from each other (this
only applies to the merged surfaces):

``` r
plot_surf3d(surf_data = allaseg_model$thresholded_tstat_map,
            surf_color = 'grey', smooth_mesh = 50)
```

![Significant clusters after RFT correction across all ASeg subcortices,
3D plot (screenshot)](allaseg_plotsurf3d.jpg)

(This is just a screenshot)

Although SubCortexMesh requires that all regions be processed for the
encompassing surface to be created, it is possible to omit specific
regions from the analysis by turning to NA their corresponding vertices.

Users can identify the exact vertex indices in the python module’s
template_data. Here is the table that applies for ASeg:

|     |                         |        |              |
|-----|-------------------------|--------|--------------|
| id  | label                   | n_vert | n_vert_cumul |
| 0   | brain-stem              | 9452   | 9452         |
| 1   | left-accumbens-area     | 1022   | 10474        |
| 2   | left-amygdala           | 1638   | 12112        |
| 3   | left-caudate            | 3440   | 15552        |
| 4   | left-cerebellum-cortex  | 19550  | 35102        |
| 5   | left-hippocampus        | 4046   | 39148        |
| 6   | left-pallidum           | 1600   | 40748        |
| 7   | left-putamen            | 4268   | 45016        |
| 8   | left-thalamus           | 3936   | 48952        |
| 9   | left-ventraldc          | 3550   | 52502        |
| 10  | right-accumbens-area    | 1022   | 53524        |
| 11  | right-amygdala          | 1792   | 55316        |
| 12  | right-caudate           | 3500   | 58816        |
| 13  | right-cerebellum-cortex | 19664  | 78480        |
| 14  | right-hippocampus       | 4086   | 82566        |
| 15  | right-pallidum          | 1600   | 84166        |
| 16  | right-putamen           | 4126   | 88292        |
| 17  | right-thalamus          | 3832   | 92124        |
| 18  | right-ventraldc         | 3594   | 95718        |

So, for example, if one does not want to count the cerebella in the
analysis, they can follow the cumulative vertex count and set them to
NA, like that:

``` r
allaseg_thickness$surf_obj[,c((15552+1):35102,(58816+1):78480)]=NA

allaseg_model_nocerebellum=RFT_vertex_analysis(
  model = as.factor(beh_data$group),
  contrast = as.factor(beh_data$group),
  surf_data = allaseg_thickness)
```

``` r
allaseg_model_nocerebellum$cluster_level_results
```

    ## $`Positive contrast`
    ##   clusid nverts      P     X     Y     Z tstat          region
    ## 1      1    363 <0.001 121.1 137.6 125.2  3.62 right-ventraldc
    ## 2      2    157 <0.001 144.8 136.5 137.8  3.31  left-ventraldc
    ## 3      3    312  0.002 113.0 135.4 155.4  2.76      brain-stem
    ## 4      4    176  0.004 147.9 131.4 158.3  3.14  left-ventraldc
    ## 5      5    111  0.005  92.6 137.9 136.2  2.89   right-putamen
    ## 6      6    160  0.007 107.3 131.5 159.0  3.18 right-ventraldc
    ## 7      7     65  0.013 100.6 133.1 136.3  3.10   right-putamen
    ## 
    ## $`Negative contrast`
    ##    clusid nverts      P     X     Y     Z tstat            region
    ## 1       1    223 <0.001 125.7 135.9 136.1 -4.36   right-ventraldc
    ## 2       2    213 <0.001 106.1 132.1 163.4 -3.05    right-thalamus
    ## 3       3    186  0.001 153.8 128.7 163.5 -3.62     left-thalamus
    ## 4       4     77  0.001 112.9 133.0 116.7 -3.44     right-caudate
    ## 5       5     99  0.004 151.4 134.2 150.0 -3.44    left-ventraldc
    ## 6       6    162  0.006 132.3 133.0 148.4 -3.55     left-thalamus
    ## 7       7     92  0.012 100.9 136.2 145.3 -3.73   right-ventraldc
    ## 8       8    127  0.013 115.3 118.5 161.0 -3.35    right-thalamus
    ## 9       9    171  0.019 136.0 112.5 141.9 -3.20     left-thalamus
    ## 10     10     74  0.025 107.3 120.9 116.1 -2.93     right-caudate
    ## 11     11     96  0.027 106.3 147.3 151.8 -3.38 right-hippocampus
    ## 12     12     68   0.03 106.2 134.8 141.4 -3.26   right-ventraldc
    ## 13     13    105  0.046 148.3 147.9 152.0 -3.17  left-hippocampus

## References:

Fischl, Bruce. 2012. “FreeSurfer.” *NeuroImage* 62 (2): 774–81.
<https://doi.org/10.1016/j.neuroimage.2012.01.021>.

Garza-Villarreal, EA, MM Chakravarty, B Hansen, SF Eskildsen, GA
Devenyi, D Castillo-Padilla, T Balducci, et al. 2017. “The Effect of
Crack Cocaine Addiction and Age on the Microstructure and Morphology of
the Human Striatum and Thalamus Using Shape Analysis and Fast Diffusion
Kurtosis Imaging.” *Translational Psychiatry* 7 (5): e1122.
<https://doi.org/10.1038/tp.2017.92>.

Patenaude, Brian, Stephen M. Smith, David N. Kennedy, and Mark
Jenkinson. 2011. “A Bayesian Model of Shape and Appearance for
Subcortical Brain Segmentation.” *NeuroImage* 56 (3): 907–22.
<https://doi.org/10.1016/j.neuroimage.2011.02.046>.

Xu, Hui, Cheng Xu, and Chenguang Guo. 2023. “Cocaine Use Disorder Is
Associated with Widespread Surface-Based Alterations of the Basal
Ganglia.” *Journal of Psychiatric Research* 158: 95–103.
