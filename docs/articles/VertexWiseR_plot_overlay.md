# Overlaying plots and transparent thresholding

## Transparent thresholding

This article presents the plot_overlay_surf() function which can be used
to overlay plots similar to those produced with surf_plot() in
VertexWiseR. We show how we can use the function to apply “transparent
thresholding”, which consists in reporting spatial neuroimaging results
at a lower, non-significant subthreshold, to contextualise significant
findings and improve broader interpretation across studies (Taylor et
al. 2025).

The demo data (~216 MB) required to run the plotting examples, and can
be downloaded from the package’s github repository with the following
function:

``` r
#This will save the demo_data directory in a temporary directory (tempdir(), but you can change it to your own path)
demodata=VertexWiseR:::dl_demo(path=tempdir(), quiet=TRUE)
```

## Running a vertex-wise linear model

In this example, we use the Spreng dataset cortical thickness and
behavioural data as in [Example
1](https://cogbrainhealthlab.github.io/VertexWiseR/articles/VertexWiseR_Example_1.html)
and the [home
example](https://cogbrainhealthlab.github.io/VertexWiseR/index.html) ,
selecting young adults:

``` r
SPRENG_CTv = readRDS(file=paste0(demodata,"/SPRENG_CTv_site1.rds"))
dat_beh = readRDS(paste0(demodata,'/SPRENG_behdata_site1.rds'))
#only young adults
SPRENG_CTv = SPRENG_CTv[dat_beh$agegroup=='Y',]
dat_beh = dat_beh[dat_beh$agegroup=='Y',]
#smoothing at 10 FMHW
SPRENG_CTv_smoothed = smooth_surf(SPRENG_CTv, 10, VWR_check = T)
```

We test the effect of sex on cortical thickness in the sample:

``` r
#preparing model data
predictors=dat_beh[,c("site","age","sex")]
#running the model
results=RFT_vertex_analysis(model=predictors, 
                            contrast=predictors$sex, 
                            surf_data = SPRENG_CTv_smoothed,
                            p = 0.05 
                            )
results$cluster_level_results
```

    ## $`Positive contrast`
    ##   clusid nverts      P     X     Y     Z tstat              region
    ## 1      1    301 <0.001 -40.7 -13.0  16.5  5.06      lh-postcentral
    ## 2      2    181  0.002  56.3  12.1 -10.1  4.43 rh-superiortemporal
    ## 3      3    150  0.003  21.9 -51.6  -1.2  4.41          rh-lingual
    ## 
    ## $`Negative contrast`
    ## [1] "No significant clusters"

According to these results, since females are coded as 0 and males as 1,
the regions colored in red are thicker in males.

## Plotting significant and subthreshold outcome

Both RFT_vertex_analysis() and TFCE_threshold() return thresholded and
unthresholded t-maps, and can be plotted individually as follows.

Unthresholded t-map:

``` r
plot_surf(surf_data = results$tstat_map,
          filename = "sexdiff_nothresh.png",
          show.plot.window=TRUE)
```

![Unthresholded t-statistic
map](VertexWiseR_plot_overlay_files/figure-html/unnamed-chunk-4-1.png)

Thresholded (significant) t-map:

``` r
plot_surf(surf_data = results$thresholded_tstat_map,
          filename = "sexdiff_thresh.png",
          show.plot.window=TRUE)
```

![Thresholded t-statistic
map](VertexWiseR_plot_overlay_files/figure-html/unnamed-chunk-5-1.png)

## Overlaying plots together

To apply transparency thresholding, and letting information about
subthreshold effects present in the plot, one may want to “merge” both
maps together. plot_overlay() is a function that will allow you to plot
an overlay surface (surf_data_2) on top of a background (surf_data_1).

To work in this context, we plot a thresholded (significant) t-stat map
on top of an unthresholded t-stat map, both of which can be
automatically taken from the ‘results’ variable saved from the
RFT_vertex_analysis() model.

We reduce the alpha opacity of the background surface to highlight the
surface on top:

``` r
plot_overlay_surf(model_output=results,
                  #To specify maps manually:
                  #surf_data_1=results$tstat_map,
                  #surf_data_2=results$thresholded_tstat_map, 
                  cmap_1='RdBu_r', cmap_2='RdBu_r',
                  colorbar_1=FALSE, colorbar_2=TRUE,
                  alpha_1=0.4, alpha_2=1,
                  overlay_boundaries=TRUE,
                  limits_2='same',
                  filename='overlay_plot.png',
                  size=c(1400,291), 
                  title="Significant effect\nof sex",
                  show.plot.window=TRUE)
```

![Overlaid plots with the same color
map](VertexWiseR_plot_overlay_files/figure-html/unnamed-chunk-6-1.png)

By default, both color bars are plotted but either can be set to FALSE
(removed), which was done for the faded background colours here.

The function also give freedom to have two separate color maps and
independent ranges of values (which is the default if the limits==‘same’
option is not provided):

``` r
plot_overlay_surf(model_output=results, 
                  cmap_1='RdBu_r', cmap_2='hot_r',
                  alpha_1=1, alpha_2=1,
                  filename='overlay_plot_twolimits.png',
                  size=c(1400,291),
                  title="Significant effect\nof sex",
                  show.plot.window=TRUE)
```

![Overlaid plots with different color maps and
limits](VertexWiseR_plot_overlay_files/figure-html/unnamed-chunk-7-1.png)

## Overlaying HippUnfold hippocampal plots

The same can be applied to HippUnfold hippocampal surfaces. We reproduce
the results from the first analysis in [Example
2](https://cogbrainhealthlab.github.io/VertexWiseR/articles/VertexWiseR_Example_2.html):

``` r
FINK_Tv_ses13 = readRDS(file=paste0(demodata,"/FINK_Tv_ses13.rds"))
dat_beh_ses13 = readRDS(paste0(demodata,"/FINK_behdata_ses13.rds"))
model2_RFT=RFT_vertex_analysis(
  model = dat_beh_ses13[,c("session","group","session_x_group")],
  contrast = dat_beh_ses13[,"session_x_group"],
  surf_data=FINK_Tv_ses13,
  random=dat_beh_ses13[,"participant_id"], 
  smooth_FWHM = 5,
  p=0.05)
```

Plotted in the same manner:

``` r
plot_overlay_surf(model_output=model2_RFT, 
                  cmap_1='RdBu_r', cmap_2='RdBu_r',
                  colorbar_1=FALSE, colorbar_2=TRUE,
                  alpha_1=0.4, alpha_2=1,
                  limits_2='same',
                  filename='overlay_plot_hippocampus.png',
                  overlay_boundaries=TRUE,
                  size=c(1400,300),
                  title="Effect of session\nand group",
                  show.plot.window=TRUE)
```

![Overlaid hippocampal
plots](VertexWiseR_plot_overlay_files/figure-html/unnamed-chunk-9-1.png)

## Overlaying SubCortexMesh subcortical plots

Likewise subcortices obtained via ASEGvextract() can also be overlaid,
but that requires downloading associated template data (the package will
prompt you automatically or you may do it directly by typing
VertexWiseR:::scm_database_check()).

We reproduce the results from the example analysis in [the dedicated
vignette](https://cogbrainhealthlab.github.io/VertexWiseR/articles/VertexWiseR_Example_3.html):

``` r
beh_data = readRDS(paste0(demodata,"/SUDMEX_CONN_behdata.rds"))
thalamus_thickness = readRDS(file=paste0(demodata,"/SUDMEX_CONN_thalamus_thickness.rds"))
model_thalamus_thickness=RFT_vertex_analysis(
  model = as.factor(beh_data$group),
  contrast = as.factor(beh_data$group),
  surf_data = thalamus_thickness)
```

Plotted in the same manner:

``` r
plot_overlay_surf(model_thalamus_thickness, 
                  colorbar_1=F, 
                  filename = 'overlay_plot_thalamus.png',
                  cmap_1 = c('darkblue','blue'), 
                  cmap_2=c('yellow','white'), 
                  overlay_boundaries = FALSE,
                  show.plot.window = TRUE
                  )
```

![Overlaid thalamic
plots](VertexWiseR_plot_overlay_files/figure-html/unnamed-chunk-11-1.png)

## References:

Taylor, Paul A, Himanshu Aggarwal, Peter Bandettini, Marco Barilari,
Molly Bright, Cesar Caballero-Gaudes, Vince Calhoun, et al. 2025. “Go
Figure: Transparency in Neuroscience Images Preserves Context and
Clarifies Interpretation.” *arXiv Preprint arXiv:2504.07824*.
