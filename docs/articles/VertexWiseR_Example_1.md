# Example analyses with VertexWiseR - Example 1

## Installing VertexWiseR and checking for requirements

The following code installs the package and makes sure all requirements,
prompting the user to install dependencies, in order to allow analyses
to work.

``` r
install.packages("VertexWiseR") 
```

``` r
library(VertexWiseR)
```

VWRfirstrun() checks all system requiremements for specific functions,
and gives the opportunity to download and install each of them.

``` r
VWRfirstrun()
```

## Example analysis 1: linear model of age and cortical thickness, and meta-analytic decoding

The analysis will use surface data already extracted in R from a
preprocessed subjects directory, which we make available so you do not
need to preprocess a sample yourself. To obtain it, we had extracted
cortical thickness (CT) data from a FreeSurfer preprocessing directory
of the [SPRENG
dataset](https://openneuro.org/datasets/ds003592/versions/1.0.13)
(Spreng et al. 2022).

The demo data (~216 MB) required to run the example analysis, and can be
downloaded from the package’s github repository with the following
function:

``` r
#This will save the demo_data directory in a temporary directory (tempdir(), but you can change it to your own path)
demodata=VertexWiseR:::dl_demo(path=tempdir(), quiet=TRUE)
```

Here is the command line which was originally used to extract the
surface:

``` r
#SURFvextract(sdirpath = SUBJECTS_DIR, filename = "SPRENG_CTv", template='fsaverage5', measure = 'thickness', subj_ID = T)
```

To load the demo thickness matrix:

``` r
SPRENG_CTv = readRDS(file=paste0(demodata,"/SPRENG_CTv_site1.rds"))
```

To smooth the surface data:

``` r
SPRENG_CTv_smoothed = smooth_surf(SPRENG_CTv, 10)
```

To load the demo behavioural data (for participants in site 1,
SPRENG_behdata_site1.csv):

``` r
dat_beh=readRDS(paste0(demodata,'/SPRENG_behdata_site1.rds'))
```

To run the vertex-wise model analysis with random field theory-based
cluster correction, testing for the effect of age, controlling for sex,
on CT:

``` r
model1_RFT=RFT_vertex_analysis(model = dat_beh[,c("sex","age")], 
  contrast = dat_beh[,"age"], surf_data = SPRENG_CTv_smoothed, p = 0.05)
```

``` r
model1_RFT$cluster_level_results
```

    ## $`Positive contrast`
    ##   clusid nverts     P     X    Y   Z tstat          region
    ## 1      1    142 0.015 -22.8 11.5 -42  6.45 lh-temporalpole
    ## 
    ## $`Negative contrast`
    ##   clusid nverts      P   X     Y     Z  tstat              region
    ## 1      1   8039 <0.001  47   4.0 -16.6 -12.64 rh-superiortemporal
    ## 2      2   7660 <0.001 -34 -25.7  16.2 -14.23           lh-insula

To run the vertex-wise model analysis with threshold-free cluster
enhancement-based cluster correction, testing for the effect of age,
controlling for sex, on CT; with 1000 permutations:

``` r
model1_TFCE=TFCE_vertex_analysis(model= dat_beh[,c("sex","age")], 
                                 contrast = dat_beh[,"age"],
                                 surf_data=SPRENG_CTv_smoothed,
                                 nperm=1000, 
                                 nthread=4) 

TFCEoutput = TFCE_threshold(model1_TFCE, p=0.05)
```

``` r
TFCEoutput$cluster_level_results
```

    ## $`Positive contrast`
    ## [1] "No significant clusters"
    ## 
    ## $`Negative contrasts`
    ##   clusid nverts      P   X     Y     Z tstat              region
    ## 1      1   8098 <0.001  47   4.0 -16.6 12.64 rh-superiortemporal
    ## 2      2   7617 <0.001 -34 -25.7  16.2 14.23           lh-insula

To plot the results of both models on an inflated fsaverage5 surface:

``` r
tmaps = rbind(model1_RFT$thresholded_tstat_map, TFCEoutput$thresholded_tstat_map)
plot_surf(surf_data = tmaps, 
          filename ='SPRENG_tstatmaps.png', 
          size=c(1400,582),
          surface = 'inflated', 
          title=c("RFT-corrected\nclusters", "TFCE-corrected\nclusters"),
          cmap='RdBu_r',
          show.plot.window=TRUE)
```

![Significant clusters after RFT correction and TFCE
correction](VertexWiseR_Example_1_files/figure-html/unnamed-chunk-13-1.png)

To run meta-analytic decoding of the significant negative clusters (the
neurosynth dataset needs to be installed as VWRfirstrun() allows):

``` r
surf_decoding=decode_surf_data(TFCEoutput$thresholded_tstat_map, contrast="negative")
```

``` r
head(surf_decoding)
```

    ##        keyword     r
    ## 538  retrieval 0.065
    ## 202   episodic 0.059
    ## 348     memory 0.053
    ## 198 engagement 0.048
    ## 332 linguistic 0.048
    ## 439      older 0.047

## References:

Spreng, R. Nathan, Roni Setton, Udi Alter, Benjamin N. Cassidy, Bri
Darboh, Elizabeth DuPre, Karin Kantarovich, et al. 2022. “Neurocognitive
Aging Data Release with Behavioral, Structural and Multi-Echo Functional
MRI Measures.” *Scientific Data* 9 (1): 119.
<https://doi.org/10.1038/s41597-022-01231-7>.
