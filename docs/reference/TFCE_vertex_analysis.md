# Vertex-wise analysis with threshold-free cluster enhancement (fixed effect)

Fits a linear model with the cortical or hippocampal surface data as the
predicted outcome, and returns t-stat and threshold-free cluster
enhancement (TFCE) statistical maps for the selected contrast.

## Usage

``` r
TFCE_vertex_analysis(
  model,
  contrast,
  formula,
  formula_dataset,
  inverse = FALSE,
  surf_data,
  nperm = 100,
  tail = 2,
  nthread = 10,
  smooth_FWHM,
  VWR_check = TRUE
)
```

## Arguments

- model:

  An N X P data.frame object containing N rows for each subject and P
  columns for each predictor included in the model

- contrast:

  A N x 1 numeric vector or object containing the values of the
  predictor of interest. Its length should equal the number of subjects
  in model (and can be a single column from model). The t-stat and TFCE
  maps will be estimated only for this predictor.

- formula:

  An optional string or formula object describing the predictors to be
  fitted against the surface data, replacing the model, contrast, or
  random arguments. If this argument is used, the formula_dataset
  argument must also be provided.

  - The dependent variable (DV) is not needed, and the formula will
    start with ~. The DV will be the surface data value by default, but
    it can be swapped with contrast as IV via the "inverse" argument.

  - The first independent variable in the formula will always be
    interpreted as the contrast of interest for which to estimate
    cluster-thresholded t-stat maps.

  - Only one random regressor can be given and must be indicated as
    '(1\|variable_name)'.

- formula_dataset:

  An optional data.frame object containing the independent variables to
  be used with the formula (the IV names in the formula must match their
  column names in the dataset).

- inverse:

  A boolean object stating whether to set the surface data as predictor
  of the contrast variable, instead of as dependent variable (default is
  FALSE). Other covariates in the model remain independent variables.
  This makes modelling slower.

- surf_data:

  A N x V matrix object containing the surface data (N row for each
  subject, V for each vertex), in fsaverage5 (20484 vertices),
  fsaverage6 (81924 vertices), fslr32k (64984 vertices) or hippocampal
  (14524 vertices) space. See also Hipvextract(), SURFvextract(),
  FSLRvextract(), or SCMvextract() output formats. Alternatively, a
  string object containing the path to the surface object (.rds file)
  outputted by extraction functions may be given.

- nperm:

  An integer specifying the number of permutations generated for the
  subsequent thresholding procedures (default = 100)

- tail:

  An integer specifying whether to test a one-sided positive (1),
  one-sided negative (-1) or two-sided (2) hypothesis

- nthread:

  An integer specifying the number of CPU threads to allocate

- smooth_FWHM:

  A numeric vector specifying the desired smoothing width in mm. It
  should not be specified if the surf_data has been smoothed previously
  with smooth_surf(), because this results in surf_data being smoothed
  twice.

- VWR_check:

  A boolean object specifying whether to check and validate system
  requirements. Default is TRUE.

## Value

A list object containing the t-stat and the TFCE statistical maps which
can then be subsequently thresholded using TFCE_threshold()

## Details

This implements the Freedman–Lane nuisance-adjustment strategy used by
FSL's `randomise`, while adapting the TFCE method from the ['Nilearn'
Python
library](https://github.com/nilearn/nilearn/blob/main/nilearn/mass_univariate/_utils.py#L7C8-L7C8).

A reduced linear model is fitted using an intercept and all nuisance
predictors, excluding the `contrast`. For each permutation, residuals
from this reduced model are shuffled across observations and added back
to its fitted values. The full linear model is then refitted to obtain
the t-statistic for the tested predictor at each vertex.

## See also

[`RFT_vertex_analysis`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/RFT_vertex_analysis.md),
[`TFCE_vertex_analysis_mixed`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/TFCE_vertex_analysis_mixed.md),
[`TFCE_threshold`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/TFCE_threshold.md)

## Examples

``` r
demodata = readRDS(system.file('demo_data/SPRENG_behdata_site1.rds',
package = 'VertexWiseR'))[1:5,]
CTv = readRDS(file = url(paste0("https://github.com",
"/CogBrainHealthLab/VertexWiseR/blob/main/inst/demo_data/",
"SPRENG_CTv_site1.rds?raw=TRUE")))[1:5,]

TFCEpos=TFCE_vertex_analysis(model=demodata[,c("sex","age")], 
contrast=demodata[,"age"], surf_data=CTv, tail=1, 
nperm=5, nthread = 2, VWR_check=FALSE)
#> The binary variable 'sex' will be recoded with F=0 and M=1 for the analysis
#> smooth_FWHM argument was not given. surf_data will not be smoothed here.
#> Estimating unpermuted TFCE image...
#> Completed in 5.5 secs
#> Estimating permuted TFCE images...
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#> Error in {    if (inverse == FALSE) {        mod.permuted = .lm.fit(y = fitted0 + surf_data[permseq[,             perm], , drop = FALSE], x = data.matrix(cbind(1,             model)))        tmap = extract_t(mod.permuted, colno + 1)    }    else {        tmap = c()        for (vert in 1:ncol(surf_data)) {            vertmodel.permuted = cbind(invmodel, surf_data[,                 vert])            invmod.permuted = .lm.fit(y = fitted0 + contrast[permseq[,                 perm]], x = data.matrix(cbind(1, vertmodel.permuted)))            tmap = c(tmap, extract_t(invmod.permuted, ncol(vertmodel.permuted) +                 1))        }    }    . <- model.permuted <- NULL    remove(mod.permuted, model.permuted)    return(max(abs(suppressWarnings(TFCE(data = tmap, tail = tail,         edgelist = edgelist)))))}: task 2 failed - "could not find function "extract_t""

#To threshold the results, you may then run:
#results=TFCE_threshold(TFCEpos, p=0.05, atlas=1)
#results$cluster_level_results

#Formula alternative:
#formula= as.formula("~ age + sex")
#TFCEpos=TFCE_vertex_analysis(formula=formula, 
#formula_dataset=demodata, surf_data=CTv, tail=1, 
#nperm=5, nthread = 2, VWR_check=FALSE) 
```
