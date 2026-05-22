# VertexWiseR system requirements installation

Helps the user verify if VertexWisrR's system requirements are present
and install them (a suitable 'Python' or 'Miniconda' environment,
'BrainStat' toolbox and libraries). If they are installed already,
nothing will be overwritten.

## Usage

``` r
VWRfirstrun(requirement = "any", n_vert = 0, promptless = FALSE)
```

## Arguments

- requirement:

  String that specifies a requirement to enquire about:

  - For only Python/Conda installation: 'python/conda only'

  - For Python/Conda and Brainstat installation: 'conda/brainstat'

  - For specific 'BrainStat' libraries: 'fsaverage5', 'fsaverage6',
    'fslr32k', 'yeo_parcels'

  - For the neurosynth database: 'neurosynth'. Default is 'any' and
    checks everything.

- n_vert:

  Numeric vector indicating the number of vertices of a given surface
  data so that only the required templates are asked for. It will modify
  the requirement argument accordingly.

- promptless:

  A boolean object specifying whether to prompt the user for action when
  system requirements are missing. If TRUE, VWRfirstrun() will simply
  inform of what is missing and will not prompt for action. Default is
  FALSE.

## Value

No returned value in interactive session. In non-interactive sessions, a
string object informing that system requirements are missing.

## Details

VertexWiseR imports and makes use of the R package 'reticulate.'
'reticulate' is a package that allows R to borrow or translate Python
functions into R. Using 'reticulate', the package calls functions from
the 'BrainStat' Python module. For 'reticulate' to work properly with
VertexWiseR, a Python environment needs to be installed with it — the
default choice offered by VWRfirstrun is to let reticulate (version
1.41.0) create an ephemeral Python virtual environment using UV and
py_require(). If for a reason this is not desirable, VWRfirstrun() also
gives the choice to install a 'Miniconda' (lightweight version of
Python) or Python environment in a reticulate default path or a
specified path. Vertex-wise statistical analyses of cortical surface
require fsaverage and parcellation templates as imported by default in
'BrainStat'. The decode_surf_data() function also requires the
'Neurosynth' database to be downloaded.

## Examples

``` r
VWRfirstrun()
```
