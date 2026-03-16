# Python troubleshooting

VertexWiseR uses Python functions from various toolboxes that can
flexibly be called and executed in R using the *reticulate* package
(Ushey, Allaire, and Tang 2023). As of version \< 1.41.0 reticulate
favours creating a virtual environment, where the python version and
specific packages needed for VertexWiseR can be stored, leaving any
other system installations intact. This is also the option we recommend.
The [reticulate
documentation](https://rstudio.github.io/reticulate/articles/python_packages.html)
explains:

> When installing Python packages it’s best practice to isolate them
> within a Python environment (a named Python installation that exists
> for a specific project or purpose). This provides a measure of
> isolation, so that updating a Python package for one project doesn’t
> impact other projects. The risk for package incompatibilities is
> significantly higher with Python packages than it is with R packages,
> because unlike CRAN, PyPI does not enforce, or even check, if the
> current versions of packages currently available are compatible.

More details on why using virtual environments can be found in the
[reticulate
documentation](https://rstudio.github.io/reticulate/articles/python_packages.html).

Note that VertexWiseR will install a specific version of Python itself
and of each package it needs to work properly, whether in the virtual
environment or a stable environment VWRfirstrun(). If users struggle to
create a virtual environment, they may opt to install a (non-virtual)
Miniconda/Python environment in a path of their choice. Those options
are given when running VWRfirstrun().

We explain in this article how to solve a number of issues that may
arise from the interaction with Python:

1.  [How to make VWRfirstrun() ignore the system’s Python environment
    when starting a fresh installation](#id1)

2.  [How to make VWRfirstrun() use a specific Python installation](#id2)

    1.  [For a new, fresh VertexWiseR set-up](#id2.a)

    2.  [If the VWRfirstrun() installations were already done
        before](#id2.b)

3.  [VWRfirstrun() shows the message “Clearing reticulate’s uv
    cache…”](#id3)

4.  [About the message: “Would you like to create a default Python
    environment for the reticulate package? (Yes/no/cancel)” (reticulate
    \< 1.41.0)](#id4)

------------------------------------------------------------------------

## 1. How to make VWRfirstrun() ignore the system’s Python environment when starting a fresh installation

VertexWiseR will automatically select and run a Python installation it
finds via reticulate, but users may want to prevent this and run
VertexWiseR on a new Python installation of its own. This section
explains how to make
[`VWRfirstrun()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/VWRfirstrun.md)
ignore a previous Python/Miniconda installation in order to install a
new one (which will *not* affect the system’s Python libraries).

Firstly, the Python installation path may already be preloaded in the
global environment, via files such as `.Renviron`, `.Rprofile`,
`.bashrc` etc. If that is the case they should be edited to remove the
path of a previous Python installation, which will avoid reticulate from
using those paths automatically.

Secondly, if the Python paths are not specified in the global
environment, reticulate might search for a Python installation in its
default locations. This can be avoided with the method below.

To make
[`VWRfirstrun()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/VWRfirstrun.md)
ignore a previous Python installation, users can set the
`RETICULATE_PYTHON` variable as `NA` (in a new R session with a clean
workspace) before running the function:

``` r
Sys.setenv(RETICULATE_PYTHON=NA)
VWRfirstrun()
```

With the above lines, reticulate should fail to load Python and respond
in the same way it would if Python was not previously installed at all.
Therefore,
[`VWRfirstrun()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/VWRfirstrun.md)
will prompt users to install a new Python environment. After that,
[`VWRfirstrun()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/VWRfirstrun.md)
will use that newly installed Python environment by default and users
should not need to do this again.

## 2. How to make VWRfirstrun() use a specific Python installation

This section explains how users can choose a different Python
installation - other than the one automatically detected - if they wish
to, instead of installing a new one.

### a. For a new, fresh VertexWiseR set-up

Reticulate allows users to predefine the Python environment to be used
by entering their Python’s directory path into the [following
functions](https://rstudio.github.io/reticulate/reference/use_python.html):

``` r
#replace each word between brackets with the relevant path
reticulate::use_python(python)
reticulate::use_virtualenv(virtualenv)
reticulate::use_condaenv(condaenv)
reticulate::use_miniconda(condaenv)
```

These functions must be executed before
[`VWRfirstrun()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/VWRfirstrun.md)
is run in the R session. The latter function should then automatically
assume that a version of Python is correctly installed and only ask for
other packages to be downloaded/installed if they are missing. To check
the Python installation that is being used, users can run:

``` r
reticulate::py_config()
```

This will print the information and path of the Python environment
reticulate is currently using. If that is still not the path specified
in the functions listed above, users should try restarting R and/or
clearing the R workspace first.

### b. If the VWRfirstrun() installations were already done before

If Python or Miniconda was already installed using
[`VWRfirstrun()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/VWRfirstrun.md),
users will not be able to modify the Python path in the above manner.
That is because the path to the user-selected default or custom
installations is saved in order to conveniently access them again later.
Every time
[`VWRfirstrun()`](https://cogbrainhealthlab.github.io/VertexWiseR/reference/VWRfirstrun.md)
is run, it will read and use such saved paths. This also depends on
whether the user created a default ephemeral environment or chose a
manual Python/Miniconda installation.

- **If the user installed a Miniconda or Python environment in the
  classic way**

The paths to these installations are saved inside the `.Renviron` file,
in the standard path generated via:

``` r
tools::R_user_dir(package='VertexWiseR')
```

Within the `.Renviron` file, the environment variables will be written
in this manner:

``` r
RETICULATE_MINICONDA_PATH="C:/path_to_miniconda/"
RETICULATE_PYTHON_FALLBACK="C:/path_to_miniconda/"
RETICULATE_PYTHON="C:/path_to_miniconda/python.exe"
```

To swap with another Python environment, users may choose to edit
`.Renviron` directly, by replacing the paths with the Python
installation they would like VertexWiseR to use instead. The
`RETICULATE_MINICONDA_PATH` may also be safely removed if a non-conda
Python library is to be used. Alternatively, users may choose to
completely remove the `.Renviron` and redefine the path as described in
section 3a.

- **If the user created an ephemeral virtual environment**

This is in the case that users chose the virtual environment
installation proposed by VWRfirstrun(), which calls
reticulate::py_require(). The path to such environment is accessible
with that command line:

``` r
tools::R_user_dir("reticulate", "cache")
```

To stop using this environment, users may choose to clear the cache
completely:

``` r
unlink(tools::R_user_dir("reticulate", "cache"), recursive = TRUE)
```

If they do not wish to get rid of the previous virtual environment,
users will have to modify and/or remove the VIRTUAL_ENV variable which
has been stored in the `.Renviron` file, in the standard path generated
via:

``` r
tools::R_user_dir(package='VertexWiseR')
```

This variable is automatically removed when the virtual environment can
no longer be found in the cache, but will automatically be used by
VWRfirstrun() after its installation to speed up the process, if it
still exists. So the VIRTUAL_ENV variable must be rewritten to the new
ephemeral environment generated by py_require() to use the latter.

## 3. VWRfirstrun() shows the message “Clearing reticulate’s uv cache…”

That is the [default behaviour of
reticulate](https://pkgs.rstudio.com/reticulate/reference/py_require.html).
The default duration of a virtual environment is 120 days, but one may
choose to extend by setting this option in \[.Rprofile\]
(<https://rstats.wtf/r-startup.html#rprofile>) manually:

`options(reticulate.max_cache_age = as.difftime(30, units = "days"))`

In the event that a cache clearing is triggered, the message “Clearing
reticulate’s uv cache…” will appear and the VertexWiseR function that
triggered it will most likely fail. Simply restart the R session and run
VWRfirstrun() again to install a fresh virtual python environment.

## 4. About the message: “Would you like to create a default Python environment for the reticulate package? (Yes/no/cancel)” (reticulate \< 1.41.0)

As explained in the [home
page](https://cogbrainhealthlab.github.io/VertexWiseR/) footnote, this
is a pop-up or prompt from the reticulate package. Choosing ‘Yes’ will
let reticulate install Python in a virtual environment. This is entirely
optional. Simply choose/click No/Cancel to ignore.

## References:

Ushey, K, J Allaire, and Y Tang. 2023. “Reticulate: Interface to
’Python’.” <https://CRAN.R-project.org/package=reticulate>.
