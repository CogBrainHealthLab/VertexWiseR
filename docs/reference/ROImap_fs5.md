# Atlas parcellations of fsaverage5

A list containing two data frames, 1) listing fsaverage5 vertices and
each parcellation number they correspond to, and 2) listing each
available atlas and their corresponding labels (1=aparc,
2=Destrieux-148, 3=Glasser-360, 4=Schaefer-100, 5=Schaefer-200,
6=Schaefer-400).

## Usage

``` r
ROImap_fs5
```

## Format

### `ROImap_fs5`

A list object with two data frame objects: ()

- vertices:

  data frame with 20484 rows (vertices), 6 columns (atlases)

- atlases:

  data frame with 400 rows (labels, not all are filled depending on
  atlas), 6 columns (atlases)
