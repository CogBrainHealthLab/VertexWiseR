# Region-of-Interest mapping object

A class for surface vertices mapping on atlas labels

## Slots

- `data`:

  A matrix object with N vertices from a template and each parcellation
  number the vertices correspond to in 6 atlases (6 columns).

- `atlases`:

  Each available of the 6 available atlases and their corresponding
  labels (1=aparc, 2=Destrieux-148, 3=Glasser-360, 4=Schaefer-100,
  5=Schaefer-200, 6=Schaefer-400).

- `name`:

  The name of the template surface
