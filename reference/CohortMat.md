# Cohort Matrix Class

Initialize the CohortMat class

## Usage

``` r
CohortMat(
  obs_df = NULL,
  dbh = NULL,
  trees = NULL,
  species = NULL,
  dims = c(50, 30, 10),
  sp = 10,
  device = "cpu"
)
```

## Arguments

- obs_df:

  A data frame containing columns "siteID", "patchID", "species", "dbh",
  and "trees". If provided, it will be used to initialize the tensors.

- dbh:

  A tensor or array representing the diameter at breast height. Defaults
  to `self$dbh`.

- trees:

  A tensor or array representing the number of trees. Defaults to
  `self$trees`.

- species:

  A tensor or array representing the species. Defaults to
  `self$species`.

- dims:

  A numeric vector representing the dimensions of the arrays (sites,
  patches, cohorts). Defaults to `self$dims`.

- sp:

  An integer representing the number of species. Defaults to `self$sp`.

- device:

  A character string specifying the device to use ('cpu' or 'cuda').
  Defaults to `self$device`. @description

  Convert the tensors to a data frame format

## Value

data.table object

## Details

An R6 class for managing cohorts of trees in forest models. This class
allows for the initialization, transformation, and manipulation of
cohorts represented by arrays of dbh, trees, and species.

## Fields

- `dbh`:

  A tensor representing the diameter at breast height for each cohort.

- `trees`:

  A tensor representing the number of trees in each cohort.

- `species`:

  A tensor representing the species of each cohort.

- `dims`:

  A vector representing the dimensions of the arrays (sites, patches,
  cohorts).

- `sp`:

  An integer representing the number of species.

- `device`:

  A character string specifying the device to use ('cpu' or 'cuda').

- `dbh_r`:

  A numeric array representing the diameter at breast height in R array
  format.

- `trees_r`:

  A numeric array representing the number of trees in R array format.

- `species_r`:

  An integer array representing the species in R array format.

- `device_r`:

  A character string specifying the device in R format ('cpu' or
  'cuda').

- `obsDF2arrays`:

  Transforms data.table into array
