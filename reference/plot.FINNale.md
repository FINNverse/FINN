# Plot ALE curves of a FINN model

Plots the accumulated local effect (ALE) curves produced by
[`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md) as a grid:
one row per demographic process (growth, mortality, regeneration) and
one column per environmental predictor. When several species are
present, their curves are overlaid in the same panel as one coloured
line each.

## Usage

``` r
# S3 method for class 'FINNale'
plot(x, process = NULL, scale = FALSE, ...)
```

## Arguments

- x:

  an object of class `FINNale` (returned by
  [`ALE()`](https://finnverse.github.io/FINN/reference/ALE.md)).

- process:

  (`character`\|`NULL`)  
  restrict the plot to one process (`"growth"`, `"mortality"` or
  `"regeneration"`); `NULL` (default) plots all three.

- scale:

  (`logical(1)`)  
  if `TRUE`, divide each curve by the SD of its process x species rate
  so the effects are dimensionless and comparable (a shared y-axis is
  then used); if `FALSE` (default) the raw ALE is shown with free
  y-axes.

- ...:

  currently ignored.

## Value

invisibly, the `ggplot` object.
