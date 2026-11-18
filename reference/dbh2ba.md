# Convert DBH to basal area

Computes basal area in m^2 from diameter at breast height in cm: \$\$BA
= \pi \left(\frac{\mathrm{DBH}}{200}\right)^2\$\$

## Usage

``` r
dbh2ba(dbh)
```

## Arguments

- dbh:

  numeric vector of diameters (cm).

## Value

numeric vector of basal areas (m^2).

## See also

[`BA_stem()`](https://finnverse.github.io/FINN/reference/BA_stem.md) and
[`BA_stand()`](https://finnverse.github.io/FINN/reference/BA_stand.md),
the versions for `torch` tensors

## Examples

``` r
dbh2ba(c(10, 30, 50))   # basal area of one stem of 10, 30 and 50 cm
#> [1] 0.007853982 0.070685835 0.196349541

# basal area per hectare of a 0.06 ha plot from its tree list
dbh <- c(12.9, 25.4, 31.0, 48.3)
sum(dbh2ba(dbh)) / 0.06
#> [1] 5.374035
```
