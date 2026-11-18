# Calculate the basal area of a tree given the diameter at breast height (dbh)

This function calculates the basal area of a tree given the diameter at
breast height (dbh): the area of a circle of that diameter,
`pi * (dbh / 200)^2`, in square metres for a `dbh` in centimetres. It is
the tensor counterpart of
[`dbh2ba()`](https://finnverse.github.io/FINN/reference/dbh2ba.md).

## Usage

``` r
BA_stem(dbh)
```

## Arguments

- dbh:

  torch.Tensor (or numeric vector) The diameter at breast height of the
  tree, in cm.

## Value

An object of the type of `dbh` (a torch.Tensor for a tensor) with the
basal area of one stem, in m^2.

## See also

[`BA_stand()`](https://finnverse.github.io/FINN/reference/BA_stand.md)
for the basal area per hectare of a cohort,
[`dbh2ba()`](https://finnverse.github.io/FINN/reference/dbh2ba.md)

## Examples

``` r
# one stem of 10, 30 and 50 cm: 0.008, 0.071 and 0.196 m^2
dbh = torch::torch_tensor(c(10, 30, 50))
basal_area = BA_stem(dbh)
print(basal_area)
#> torch_tensor
#>  0.0079
#>  0.0707
#>  0.1963
#> [ CPUFloatType{3} ]

# the same numbers without torch
BA_stem(c(10, 30, 50))
#> [1] 0.007853982 0.070685835 0.196349541
```
