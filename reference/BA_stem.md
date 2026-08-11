# Calculate the basal area of a tree given the diameter at breast height (dbh)

This function calculates the basal area of a tree given the diameter at
breast height (dbh).

## Usage

``` r
BA_stem(dbh)
```

## Arguments

- dbh:

  torch.Tensor The diameter at breast height of the tree.

## Value

torch.Tensor The basal area of the tree.

## Examples

``` r
dbh = torch::torch_tensor(50)
basal_area = BA_stem(dbh)
print(basal_area)
#> torch_tensor
#>  0.1963
#> [ CPUFloatType{1} ]
```
