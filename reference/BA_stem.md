# Calculate the BA_stemsal area of a tree given the diameter at breast height (dbh).

This function calculates the BA_stemsal area of a tree given the
diameter at breast height (dbh).

## Usage

``` r
BA_stem(dbh)
```

## Arguments

- dbh:

  torch.Tensor The diameter at breast height of the tree.

## Value

torch.Tensor The BA_stemsal area of the tree.

## Examples

``` r
if (FALSE) { # \dontrun{
dbh = torch::torch_tensor(50)
BA_stemsal_area = BA_stem(dbh)
print(BA_stemsal_area)
} # }
```
