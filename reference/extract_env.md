# Extract Environmental Data for a Process

This function extracts and transforms environmental data according to a
specified process object. The environmental data is formatted into an
array suitable for input into the simulation model.

## Usage

``` r
extract_env(formula, env)
```

## Arguments

- formula:

  formula object `createProcess`.

- env:

  A data frame containing environmental data, which must include columns
  "siteID" and "year".

## Value

An array of environmental data formatted according to the process
object.

## Examples

``` r
if (FALSE) { # \dontrun{
env_array <- extract_env(growth_process, env_data)
} # }
```
