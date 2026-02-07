# Estimates the size of phenotypic spaces

`space_size`

## Usage

``` r
space_size(
  formula,
  data,
  cores = 1,
  method = "mcp",
  pb = TRUE,
  outliers = 0.95,
  ...
)
```

## Arguments

- formula:

  an object of class "formula" (or one that can be coerced to that
  class).Must follow the form `group ~ dim1 + dim2` where dim1 and dim2
  are the dimensions of the phenotype space and `group` refers to the
  group labels.

- data:

  Data frame containing columns for the dimensions of the phenotypic
  space (numeric) and a categorical or factor column with group labels.

- cores:

  Numeric vector of length 1. Controls whether parallel computing is
  applied by specifying the number of cores to be used. Default is 1
  (i.e. no parallel computing).

- method:

  Character vector of length 1. Controls the method to be used for
  quantifying space size. Three metrics are available:

  - `mcp`: minimum convex polygon area using the function
    [`mcp`](https://rdrr.io/pkg/adehabitatHR/man/mcp.html). The minimum
    sample size (per group) must be 2 observations. Only works on
    2-dimensional spaces.

  - `density`: kernel density area using the function
    [`kernelUD`](https://rdrr.io/pkg/adehabitatHR/man/kernelUD.html).
    The minimum sample size (per group) must be 6 observations. Only
    works on 2-dimensional spaces.

  - `mst`: minimum spanning tree using the function
    [`spantree`](https://vegandevs.github.io/vegan/reference/spantree.html).
    The minimum sample size (per group) must be 2 observations. This
    method is expected to be more robust to the influence of outliers.
    Note that mst is not a actually measuring area but distance between
    observations. However, it still help to quantify the size of the
    sub-region in trait space. Any number of dimensions can be used with
    this method.

  - `ellipse`: Calculate the size of an sub-region assuming an
    elliptical shape. The axes of the ellipse are estimated from the
    covariance matrix of the data points in the sub-region. Estimated
    with the function
    [`niche.size`](https://rdrr.io/pkg/nicheROVER/man/niche.size.html)
    from the package 'nicheROVER'. The minimum sample size (per group)
    must be 1 observation. Any number of dimensions can be used with
    this method.

- pb:

  Logical argument to control if progress bar is shown. Default is
  `TRUE`.

- outliers:

  Numeric vector of length 1. A value between 0 and 1 controlling the
  proportion of outlier observations to be excluded. Outliers are
  determined as those farthest away from the sub-space centroid.

- ...:

  Additional arguments to be passed to
  [`kernelUD`](https://rdrr.io/pkg/adehabitatHR/man/kernelUD.html) for
  kernel density estimation (when `method = 'density'`.

## Value

A data frame containing the phenotypic space size for each group.

## Details

The function quantifies the size of the phenotypic sub-spaces.

## References

Araya-Salas, M, K. Odom. & A. Rico-Guevara. 2022, PhenotypeSpace: an R
package to quantify and compare phenotypic trait spaces R package
version 0.1.0.

## See also

[`rarefact_space_size`](https://marce10.github.io/PhenotypeSpace/reference/rarefact_space_size.md),
[`space_size_difference`](https://marce10.github.io/PhenotypeSpace/reference/space_size_difference.md),
[`rarefact_space_size_difference`](https://marce10.github.io/PhenotypeSpace/reference/rarefact_space_size_difference.md)

## Author

Marcelo Araya-Salas <marcelo.araya@ucr.ac.cr>)

## Examples

``` r
{
# load data
data("example_space")

# plot data
xs <- tapply(example_space$dimension_1, example_space$group, mean)
ys <- tapply(example_space$dimension_2, example_space$group, mean)
plot(example_space[, c("dimension_1", "dimension_2")], 
col = example_space$color, pch = 20, cex = 1.8)
text(xs, ys, labels = names(xs), cex = 2.5)

# MCP spaces
space_size(
 formula = group ~ dimension_1 + dimension_2,
 data = example_space,
 method = "mcp")

# MST 
space_size(
 formula = group ~ dimension_1 + dimension_2,
 data = example_space,
 method = "mst")
}

#>   group   n      size
#> 1    G1 150  6.062202
#> 2    G2 200  9.179167
#> 3    G3 500 34.469612
#> 4    G4 300 16.668042
#> 5    G5 400 26.412904
```
