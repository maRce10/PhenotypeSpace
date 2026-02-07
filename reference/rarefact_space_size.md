# Estimates rarefacted size of phenotypic spaces

`rarefact_space_size`

## Usage

``` r
rarefact_space_size(
  formula,
  data,
  n = NULL,
  replace = FALSE,
  seed = NULL,
  cores = 1,
  pb = TRUE,
  iterations = 30,
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

- n:

  Integer vector of length 1 indicating the number of samples to be use
  for rarefaction (i.e. how many samples per group will be gather at
  each iteration). Default is the minimum sample size across groups.
  Integer vector of length 1 indicating the number of samples to be use
  for rarefaction (i.e. how many samples per group will be gather at
  each iteration). Default is the minimum sample size across groups.

- replace:

  Logical argument to control if sampling is done with replacement.
  Default is `FALSE`.

- seed:

  Integer vector of length 1 setting the seed (see
  [`set.seed`](https://rdrr.io/r/base/Random.html)). If used results
  should be the same on different runs, so it makes them replicable.

- cores:

  Numeric vector of length 1. Controls whether parallel computing is
  applied by specifying the number of cores to be used. Default is 1
  (i.e. no parallel computing).

- pb:

  Logical argument to control if progress bar is shown. Default is
  `TRUE`.

- iterations:

  Integer vector of length 1. Controls how the number of times the
  rarefaction routine is iterated. Default is 30.

- ...:

  Additional arguments to be passed to
  [`space_size`](https://marce10.github.io/PhenotypeSpace/reference/space_size.md).

## Value

A data frame containing the mean, minimum, maximum and standard
deviation of the size difference across iterations for each pair of
groups.

## Details

The function applies a rarefaction sub-sampling procedure for evaluating
pairwise space size differences (internally using
[`space_size`](https://marce10.github.io/PhenotypeSpace/reference/space_size.md)).
The size of a phenotypic space might change as a function of number of
samples. Hence, ideally, spaces should be compared between groups of
similar sample sizes. Rarefaction allows to compare groups of unbalanced
sample sizes by randomly re-sampling observations using the same number
samples across groups iteratively.

## References

Araya-Salas, M, K. Odom. & A. Rico-Guevara. 2022, PhenotypeSpace: an R
package to quantify and compare phenotypic trait spaces R package
version 0.1.0.

## See also

[`rarefact_space_similarity`](https://marce10.github.io/PhenotypeSpace/reference/rarefact_space_similarity.md),
[`space_size_difference`](https://marce10.github.io/PhenotypeSpace/reference/space_size_difference.md)

## Author

Marcelo Araya-Salas <marcelo.araya@ucr.ac.cr>)

## Examples

``` r
{
# load data
data("example_space")

# get rarefacted MCP space size 
# (try with more iterations on your own data)
rarefact_space_size(
 formula = group ~ dimension_1 + dimension_2,
 data = example_space,
 method = "mcp")

# mst rarefacted
rarefact_space_size(
 formula = group ~ dimension_1 + dimension_2,
 data = example_space,
 method = "mst")
}
#> Rarefaction sample size = 150
#> Rarefaction sample size = 150
#>   group   n mean.size  min.size  max.size   sd.size
#> 1    G1 150  6.062202  6.062202  6.062202 0.0000000
#> 2    G2 150  8.000190  7.602379  8.377236 0.1661115
#> 3    G3 150 18.474301 17.351547 19.841210 0.6339556
#> 4    G4 150 11.619925 10.471626 12.544382 0.4116350
#> 5    G5 150 15.768072 15.034279 17.010687 0.5025562
```
