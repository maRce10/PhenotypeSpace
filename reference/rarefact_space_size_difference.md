# Calculates rarefacted space size differences

`rarefact_space_size_difference`

## Usage

``` r
rarefact_space_size_difference(
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
  [`space_size_difference`](https://marce10.github.io/PhenotypeSpace/reference/space_size_difference.md).

## Value

A data frame containing the mean, minimum, maximum and standard
deviation of the space size difference across iterations for each pair
of groups.

## Details

The function applies a rarefaction sub-sampling procedure for evaluating
pairwise space size differences (internally using
[`space_size_difference`](https://marce10.github.io/PhenotypeSpace/reference/space_size_difference.md)).
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

# get rarefied size difference using MCP (try with more iterations on your own data)
mcp_size_diff <- rarefact_space_size_difference(
 formula = group ~ dimension_1 + dimension_2,
 data = example_space,
 method = "mcp", 
 iterations = 5)

# convert to non-symmetric triangular matrix
rectangular_to_triangular(mcp_size_diff, symmetric = FALSE)
}
#>               G1            G2            G3            G4            G5
#> G1  0.000000e+00 -2.442456e-05 -0.0004072425 -1.130053e-04 -0.0002651856
#> G2 -2.754101e-05  0.000000e+00 -0.0003828180 -8.858071e-05 -0.0002407610
#> G3 -4.526638e-04 -4.310886e-04  0.0000000000  2.942373e-04  0.0001420570
#> G4 -1.240426e-04 -1.006927e-04  0.0002529812  0.000000e+00 -0.0001521803
#> G5 -2.765404e-04 -2.535265e-04  0.0001001474 -1.728797e-04  0.0000000000
```
