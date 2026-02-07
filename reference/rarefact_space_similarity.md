# Calculates rarefacted space overlaps

`rarefact_space_similarity`

## Usage

``` r
rarefact_space_similarity(
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
  [`space_similarity`](https://marce10.github.io/PhenotypeSpace/reference/space_similarity.md)
  for customizing similarity measurements.

## Value

A data frame containing the mean, minimum, maximum and standard
deviation of the similarity metric across iterations for each pair of
groups. If the similarity metric is not symmetric (e.g. the proportional
area of A that overlaps B is not necessarily the same as the area of B
that overlaps A, see
[`space_similarity`](https://marce10.github.io/PhenotypeSpace/reference/space_similarity.md))
separated columns are supplied for the two comparisons.

## Details

The function applies a rarefaction sub-sampling procedure for evaluating
pairwise space similarity (internally using
[`space_similarity`](https://marce10.github.io/PhenotypeSpace/reference/space_similarity.md)).
The spread and shape of a phenotypic space might change as a function of
number of samples. Hence, ideally, spaces should be compared between
groups of similar sample sizes. Rarefaction allows to compare groups of
unbalanced sample sizes by randomly re-sampling observations using the
same number samples across groups iteratively.

## References

Araya-Salas, M, K. Odom. & A. Rico-Guevara. 2022, PhenotypeSpace: an R
package to quantify and compare phenotypic trait spaces R package
version 0.1.0.

## See also

[`space_similarity`](https://marce10.github.io/PhenotypeSpace/reference/space_similarity.md),
[`rectangular_to_triangular`](https://marce10.github.io/PhenotypeSpace/reference/rectangular_to_triangular.md)

## Author

Marcelo Araya-Salas <marcelo.araya@ucr.ac.cr>)

## Examples

``` r
{
# load data
data("example_space")

# get proportion of space that overlaps (try with more iterations on your own data)
prop_overlaps <- rarefact_space_similarity(
 formula = group ~ dimension_1 + dimension_2,
 data = example_space,
method = "proportional.overlap", 
iterations = 5)

# get minimum convex polygon overlap for each group (non-symmetric)
mcp_overlaps <- rarefact_space_similarity(
 formula = group ~ dimension_1 + dimension_2,
 data = example_space,
 iterations = 5)

# convert to non-symmetric triangular matrix
rectangular_to_triangular(mcp_overlaps, symmetric = FALSE)
}
#> Rarefaction sample size = 150
#> Registered S3 methods overwritten by 'adehabitatMA':
#>   method                       from
#>   print.SpatialPixelsDataFrame sp  
#>   print.SpatialPixels          sp  
#> Rarefaction sample size = 150
#>           G1        G2        G3        G4        G5
#> G1 0.0000000 0.0000000 0.3155386 0.0000000 0.0000000
#> G2 0.0000000 0.0000000 0.6089176 0.0000000 0.0000000
#> G3 0.1851194 0.5224015 0.0000000 0.0000000 0.0000000
#> G4 0.0000000 0.0000000 0.0000000 0.0000000 0.3054063
#> G5 0.0000000 0.0000000 0.0000000 0.2798535 0.0000000
```
