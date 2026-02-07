# Convert rectangular pairwise matrices to triangular matrices

`rectangular_to_triangular` converts rectangular pairwise matrices as
those output by many PhenotypeSpace functions into triangular pairwise
matrices.

## Usage

``` r
rectangular_to_triangular(X, distance = TRUE, symmetric = TRUE)
```

## Arguments

- X:

  Data frame containing three columns. The first two columns must
  contain group labels which will appear as rows (1 column) and column
  names (2 column) in the output triangular matrix. The third column
  (and fourth column if `symmetric = FALSE`) must have the numeric
  values to be included in the output triangular matrix.

- distance:

  Logical argument to control if the input data contains pairwise
  distances (dissimilarities) or similarities. If `TRUE` then diagonal
  values are filled with 0, otherwise they are filled with 1. Note that
  diagonal values can be set with
  [`diag`](https://rdrr.io/r/base/diag.html).

- symmetric:

  Logical argument to define if values are duplicated on both
  off-diagonal triangles (a symmetric triangular matrix,
  `symmetric = TRUE`, default) or each triangle is filled with values
  from different columns (a non-symmetric triangular matrix,
  `symmetric = FALSE`). In the latter the upper triangle is filled with
  the first column and the lower triangle with the second column. In
  this case, a fourth column with numeric values should be supplied.

## Value

A pairwise triangular matrix in which labels from the first group column
in 'X' are shown in the columns and labels from the second group are
shown in the rows. If `symmetric = FALSE` the same information is shown
below and above the diagonal.

## Details

The function converts rectangular pairwise matrices as those output by
many PhenotypeSpace functions into triangular pairwise matrices. It
takes a data frame in which each observation (row) contains the pairwise
value and related labels of the 'groups' being compared. The first two
columns must contain group labels which will appear as rows (1 column)
and column names (2 column) in the output triangular matrix. The third
column (and fourth column if `symmetric = FALSE`) must have the numeric
values to be included in the output triangular matrix.

## References

Araya-Salas, M, K. Odom. & A. Rico-Guevara. 2022, PhenotypeSpace: an R
package to quantify and compare phenotypic trait spaces R package
version 0.1.0.

## See also

[`distance_to_rectangular`](https://marce10.github.io/PhenotypeSpace/reference/distance_to_rectangular.md),
[`binary_triangular_matrix`](https://marce10.github.io/PhenotypeSpace/reference/binary_triangular_matrix.md)

## Author

Marcelo Araya-Salas <marcelo.araya@ucr.ac.cr>)

## Examples

``` r
{
# load data
data("example_space")

# get proportion of space that overlaps 
prop_overlaps <- space_similarity(
 formula = group ~ dimension_1 + dimension_2,
 data = example_space,
 method = "proportional.overlap")

# get symmetric triangular matrix
rectangular_to_triangular(prop_overlaps)

# get minimum convex polygon overlap for each group (non-symmetric)
mcp_overlaps <- space_similarity(
 formula = group ~ dimension_1 + dimension_2,
 data = example_space,
 method = "mcp.overlap")

# get a non-symmetric triangular matrix
rectangular_to_triangular(mcp_overlaps, symmetric = FALSE)
}
#>            G1        G2        G3        G4        G5
#> G1 0.00000000 0.0000000 0.6043801 0.0000000 0.0000000
#> G2 0.00000000 0.0000000 1.0000000 0.0000000 0.0000000
#> G3 0.05339738 0.1419746 0.0000000 0.0000000 0.0000000
#> G4 0.00000000 0.0000000 0.0000000 0.0000000 0.4057546
#> G5 0.00000000 0.0000000 0.0000000 0.2076138 0.0000000
```
