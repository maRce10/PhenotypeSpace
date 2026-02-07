# Get binary triangular matrices

`binary_triangular_matrix` creates binary triangular matrices
representing categorical data in a distance matrix form

## Usage

``` r
binary_triangular_matrix(group, labels = NULL)
```

## Arguments

- group:

  Character vector or factor containing categories to be represented as
  a pairwise binary matrix. Several observations per categories (at
  least some categories) are required.

- labels:

  Character vector or factor containing labels to be used for
  rows/columns in the output matrix. Optional. Default is `NULL`.

## Value

A pairwise distance matrix that represents group membership. See
details.

## Details

The function creates binary triangular matrices representing categorical
data in a pairwise distance matrix form. Such matrices represent group
membership by assigning 0 to pairs of observations that belong to the
same category (individual, group, population) and 1 to those belonging
to different categories. Binary pairwise matrices can be useful to
evaluate association between a categorical and continuous variable
(represented as pairwise distances) using Mantel test (as in Araya-Salas
et al. 2019).

## References

Araya-Salas, M, K. Odom. & A. Rico-Guevara. 2022, PhenotypeSpace: an R
package to quantify and compare phenotypic trait spaces R package
version 0.1.0.

Araya-Salas M, G Smith-vidaurre, D Mennill, P González-Gómez, J Cahill &
T Wright. 2019. Social group signatures in hummingbird displays provide
evidence of co-occurrence of vocal and visual learning. Proceedings of
the Royal Society B. 286: 20190666.

Smouse PE, Long JC, Sokal RR. 1986 Multiple regression and correlation
extensions of the Mantel test of matrix correspondence. Syst. Zool. 35,
627–632.

## See also

[`distance_to_rectangular`](https://marce10.github.io/PhenotypeSpace/reference/distance_to_rectangular.md),
[`rectangular_to_triangular`](https://marce10.github.io/PhenotypeSpace/reference/rectangular_to_triangular.md)

## Author

Marcelo Araya-Salas <marcelo.araya@ucr.ac.cr>)

## Examples

``` r
{
# create 3 groups each one with 2 observations
groups <- paste0("G", rep(1:3, each = 2))
# create binary matrix
binary_triangular_matrix(group = groups)

# create binary matrix using labels
binary_triangular_matrix(group = groups, labels = paste(groups, 1:6, sep = "-"))
}
#>      G1-1 G1-2 G2-3 G2-4 G3-5 G3-6
#> G1-1    0    0    1    1    1    1
#> G1-2    0    0    1    1    1    1
#> G2-3    1    1    0    0    1    1
#> G2-4    1    1    0    0    1    1
#> G3-5    1    1    1    1    0    0
#> G3-6    1    1    1    1    0    0
```
