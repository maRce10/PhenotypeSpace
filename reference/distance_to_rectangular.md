# Convert pairwise distance matrices to rectangular matrices

`distance_to_rectangular` converts binary triangular matrices to
rectangular matrices using Multidimensional Scaling.

## Usage

``` r
distance_to_rectangular(
  distance.matrix,
  labels = names(distance.matrix),
  n.dimensions = 2,
  metric = TRUE,
  ...
)
```

## Arguments

- distance.matrix:

  Distance matrix (i.e. object of class 'dist'). Can be created using
  the function [`dist`](https://rdrr.io/r/stats/dist.html) or converted
  to using [`as.dist`](https://rdrr.io/r/stats/dist.html).

- labels:

  Character vector or factor containing labels to be used for
  rows/columns in the output data frame. Default is
  `names(distance.matrix)`. Must be the same length as the number of
  observations in 'distance.matrix'.

- n.dimensions:

  Integer vector of length 1 indicating the number of of dimensions to
  represent distances in a new space. Default is 2.

- metric:

  Logical argument to control if Metric (a.k.a. Classical, `TRUE`,
  default) or Non-Metric MUltidimensional Scaling (`FALSE`) is used to
  project in a new n-dimension space. Non-Metric MDS is conducted using
  the function [`isoMDS`](https://rdrr.io/pkg/MASS/man/isoMDS.html)
  while Classical MDS uses the function
  [`cmdscale`](https://rdrr.io/r/stats/cmdscale.html). So yes, it is a
  silly wrapper over those 2 functions.

- ...:

  Additional arguments to be passed to
  [`train`](https://rdrr.io/pkg/caret/man/train.html) (only used for
  machine learning models).

## Value

A data frame with the new dimensions representing the position of
observations in a new n-dimension space. If `metric = FALSE` the output
data frame is embedded in a list that also includes the stress value.

## Details

It is a silly wrapper over 2 multidimensional scaling functions
([`isoMDS`](https://rdrr.io/pkg/MASS/man/isoMDS.html) and
[`cmdscale`](https://rdrr.io/r/stats/cmdscale.html)) that simplifies the
calculation of Multidimensional Scaling and formatting of its output to
be used with other functions in the package.

## References

Araya-Salas, M, K. Odom. & A. Rico-Guevara. 2022, PhenotypeSpace: an R
package to quantify and compare phenotypic trait spaces R package
version 0.1.0.

## See also

`distance_to_rectangular`,
[`rectangular_to_triangular`](https://marce10.github.io/PhenotypeSpace/reference/rectangular_to_triangular.md)

## Author

Marcelo Araya-Salas <marcelo.araya@ucr.ac.cr>)

## Examples

``` r
{
data("example_space")

dist_example <- dist(example_space[example_space$group %in% c("G1", "G2"), 
c("dimension_1", "dimension_2")])

# convert into a 2-dimension space
rect_example <- distance_to_rectangular(distance.matrix = dist_example,
metric = TRUE)

head(rect_example)

# \donttest{
# convert into a 2-dimension space with non-metric MDS
rect_example <- distance_to_rectangular(distance.matrix = dist_example, 
metric = FALSE, maxit = 3)
 # }
}
#> Error: 'labels' should have the same length as the number of observations in 'distance.matrix' (see 'attributes(distance.matrix)$Size')
```
