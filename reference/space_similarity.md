# Pairwise similarities of phenotype spaces

`space_similarity` estimate pairwise similarities of phenotype spaces

## Usage

``` r
space_similarity(
  formula,
  data,
  cores = 1,
  method = "mcp.overlap",
  pb = TRUE,
  outliers = 0.95,
  pairwise.scale = FALSE,
  distance.method = "Euclidean",
  seed = NULL,
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

  Character vector of length 1. Controls the method of (di)similarity
  metric to be compare the phenotypic sub-spaces of two groups at the
  time. Seven built-in metrics are available which quantify as pairwise
  sub-space overlap ('similarity') or pairwise distance between
  bi-dimensional sub-spaces ('dissimilarity'):

  - `density.overlap`: proportion of the phenotypic sub-spaces area that
    overlap, taking into account the irregular densities of the
    sub-spaces. Two groups that share their higher density areas will be
    more similar than similar sub-spaces that only share their lower
    density areas. Two values are supplied as the proportion of the
    space of A that overlaps B is not necessarily the same as the
    proportion of B that overlaps A. Similarity metric (higher values
    means more similar). The minimum sample size (per group) must be 6
    observations.

  - `mean.density.overlap`: similar to 'density.overlap' but the two
    values are merged into a single pairwise mean overlap. Similarity
    metric (higher values means more similar). The minimum sample size
    (per group) must be 6 observations.

  - `mcp.overlap`: proportion of the phenotypic sub-spaces area that
    overlap, in which areas are calculated as the minimum convex polygon
    of all observations for each sub-space. Two values are supplied as
    the proportion of the space of A that overlaps B is not necessarily
    the same as the proportion of B that overlaps A. Similarity metric
    (higher values means more similar). The minimum sample size (per
    group) must be 5 observations.

  - `mean.mcp.overlap`: similar to 'mcp.overlap' but the two values are
    merged into a single pairwise mean overlap. Similarity metric
    (higher values means more similar). The minimum sample size (per
    group) must be 5 observations.

  - `proportional.overlap`: proportion of the joint area of both
    sub-spaces that overlaps (overlapped area / total area of both
    groups). Sub-space areas are calculated as the minimum convex
    polygon. Similarity metric (higher values means more similar). The
    minimum sample size (per group) must be 5 observations.

  - `distance`: mean euclidean pairwise distance between all
    observations of the compared sub-spaces. Dissimilarity metric
    (higher values means less similar). The minimum sample size (per
    group) must be 1 observation.

  - `centroid.distance`: euclidean distance between the centroid of the
    compared sub-spaces. Dissimilarity metric (higher values means less
    similar). The minimum sample size (per group) must be 1 observation.

  - `probability`: Bayesian probability of observations of one group
    being classified as belonging to the other group. Similarity metric
    (higher values means less similar). The minimum sample size (per
    group) must be higher the number of dimensions. Probabilities are
    calculated using the function
    [`overlap`](https://rdrr.io/pkg/nicheROVER/man/overlap.html) from
    the nicheROVER package. The following values are used internally by
    [`overlap`](https://rdrr.io/pkg/nicheROVER/man/overlap.html): nreps
    = 1000, nprob = 1000, kappa = 0, Psi = 0, nu = number of
    predictors + 1. Random draws are taken from the posterior
    distribution with Normal-Inverse-Wishart (NIW) prior using the
    function
    [`niw.post`](https://rdrr.io/pkg/nicheROVER/man/niw.post.html). Take
    a look at the nicheROVER package for further details on this method.

  In addition, machine learning classification models can also be used
  for quantify dissimilarity as a measured of how discriminable two
  groups are. These models can use more than two dimensions to represent
  phenotyypic spaces. The following classification models can be used:
  "AdaBag", "avNNet", "bam", "C5.0", "C5.0Cost", "C5.0Rules",
  "C5.0Tree", "gam", "gamLoess", "glmnet", "glmStepAIC", "kernelpls",
  "kknn", "lda", "lda2", "LogitBoost", "msaenet", "multinom", "nnet",
  "null", "ownn", "parRF", "pcaNNet", "pls", "plsRglm", "pre", "qda",
  "randomGLM", "rf", "rFerns", "rocc", "rotationForest",
  "rotationForestCp", "RRF", "RRFglobal", "sda", "simpls", "slda",
  "smda", "snn", "sparseLDA", "svmLinear2", "svmLinearWeights",
  "treebag", "widekernelpls" and "wsrf". See
  <https://topepo.github.io/caret/train-models-by-tag.html> for details
  on each of these models. Additional arguments can be pased using
  `...`. Note that some machine learning methods can significantly
  affect com

- pb:

  Logical argument to control if progress bar is shown. Default is
  `TRUE`.

- outliers:

  Numeric vector of length 1. A value between 0 and 1 controlling the
  proportion of outlier observations to be excluded. Outliers are
  determined as those farthest away from the sub-space centroid. Ignored
  when using machine learning methods.

- pairwise.scale:

  Logical argument to control if pairwise phenotypic spaces are scaled
  (i.e. z-transformed) prior to similarity estimation. If so (`TRUE`)
  similarities are decoupled from the size of the global phenotypic
  space. Useful to compare similarities coming from different phenotypic
  spaces. Default is `FALSE`. Not available for 'density.overlap',
  'mean.density.overlap' or any machine learning model.

- distance.method:

  Character vector of length 1 indicating the method to be used for
  measuring distances (hence only applicable when distances are
  calculated). Available distance measures are: "Euclidean" (default),
  "Manhattan", "supremum", "Canberra", "Wave", "divergence", "Bray",
  "Soergel", "Podani", "Chord", "Geodesic" and "Whittaker". If a
  similarity measure is used similarities are converted to distances.

- seed:

  Integer number containing the random number generator (RNG) state for
  random number generation in order to make results from the machine
  learning stochastic methods replicable.

- ...:

  Additional arguments to be passed to
  [`train`](https://rdrr.io/pkg/caret/man/train.html).

## Value

A data frame containing the similarity metric for each pair of groups.
If the similarity metric is not symmetric (e.g. the proportional area of
A that overlaps B is not necessarily the same as the area of B that
overlaps A, see `space_similarity`) separated columns are supplied for
the two comparisons.

## Details

The function quantifies pairwise similarity between phenotypic
sub-spaces. The built-in methods quantify similarity as the overlap
(similarity, or machine learning based discriminability) or distance
(dissimilarity) between group. Machine learning methods implemented in
the caret package function
[`train`](https://rdrr.io/pkg/caret/man/train.html) are available to
assess the similarity of spaces as the proportion of observations that
are incorrectly classified. In this case group overlaps are the
class-wise errors (if available) while the mean overlap is calculated as
`1- model accuracy`.

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

# get proportion of space that overlaps
prop_overlaps <- space_similarity(
 formula = group ~ dimension_1 + dimension_2,
 data = example_space,
 method = "proportional.overlap")

#' # get symmetric triangular matrix
rectangular_to_triangular(prop_overlaps)

# get minimum convex polygon overlap for each group (non-symmetric)
mcp_overlaps <- space_similarity(
 formula = group ~ dimension_1 + dimension_2,
 data = example_space,
 method = "mcp.overlap")

# convert to non-symmetric triangular matrix
rectangular_to_triangular(mcp_overlaps, symmetric = FALSE)

# check available distance measures
summary(proxy::pr_DB)

# get eculidean distances (default)
area_dist <- space_similarity(
 formula = group ~ dimension_1 + dimension_2,
 data = example_space,
 method = "distance",
 distance.method = "Euclidean")

# get Canberra distances
area_dist <- space_similarity(
 formula = group ~ dimension_1 + dimension_2,
 data = example_space,
 method = "distance",
 distance.method = "Canberra")

## using machine learning classification methods

# check if caret package and needed dependencies are available
 rlang::check_installed("caret")
 rlang::check_installed("randomForest")
 
# random forest 3 dimension data, using 5 repeats and repeated CV resampling
# extract data subset
sub_data <- example_space[example_space$group %in% c("G1", "G2", "G3"), ]

# set method parameters
ctrl <- caret::trainControl(method = "repeatedcv", repeats = 5)

# get similarities ("overlap")
space_similarity(
 formula = group ~ dimension_1 + dimension_2 + dimension_3,
 data = sub_data,
 method = "rf",
 trControl = ctrl, 
 tuneLength = 4,
 seed = 123
)

# Single C5.0 Tree using boot resampling
ctrl <- caret::trainControl(method = "boot")

space_similarity(
 formula = group ~ dimension_1 + dimension_2,
 data = sub_data,
 method = "C5.0Tree",
 trControl = ctrl,
 tuneLength =  3
)
}
#> Loading required package: lattice
#>   group.1 group.2     overlap
#> 1      G1      G2 0.009488302
#> 2      G1      G3 0.046777361
#> 3      G2      G3 0.136665576
```
