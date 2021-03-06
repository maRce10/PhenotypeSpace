% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/space_similarity.R
\name{space_similarity}
\alias{space_similarity}
\title{Pairwise similarities of phenotype spaces}
\usage{
space_similarity(X, dimensions, group, parallel = 1, type = "mcp.overlap", 
pb = TRUE, outliers = 0.95, pairwise.scale = FALSE, distance.method = "Euclidean")
}
\arguments{
\item{X}{Data frame containing columns for the dimensions of the phenotypic space (numeric) and a categorical or factor column with group labels.}

\item{dimensions}{Character vector with the names of the columns containing the dimensions of the phenotypic space. For similarity metrics using overlap 2 dimensions are required. Distance metrics can used any number of dimensions.}

\item{group}{Character vector with the name of the column (character or factor) containing group labels.}

\item{parallel}{Integer vector of length 1. Controls whether parallel computing is applied. It specifies the number of cores to be used. Default is 1 (i.e. no parallel computing).}

\item{type}{Character vector of length 1. Controls the type of (di)similarity metric to be compare the phenotypic sub-spaces of two groups at the time. Seven metrics are available which quantify as pairwise sub-space overlap ('similarity') or pairwise distance between sub-spaces ('dissimilarity'):
\itemize{
 \item \code{density.overlap}: proportion of the phenotypic sub-spaces area that overlap, taking into account the irregular densities of the sub-spaces. Two groups that share their higher density areas will be more similar than similar sub-spaces that only share their lower density areas. Two values are supplied as the proportion of the space of A that overlaps B is not necessarily the same as the proportion of B that overlaps A. Similarity metric (higher values means more similar). The minimum sample size (per group) must be 6 observations. 
 \item \code{mean.density.overlap}: similar to 'density.overlap' but the two values are merged into a single pairwise mean overlap. Similarity metric (higher values means more similar). The minimum sample size (per group) must be 6 observations.
 \item \code{mcp.overlap}: proportion of the phenotypic sub-spaces area that overlap, in which areas are calculated as the minimum convex polygon of all observations for each su-space. Two values are supplied as the proportion of the space of A that overlaps B is not necessarily the same as the proportion of B that overlaps A. Similarity metric (higher values means more similar). The minimum sample size (per group) must be 5 observations.
 \item \code{mean.mcp.overlap}: similar to 'mcp.overlap' but the two values are merged into a single pairwise mean overlap. Similarity metric (higher values means more similar). The minimum sample size (per group) must be 5 observations.
 \item \code{proportional.overlap}: proportion of the joint area of both sub-spaces that overlaps (overlapped area / total area of both groups). Sub-space areas are calculated as the minimum convex polygon. Similarity metric (higher values means more similar). The minimum sample size (per group) must be 5 observations.
 \item \code{distance}: mean euclidean pairwise distance between all observations of the compared sub-spaces. Dissimilarity metric (higher values means less similar). The minimum sample size (per group) must be 1 observation.
 \item \code{centroid.distance}: euclidean distance between the centroid of the compared sub-spaces. Dissimilarity metric (higher values means less similar). The minimum sample size (per group) must be 1 observation.
 }}

\item{pb}{Logical argument to control if progress bar is shown. Default is \code{TRUE}.}

\item{outliers}{Numeric vector of length 1. A value between 0 and 1 controlling the proportion of outlier observations to be excluded. Outliers are determined as those farthest away from the sub-space centroid.}

\item{pairwise.scale}{Logical argument to control if pairwise phenotypic spaces are scaled (i.e. z-transformed) prior to similarity estimation. If so (\code{TRUE}) similarities are decoupled from the size of the global phenotypic space. Useful to compare similarities coming from different phenotypic spaces. Default is \code{FALSE}. Not available for 'density.overlap' and 'mean.density.overlap'.}

\item{distance.method}{Character vector of length 1 indicating the method to be used for measuring distances (hence only applicable when distances are calculated). Default is 'Euclidean'. All distance and similarity measures available in \code{\link[proxy]{dist}} can be used (but note that not all of them apply to continuous data). Check available metrics by running \code{summary(proxy::pr_DB)}. If a similarity measure is used similarities are converted to distances.}
}
\value{
A data frame containing the similarity metric for each pair of groups. If the similarity metric is not symmetric (e.g. the proportional area of A that overlaps B is not necessarily the same as the area of B that overlaps A, see \code{\link{space_similarity}}) separated columns are supplied for the two comparisons.
}
\description{
\code{space_similarity} estimate pairwise similarities of phenotype spaces
}
\details{
The function quantifies pairwise similarity between phenotypic sub-spaces. Similarity is evaluated as the overlap (similarity) or distance (dissimilarity) between group.
}
\examples{
{
# load data
data("example_space")

# get proportion of space that overlaps 
prop_overlaps <- space_similarity(
X = example_space,
dimensions =  c("Dimension_1", "Dimension_2"),
group = "ID",
type = "proportional.overlap")

#' # get symmetric triangular matrix
rectangular_to_triangular(prop_overlaps)

# get minimum convex polygon overlap for each group (non-symmetric)
mcp_overlaps <- space_similarity(
X = example_space,
dimensions =  c("Dimension_1", "Dimension_2"),
group = "ID",
type = "mcp.overlap")

# convert to non-symmetric triangular matrix
rectangular_to_triangular(mcp_overlaps, symmetric = FALSE)

# check available distance measures 
summary(proxy::pr_DB)

# get eculidean distances (default)
area_dist <- space_similarity(
X = example_space,
dimensions =  c("Dimension_1", "Dimension_2"),
group = "ID",
type = "distance",
distance.method = "Euclidean")

# get Canberra distances 
area_dist <- space_similarity(
X = example_space,
dimensions =  c("Dimension_1", "Dimension_2"),
group = "ID",
type = "distance",
distance.method = "Canberra")
}
}
\references{
{
Araya-Salas, M, & K. Odom. 2022, PhenotypeSpace: an R package to quantify and compare phenotypic trait spaces R package version 0.1.0.
}
}
\seealso{
\code{\link{rarefact_space_similarity}}, \code{\link{space_size_difference}}
}
\author{
Marcelo Araya-Salas \email{marcelo.araya@ucr.ac.cr})
}
