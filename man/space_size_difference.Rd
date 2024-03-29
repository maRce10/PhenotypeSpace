% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/space_size_difference.R
\name{space_size_difference}
\alias{space_size_difference}
\title{Estimates pairwise size differences of phenotypic spaces}
\usage{
space_size_difference(
  X,
  dimensions,
  group,
  cores = 1,
  type = "mcp",
  pb = TRUE,
  outliers = 0.95,
  ...
)
}
\arguments{
\item{X}{Data frame containing columns for the dimensions of the phenotypic space (numeric) and a categorical or factor column with group labels.}

\item{dimensions}{Character vector with the names of columns containing the dimensions of the phenotypic space.}

\item{group}{Character vector with the name of the column (character or factor) containing group labels.}

\item{cores}{Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).}

\item{type}{Character vector of length 1. Controls the type of metric to be used for quantifying space size. Three metrics are available:
\itemize{
 \item \code{mcp}: minimum convex polygon area using the function  \code{\link[adehabitatHR]{mcp}}. The minimum sample size (per group) must be 2 observations.
 \item \code{density}: kernel density area using the function \code{\link[adehabitatHR]{kernelUD}}. The minimum sample size (per group) must be 6 observations.
 \item \code{mst}: minimum spanning tree using the function \code{\link[vegan]{spantree}}. The minimum sample size (per group) must be 5 observations. This method is expected to be more robust to the influence of outliers.
 }}

\item{pb}{Logical argument to control if progress bar is shown. Default is \code{TRUE}.}

\item{outliers}{Numeric vector of length 1. A value between 0 and 1 controlling the proportion of outlier observations to be excluded. Outliers are determined as those farthest away from the sub-space centroid.}

\item{...}{Additional arguments to be passed to \code{\link{space_size}} for customizing space size calculation.}
}
\value{
A data frame containing the space size difference for each pair of groups.
}
\description{
\code{space_size_difference}
}
\details{
The function estimates the pairwise size difference in phenotypic space as a simple subtraction between the sizes of two spaces. As such it can be seen as an additional metric of similarity complementing those found in \code{\link{space_similarity}}.
}
\examples{
{
# load data
data("example_space")

# MCP size (try with more iterations on your own data)
mcp_size <- space_size_difference(
X = example_space,
dimensions =  c("dimension_1", "dimension_2"),
group = "group",
type = "mcp")

# MST size
mcp_size <- space_size_difference(
X = example_space,
dimensions =  c("dimension_1", "dimension_2"),
group = "group",
type = "mst")
}
}
\references{
{
Araya-Salas, M, & K. Odom. 2022, PhenotypeSpace: an R package to quantify and compare phenotypic trait spaces R package version 0.1.0.
}
}
\seealso{
\code{\link{space_size}}, \code{\link{space_similarity}}, \code{\link{rarefact_space_size_difference}}
}
\author{
Marcelo Araya-Salas \email{marcelo.araya@ucr.ac.cr})
}
