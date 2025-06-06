% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/space_size.R
\name{space_size}
\alias{space_size}
\title{Estimates the size of phenotypic spaces}
\usage{
space_size(
  formula,
  data,
  cores = 1,
  method = "mcp",
  pb = TRUE,
  outliers = 0.95,
  ...
)
}
\arguments{
\item{formula}{an object of class "formula" (or one that can be coerced to that class).Must follow the form \code{group ~ dim1 + dim2} where dim1 and dim2 are the dimensions of the phenotype space and \code{group} refers to the group labels.}

\item{data}{Data frame containing columns for the dimensions of the phenotypic space (numeric) and a categorical or factor column with group labels.}

\item{cores}{Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).}

\item{method}{Character vector of length 1. Controls the method to be used for quantifying space size. Three metrics are available:
\itemize{
 \item \code{mcp}: minimum convex polygon area using the function  \code{\link[adehabitatHR]{mcp}}. The minimum sample size (per group) must be 2 observations. Only works on 2-dimensional spaces.
 \item \code{density}: kernel density area using the function \code{\link[adehabitatHR]{kernelUD}}. The minimum sample size (per group) must be 6 observations. Only works on 2-dimensional spaces.
 \item \code{mst}: minimum spanning tree using the function \code{\link[vegan]{spantree}}. The minimum sample size (per group) must be 2 observations. This method is expected to be more robust to the influence of outliers. Note that mst is not a actually measuring area but distance between observations. However, it still help to quantify the size of the sub-region in trait space. Any number of dimensions can be used with this method.
 \item \code{ellipse}: Calculate the size of an sub-region assuming an elliptical shape. The axes of the ellipse are estimated from the covariance matrix of the data points in the sub-region. Estimated with the function \code{\link[nicheROVER]{niche.size}} from the package 'nicheROVER'. The minimum sample size (per group) must be 1 observation. Any number of dimensions can be used with this method.
 }}

\item{pb}{Logical argument to control if progress bar is shown. Default is \code{TRUE}.}

\item{outliers}{Numeric vector of length 1. A value between 0 and 1 controlling the proportion of outlier observations to be excluded. Outliers are determined as those farthest away from the sub-space centroid.}

\item{...}{Additional arguments to be passed to \code{\link[adehabitatHR]{kernelUD}} for kernel density estimation (when  \code{method = 'density'}.}
}
\value{
A data frame containing the phenotypic space size for each group.
}
\description{
\code{space_size}
}
\details{
The function quantifies the size of the phenotypic sub-spaces.
}
\examples{
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
}
\references{
Araya-Salas, M, & K. Odom. 2022, PhenotypeSpace: an R package to quantify and compare phenotypic trait spaces R package version 0.1.0.
}
\seealso{
\code{\link{rarefact_space_size}}, \code{\link{space_size_difference}}, \code{\link{rarefact_space_size_difference}}
}
\author{
Marcelo Araya-Salas \email{marcelo.araya@ucr.ac.cr})
}
