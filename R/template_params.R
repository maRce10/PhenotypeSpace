#' @param formula an object of class "formula" (or one that can be coerced to that class).Must follow the form \code{dim1 + dim2 ~ group} where dim1 and dim2 are the dimensions of the phenotype space and \code{group} refers to the group labels.
#' @param data Data frame containing columns for the dimensions of the phenotypic space (numeric) and a categorical or factor column with group labels.
#' @param cores Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param seed Integer vector of length 1 setting the seed (see \code{\link[base]{set.seed}}). If used results should be the same on different runs, so it makes them replicable.
#' @name template_params
NULL