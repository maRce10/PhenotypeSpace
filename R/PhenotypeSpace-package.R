#' PhenotypeSpace: Quantifying phenotypic spaces
#'
#' Facilitates quantifying phenotypic trait spaces and comparing spaces across groups. 
#'
#'
#' The package offers functions for:
#'   \itemize{
#'   \item Estimate absolute values of phenotypic trait spaces
#'   \item Compare phenotypic trait spaces across groups
#'   \item Modify output formats to facilitate statistical analysis
#'   }
#'
#' All functions allow the parallelization of tasks, which distributes the tasks among several processors to improve computational efficiency. 
#' 
#' @import warbleR
#' @import pbapply
#' @import parallel
#' @import utils
#' @import stats
#' @import proxy
#' @importFrom vegan spantree
#' @importFrom utils combn
#' @importFrom spatstat.core density.ppp
#' @importFrom spatstat.geom as.ppp
#' @importFrom raster extend raster values extent coordinates area resample
#' @importFrom viridis viridis inferno
#' @importFrom stats as.formula dist sd
#' @author Marcelo Araya-Salas
#'
#'   Maintainer: Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#'
#' @docType package
#' @name PhenotypeSpace
#' @details License: GPL (>= 2)
NULL
#> NULL
#'
