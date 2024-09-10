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
#' @import pbapply
#' @import parallel
#' @import utils
#' @importFrom proxy pr_DB dist
#' @importFrom methods is
#' @importFrom rlang check_installed
#' @importFrom sf st_as_sf st_intersection st_area
#' @importFrom vegan spantree
#' @importFrom utils combn
#' @importFrom spatstat.explore density.ppp
#' @importFrom spatstat.geom as.ppp
#' @importFrom raster extend raster values extent coordinates area resample
#' @importFrom viridis viridis inferno
#' @importFrom stats as.formula sd terms
#' @importFrom grDevices adjustcolor colorRampPalette
#' @importFrom graphics image legend par points
#' @author Marcelo Araya-Salas
#'
#'   Maintainer: Marcelo Araya-Salas (\email{marcelo.araya@@ucr.ac.cr})
#'
#' @docType package
#' @details License: GPL (>= 2)
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

