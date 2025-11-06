#' @title Estimates pairwise size differences of phenotypic spaces
#'
#' @description \code{space_size_difference}
#' @inheritParams template_params
#' @param method Character vector of length 1. Controls the method to be used for quantifying space size. Three metrics are available:
#' \itemize{
#'  \item \code{mcp}: minimum convex polygon area using the function  \code{\link[adehabitatHR]{mcp}}. The minimum sample size (per group) must be 2 observations. Only works on 2-dimensional spaces. 
#'  \item \code{density}: kernel density area using the function \code{\link[adehabitatHR]{kernelUD}}. The minimum sample size (per group) must be 6 observations. Only works on 2-dimensional spaces.
#'  \item \code{mst}: minimum spanning tree using the function \code{\link[vegan]{spantree}}. The minimum sample size (per group) must be 5 observations. This method is expected to be more robust to the influence of outliers. Any number of dimensions can be used with this method.
#'  \item \code{ellipse}: Calculate the size of an sub-region assuming an elliptical shape. The axes of the ellipse are estimated from the covariance matrix of the data points in the sub-region. Estimated with the function \code{\link[nicheROVER]{niche.size}} from the package 'nicheROVER'. The minimum sample size (per group) must be 1 observation. Any number of dimensions can be used with this method.
#'  }
#' @param outliers Numeric vector of length 1. A value between 0 and 1 controlling the proportion of outlier observations to be excluded. Outliers are determined as those farthest away from the sub-space centroid.
#' @param ... Additional arguments to be passed to \code{\link{space_size}} for customizing space size calculation.
#' @return A data frame containing the space size difference for each pair of groups. 
#' @export
#' @name space_size_difference
#' @details The function estimates the pairwise size difference in phenotypic space as a simple subtraction between the sizes of two spaces. As such it can be seen as an additional metric of similarity complementing those found in \code{\link{space_similarity}}.  
#' @examples {
#' # load data
#' data("example_space")
#' 
#' # MCP size (try with more iterations on your own data)
#' mcp_size <- space_size_difference(
#'  formula = group ~ dimension_1 + dimension_2,
#'  data = example_space,
#'  method = "mcp")
#' 
#' # MST size
#' mcp_size <- space_size_difference(
#'  formula = group ~ dimension_1 + dimension_2,
#'  data = example_space,
#'  method = "mst")
#' }
#' @seealso \code{\link{space_size}}, \code{\link{space_similarity}}, \code{\link{rarefact_space_size_difference}}
#' @author Marcelo Araya-Salas \email{marcelo.araya@@ucr.ac.cr})
#'
#' @references
#' Araya-Salas, M, K. Odom. & A. Rico-Guevara. 2022, PhenotypeSpace: an R package to quantify and compare phenotypic trait spaces R package version 0.1.0.

space_size_difference <-
  function(formula,
           data,
           cores = 1,
           method = "mcp",
           pb = TRUE,
           outliers = 0.95,
           ...
  ){
    
    # get term names from formula 
    form_terms <- terms(formula)
    dimensions <- attr(form_terms, "term.labels")
    group <- as.character(form_terms[[2]])
    
    group_combs <- t(utils::combn(sort(unique(data[, group])), 2))
    
    propspace_size_diff_list <- lapply(1:nrow(group_combs), function(x){
      
      W <- data[data[, group] == group_combs[x, 1], ]
      Z <- data[data[, group] == group_combs[x, 2], ]     
      
      # combined both
      both <- rbind(W, Z)
      both[, group] <- as.character(both[, group])
      both[, group] <- "both"
      
      Y <- rbind(W, Z, both)
    
    sizes <- space_size(formula = formula, data = Y, cores = 1, method = method, pb = FALSE, outliers = outliers)
    
    size.diff <- sizes$size[sizes$group == group_combs[x, 1]] - sizes$size[sizes$group == group_combs[x, 2]]
    
    out <- data.frame(group.1 = group_combs[x, 1], group.2 = group_combs[x, 2], size.difference = size.diff)
    
    return(out)
    
  })
  
  result <- do.call(rbind, propspace_size_diff_list)

  return(result)  
}
