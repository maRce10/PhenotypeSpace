
#' @title Estimates the size of phenotypic spaces
#'
#' @description \code{space_size}
#' @inheritParams template_params
#' @param method Character vector of length 1. Controls the method to be used for quantifying space size. Three metrics are available:
#' \itemize{
#'  \item \code{mcp}: minimum convex polygon area using the function  \code{\link[adehabitatHR]{mcp}}. The minimum sample size (per group) must be 2 observations. Only works on 2-dimensional spaces.
#'  \item \code{density}: kernel density area using the function \code{\link[adehabitatHR]{kernelUD}}. The minimum sample size (per group) must be 6 observations. Only works on 2-dimensional spaces.
#'  \item \code{mst}: minimum spanning tree using the function \code{\link[vegan]{spantree}}. The minimum sample size (per group) must be 2 observations. This method is expected to be more robust to the influence of outliers. Note that mst is not a actually measuring area but distance between observations. However, it still help to quantify the size of the sub-region in trait space. Any number of dimensions can be used with this method.
#'  \item \code{ellipse}: Calculate the size of an sub-region assuming an elliptical shape. The axes of the ellipse are estimated from the covariance matrix of the data points in the sub-region. Estimated with the function \code{\link[nicheROVER]{niche.size}} from the package 'nicheROVER'. The minimum sample size (per group) must be 1 observation. Any number of dimensions can be used with this method.
#'  }
#' @param outliers Numeric vector of length 1. A value between 0 and 1 controlling the proportion of outlier observations to be excluded. Outliers are determined as those farthest away from the sub-space centroid.
#' @param ... Additional arguments to be passed to \code{\link[adehabitatHR]{kernelUD}} for kernel density estimation (when  \code{method = 'density'}.
#' @return A data frame containing the phenotypic space size for each group.
#' @export   
#' @name space_size
#' @details The function quantifies the size of the phenotypic sub-spaces. 
#' @examples {
#' # load data
#' data("example_space")
#' 
#' # plot data
#' xs <- tapply(example_space$dimension_1, example_space$group, mean)
#' ys <- tapply(example_space$dimension_2, example_space$group, mean)
#' plot(example_space[, c("dimension_1", "dimension_2")], 
#' col = example_space$color, pch = 20, cex = 1.8)
#' text(xs, ys, labels = names(xs), cex = 2.5)
#' 
#' # MCP spaces
#' space_size(
#'  formula = group ~ dimension_1 + dimension_2,
#'  data = example_space,
#'  method = "mcp")
#' 
#' # MST 
#' space_size(
#'  formula = group ~ dimension_1 + dimension_2,
#'  data = example_space,
#'  method = "mst")
#' }
#' @seealso \code{\link{rarefact_space_size}}, \code{\link{space_size_difference}}, \code{\link{rarefact_space_size_difference}}
#' @author Marcelo Araya-Salas \email{marcelo.araya@@ucr.ac.cr})
#'
#' @references
#' Araya-Salas, M, & K. Odom. 2022, PhenotypeSpace: an R package to quantify and compare phenotypic trait spaces R package version 0.1.0.

space_size <-
  function(formula,
           data,
           cores = 1,
           method = "mcp",
           pb = TRUE,
           outliers = 0.95,
           ...)
  {
  
    # get term names from formula 
    form_terms <- terms(formula)
    dimensions <- attr(form_terms, "term.labels")
    group <- as.character(form_terms[[2]])
    
  if (!method %in% c("mst", "mcp", "density", "ellipse"))
    stop2("'method' must be either 'mcp', 'mst', 'density' or 'ellipse'")
  
  # split in a list of data frames
  data_list <- split(x = data, f = data[, group])
  
  # minimum sample size
  n <- min(sapply(data_list, nrow))
  
  # stop if too small sample sizes
  if (n < 6 & method == "density")
    warning("There is at least one group with less than 6 observations which is the minimum needed for kernel density area estimation (method = 'density')")
  
  if (n < 5 & method == "mcp")
    warning("There is at least one group with less than 5 observations which is the minimum needed for 'mcp' area estimation (method = 'mcp'). NA's will be returned in those cases")

  # check dimensions by method
  if ((length(dimensions) > 2 | length(dimensions) < 2) & method %in% c("density", "mcp"))
    stop2("Methods 'density' and 'mcp' only work with 2 dimensions. For spaces with other dimensionality, use method = 'mst' or 'ellipse'")
  
  if (n < 5 & method == "mcp")
    warning("There is at least one group with less than 5 observations which is the minimum needed for 'mcp' area estimation (method = 'mcp'). NA's will be returned in those cases")
  
  if (n < 2 & method == "mst")
    warning("There is at least one group with less than 2 observations which is the minimum needed for 'mst' 'area' estimation (method = 'mst'). NA's will be returned in those cases")
  
  
  # function to calculate areas
  area_fun <- function(W, mth, ol, obs) {
    
    if (obs >= 6 & mth == "density" | obs >= 5 & mth == "mcp" | obs >= 2 & mth == "mst" | obs >= 0 & mth == "ellipse"){
    ## get area
    # kernel area
    if (mth == "density")
      area <- try(adehabitatHR::kernel.area(adehabitatHR::kernelUD(W, ...), percent = ol * 100)[[1]], silent = TRUE)
    
    # mcp area
    if (mth == "mcp")
      area <- try(adehabitatHR::mcp(xy = W, percent = ol * 100)$area, silent = TRUE)
    
    # mst "area"
    if (mth == "mst")
      area <- sum(vegan::spantree(stats::dist(W))$dist)
    
    # bayesian probalistic area 
    if (mth == "ellipse")
      area <- nicheROVER::niche.size(Sigma = stats::var(W), alpha = ol)
    
    # convert to NA if error
    if (is(area, "try-error")) area <- NA
    } else
      area <- NA
    
    return(area)
  }
  
  # calculate all areas
  areas_l <- pblapply_phtpspc_int(data_list, pbar = pb, cl = cores, function(Y, dims = dimensions, meth = method, outlrs = outliers) {
    
    # subset for each individual 
    W <- Y[, dims]
    
    # add dimensions as coordinates to calculate space area
    if (meth %in% c("mcp", "density"))
    sp::coordinates(W) <- stats::as.formula(paste("~ ", dims[1], "+", dims[2]))
    
    size <- area_fun(W, meth, outlrs, obs = nrow(Y))
    
    # put in a data frame
    out_df <- data.frame(group = Y[1, group], n = nrow(Y), size = size)
    
    return( out_df)
  }) 
  
  # put in a data frame
  sizes <- do.call(rbind, areas_l)
  
  # rename rows
  row.names(sizes) <- 1:nrow(sizes)
  
  return(sizes)
}
