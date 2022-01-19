# mst miminum spanning tree package emstreeR useful for low sample sizes (minimum is 2 observations). Note that mst is not a actually measuring area but distance between observations. However, it still help to quantify the size of the sub-region in trait space.   
# each level must have at least 5 data points
# ... additional parameters to be passed to kernelUD() for kernel area estimation
#' @title Estimates the size of phenotypic spaces
#'
#' @description \code{space_size}
#' @usage space_size(X, dimensions = NULL, group, parallel = 1, type = "mcp", 
#' pb = TRUE, outliers = 0.95, ...)
#' @param X Data frame containing columns for the dimensions of the phenotypic space (numeric) and a categorical or factor column with group labels. 
#' @param dimensions Character vector with the names of columns containing the dimensions of the phenotypic space.
#' @param group Character vector with the name of the column (character or factor) containing group labels.
#' @param parallel Integer vector of length 1. Controls whether parallel computing is applied. It specifies the number of cores to be used. Default is 1 (i.e. no parallel computing).
#' @param type Character vector of length 1. Controls the type of metric to be used for quantifying space size. Three metrics are available:
#' \itemize{
#'  \item \code{mcp}: minimum convex polygon area using the function  \code{\link[adehabitatHR]{mcp}}. The minimum sample size (per group) must be 2 observations.
#'  \item \code{density}: kernel density area using the function \code{\link[adehabitatHR]{kernelUD}}. The minimum sample size (per group) must be 6 observations.
#'  \item \code{mst}: minimum spanning tree using the function \code{\link[vegan]{spantree}}. The minimum sample size (per group) must be 5 observations. This method is expected to be more robust to the influence of outliers.
#'  }
#' @param pb Logical argument to control if progress bar is shown. Default is \code{TRUE}.
#' @param outliers Numeric vector of length 1. A value between 0 and 1 controlling the proportion of outlier observations to be excluded. Outliers are determined as those farthest away from the sub-space centroid.
#' @param ... Additional arguments to be passed to \code{\link[adehabitatHR]{kernelUD}} for kernel density estimation (when  \code{type = 'density'}.
#' @return A data frame containing the phenotypic space size for each group.
#' @export   
#' @name space_size
#' @details The function quantifies the size of the phenotypic sub-spaces.
#' @examples {
#' # load data
#' data("example_space")
#' 
#' # plot data
#' xs <- tapply(example_space$Dimension_1, example_space$ID, mean)
#' ys <- tapply(example_space$Dimension_2, example_space$ID, mean)
#' plot(example_space[, c("Dimension_1", "Dimension_2")], 
#' col = example_space$color, pch = 20, cex = 1.8)
#' text(xs, ys, labels = names(xs), cex = 2.5)
#' 
#' # MCP spaces
#' space_size(
#' X = example_space,
#' dimensions =  c("Dimension_1", "Dimension_2"),
#' group = "ID",
#' type = "mcp")
#' 
#' # MST 
#' space_size(
#'   X = example_space,
#'   dimensions =  c("Dimension_1", "Dimension_2"),
#'   group = "ID",
#'   type = "mst")
#' }
#' @seealso \code{\link{rarefact_space_size}}, \code{\link{space_size_difference}}, \code{\link{rarefact_space_size_difference}}
#' @author Marcelo Araya-Salas \email{marcelo.araya@@ucr.ac.cr})
#'
#' @references {
#' Araya-Salas, M, & K. Odom. 2022, PhenotypeSpace: quantifying phenotypic trait spaces. R package version 0.1.0.
#' }
# last modification on jan-2022 (MAS)
space_size <- function(X, dimensions = NULL, group, parallel = 1, type = "mcp", pb = TRUE, outliers = 0.95, ...) {
  
  if (!type %in% c("mst", "mcp", "density"))
    stop("'type' must be either 'mcp', 'mst' or 'density'")
  
  # split in a list of vectors
  X_l <- split(x = X, f = X[, group])
  
  # minimum sample size
  n <- min(sapply(X_l, nrow))
  
  # stop if too small sample sizes
  if (n < 6 & type == "density")
    stop("There is at least one group with less than 6 observations which is the minimum needed for kernel density area estimation (type = 'density')")
  
  if (n < 5 & type == "mcp")
    stop("There is at least one group with less than 5 observations which is the minimum needed for 'mcp' area estimation (type = 'mcp')")
  
  if (n < 2 & type == "mst")
    stop("There is at least one group with less than 2 observations which is the minimum needed for 'mst' 'area' estimation (type = 'mst')")
  
  
  # function to calculate areas
  area_fun <- function(W, tp, ol) {
    
    # get area
    # kernel area
    if (tp == "density")
      area <- adehabitatHR::kernel.area(adehabitatHR::kernelUD(W, extent = 1.5, ...), percent = ol * 100)[[1]]
    
    # mcp area
    if (tp == "mcp")
      area <- adehabitatHR::mcp(xy = W, percent = ol * 100)$area
    
    
    # mst "area"
    if (tp == "mst")
      area <- sum(spantree(stats::dist(W))$dist)        
    
    return(area)
  }
  
  # calculate all areas
  areas_l <- warbleR:::pblapply_wrblr_int(X_l, pbar = pb, cl = parallel, function(Y, dims = dimensions, typ = type, outlrs = outliers) {
    
    # subset for each individual 
    W <- Y[, dims]
    
    # add dimensions as coordinates to calculate acoustic space area
    if (typ != "mst")
    sp::coordinates(W) <- stats::as.formula(paste("~ ", dims[1], "+", dims[2]))
    
    size <- area_fun(W, typ, outlrs)
    
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
